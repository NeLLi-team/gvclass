"""Protein similarity search and reference-reduction utilities."""

import logging
import math
import os
import tempfile
from collections import defaultdict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, TypedDict

import click
import pyswrd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.utils.common import validate_file_path
from src.utils.error_handling import ErrorHandler, ProcessingError

# Set up logging
logger = logging.getLogger(__name__)
error_handler = ErrorHandler("blast")


class _SearchResult(Protocol):
    query_start: int
    query_end: int
    target_start: int
    target_end: int

    def identity(self) -> float:
        """Return the aligned-residue identity fraction."""
        ...


class _SearchHit(Protocol):
    query_index: int
    target_index: int
    evalue: float
    score: float
    result: _SearchResult


class _BlastRecord(TypedDict):
    query_id: str
    target_id: str
    identity: float
    aln_len: int
    mismatches: int
    gaps: int
    q_start: int
    q_end: int
    s_start: int
    s_end: int
    evalue: float
    score: float


@dataclass(frozen=True, slots=True)
class BlastM8Record:
    """Fields consumed from one validated BLAST m8 record."""

    query_id: str
    subject_id: str
    identity: float
    score: float


def iter_blastp_records(blastpout: str | Path) -> Iterator[BlastM8Record]:
    """Yield validated records, accepting a zero-byte file as no hits."""
    blast_path = Path(blastpout)
    try:
        with blast_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                yield _parse_blastp_record(line, blast_path, line_number)
    except OSError as exc:
        raise OSError(f"Failed to read BLAST m8 file {blast_path}: {exc}") from exc


def _parse_blastp_record(
    line: str,
    blast_path: Path,
    line_number: int,
) -> BlastM8Record:
    parts = line.rstrip("\r\n").split("\t")
    if len(parts) < 12:
        raise ValueError(
            f"Malformed BLAST m8 record in {blast_path} at line {line_number}: "
            f"expected at least 12 tab-separated fields, found {len(parts)}"
        )
    query_id = parts[0]
    subject_id = parts[1]
    if not query_id or not subject_id:
        raise ValueError(
            f"Malformed BLAST m8 record in {blast_path} at line {line_number}: "
            "query and subject IDs must be nonempty"
        )
    try:
        identity = float(parts[2])
        score = float(parts[11])
    except ValueError as exc:
        raise ValueError(
            f"Malformed BLAST m8 record in {blast_path} at line {line_number}: "
            "identity and score must be numeric"
        ) from exc
    if not math.isfinite(identity) or not 0.0 <= identity <= 100.0:
        raise ValueError(
            f"Malformed BLAST m8 record in {blast_path} at line {line_number}: "
            f"identity {identity} is outside 0--100"
        )
    if not math.isfinite(score):
        raise ValueError(
            f"Malformed BLAST m8 record in {blast_path} at line {line_number}: "
            f"score {score} is not finite"
        )
    return BlastM8Record(query_id, subject_id, identity, score)


def _write_blast_output(
    hits_by_query: dict[str, list[_BlastRecord]],
    output_file: Path,
    top_per_query: int,
) -> int:
    """Atomically publish deterministically ordered BLAST m8 records."""
    temp_path: Path | None = None
    total_hits_written = 0
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_file.parent,
            prefix=f".{output_file.name}.",
            suffix=".tmp",
            delete=False,
        ) as out:
            temp_path = Path(out.name)
            for query_id in sorted(hits_by_query):
                query_hits = hits_by_query[query_id]
                query_hits.sort(
                    key=lambda hit: (
                        -hit["score"],
                        hit["evalue"],
                        hit["target_id"],
                    )
                )
                for hit_data in query_hits[:top_per_query]:
                    out.write(
                        f"{hit_data['query_id']}\t{hit_data['target_id']}\t"
                        f"{hit_data['identity']:.1f}\t{hit_data['aln_len']}\t"
                        f"{hit_data['mismatches']}\t{hit_data['gaps']}\t"
                        f"{hit_data['q_start']}\t{hit_data['q_end']}\t"
                        f"{hit_data['s_start']}\t{hit_data['s_end']}\t"
                        f"{hit_data['evalue']:.2e}\t{hit_data['score']:.1f}\n"
                    )
                    total_hits_written += 1
        temp_path.replace(output_file)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)

    return total_hits_written


def _validate_search_index(index: int, count: int, label: str) -> None:
    if not 0 <= index < count:
        raise IndexError(f"{label} index {index} out of bounds (count: {count})")


def _parse_search_hit(
    hit: _SearchHit,
    query_seqs: Sequence[SeqRecord],
    ref_seqs: Sequence[SeqRecord],
    queries: Sequence[str],
) -> tuple[_BlastRecord, bool, list[str]]:
    """Convert one pyswrd hit into a validated m8 record."""
    query_idx = hit.query_index
    target_idx = hit.target_index
    _validate_search_index(query_idx, len(query_seqs), "Query")
    _validate_search_index(target_idx, len(ref_seqs), "Target")

    query_id = query_seqs[query_idx].id
    target_id = ref_seqs[target_idx].id
    result = hit.result
    identity_method = getattr(result, "identity", None)
    if not callable(identity_method):
        raise TypeError(f"Hit has no identity method for {query_id} vs {target_id}")

    identity = identity_method()
    identity_pct = identity * 100.0
    if not math.isfinite(identity_pct) or not 0.0 <= identity_pct <= 100.0:
        raise ValueError(
            f"Invalid identity {identity_pct} for {query_id} vs {target_id}"
        )

    query_start = result.query_start + 1
    query_end = result.query_end + 1
    target_start = result.target_start + 1
    target_end = result.target_end + 1
    alignment_length = result.query_end - result.query_start + 1
    gaps = 0
    mismatches = alignment_length - int(alignment_length * identity)
    query_alignment = getattr(hit, "query_aln", None)
    target_alignment = getattr(hit, "target_aln", None)
    if query_alignment is not None and target_alignment is not None:
        gaps = query_alignment.count("-") + target_alignment.count("-")
        mismatches = sum(
            1
            for query_base, target_base in zip(
                query_alignment, target_alignment, strict=False
            )
            if query_base not in (target_base, "-") and target_base != "-"
        )
    else:
        alignment = getattr(result, "alignment", "")
        gaps = alignment.count("I") + alignment.count("D")

    if alignment_length <= 0:
        raise ValueError(
            f"Invalid alignment length {alignment_length} for {query_id} vs {target_id}"
        )
    if mismatches < 0:
        raise ValueError(
            f"Invalid mismatch count {mismatches} for {query_id} vs {target_id}"
        )

    warnings = []
    if identity_pct == 100.0 and alignment_length < 10:
        warnings.append(
            f"Perfect identity with very short alignment ({alignment_length}bp) "
            f"for {query_id} vs {target_id}"
        )
    if alignment_length > len(queries[query_idx]) * 1.5:
        warnings.append(
            f"Alignment length ({alignment_length}) exceeds query length for {query_id}"
        )

    return (
        {
            "query_id": query_id,
            "target_id": target_id,
            "identity": identity_pct,
            "aln_len": alignment_length,
            "mismatches": mismatches,
            "gaps": gaps,
            "q_start": query_start,
            "q_end": query_end,
            "s_start": target_start,
            "s_end": target_end,
            "evalue": hit.evalue,
            "score": hit.score,
        },
        identity_pct == 100.0 and gaps == 0,
        warnings,
    )


def run_blastp(
    queryfaa: str,
    refdb: str,
    blastpout: str,
    threads: int = 8,
    top_per_query: int = 100,
) -> None:
    """
    Run pyswrd for sequence similarity search (replacement for diamond blastp).

    Args:
        queryfaa (str): Path to the query FASTA file.
        refdb (str): Path to the reference database (for pyswrd, this should be the .faa file).
        blastpout (str): Path to the output file for search results.
        threads (int): Number of threads to use.
        top_per_query (int): Maximum number of hits to keep per query (default: 100).
    """
    # Validate input files
    query_file = validate_file_path(queryfaa, must_exist=True)
    output_file = validate_file_path(blastpout, must_exist=False)

    # Validate threads parameter
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("Threads must be a positive integer")

    # For pyswrd, we need the actual FASTA file, not a diamond database
    # Convert refdb path from .dmnd to .faa
    if refdb.endswith(".dmnd"):
        ref_faa = refdb.replace(".dmnd", ".faa")
        # Check in the faa directory
        if not os.path.exists(ref_faa):
            faa_dir = os.path.join(os.path.dirname(os.path.dirname(refdb)), "faa")
            ref_faa = os.path.join(faa_dir, os.path.basename(ref_faa))
    else:
        ref_faa = refdb

    ref_faa = validate_file_path(ref_faa, must_exist=True)

    logger.info(f"Running pyswrd search: {query_file} vs {ref_faa}")
    print(f"Running pyswrd search: {query_file} vs {ref_faa}")

    try:
        # Load and validate sequences
        logger.debug(f"Loading query sequences from {query_file}")
        try:
            query_seqs = list(SeqIO.parse(query_file, "fasta"))
        except Exception as e:
            error_msg = f"Failed to parse query FASTA file {query_file}: {e}"
            logger.error(error_msg)
            error_handler.handle_error(e, {"file": query_file, "type": "query"})
            raise ValueError(error_msg)

        logger.debug(f"Loading reference sequences from {ref_faa}")
        try:
            ref_seqs = list(SeqIO.parse(ref_faa, "fasta"))
        except Exception as e:
            error_msg = f"Failed to parse reference FASTA file {ref_faa}: {e}"
            logger.error(error_msg)
            error_handler.handle_error(e, {"file": ref_faa, "type": "reference"})
            raise ValueError(error_msg)

        # Validate sequences exist
        if not query_seqs:
            logger.warning(f"No sequences found in query file: {query_file}")
            print(f"No sequences found in query file: {query_file}")
            _write_blast_output({}, output_file, top_per_query)
            return

        if not ref_seqs:
            logger.warning(f"No sequences found in reference file: {ref_faa}")
            print(f"No sequences found in reference file: {ref_faa}")
            _write_blast_output({}, output_file, top_per_query)
            return

        # Validate sequences are valid
        for i, seq in enumerate(query_seqs[:5]):  # Check first 5
            if not seq.seq or len(seq.seq) == 0:
                logger.warning(f"Empty sequence found in query: {seq.id}")

        logger.info(
            f"Loaded {len(query_seqs)} query sequences and {len(ref_seqs)} reference sequences"
        )

        # Prepare sequences for pyswrd
        # Create lists of sequence strings
        queries = [str(seq.seq) for seq in query_seqs]
        targets = [str(seq.seq) for seq in ref_seqs]

        # Group hits by query for top-N filtering
        hits_by_query = defaultdict(list)

        # Track statistics for validation
        perfect_hits = 0
        total_hits_processed = 0
        suspicious_patterns = []

        # Process hits from pyswrd.search with error handling
        logger.info("Starting pyswrd search...")
        try:
            for hit in pyswrd.search(
                queries,
                targets,
                scorer_name="BLOSUM62",
                score_threshold=0,  # We'll filter by e-value later
                # pyswrd defaults max_alignments to 10, which silently caps the
                # per-query hits below the documented top_per_query (100) before
                # our own slice. Pass it explicitly so the cap is the intended
                # value, not the library default.
                max_alignments=top_per_query,
                threads=threads,
            ):
                # Apply e-value filter (1e-5)
                if hit.evalue > 1e-5:
                    continue

                hit_data, is_perfect, hit_warnings = _parse_search_hit(
                    hit, query_seqs, ref_seqs, queries
                )
                total_hits_processed += 1
                perfect_hits += int(is_perfect)
                suspicious_patterns.extend(hit_warnings)
                hits_by_query[hit_data["query_id"]].append(hit_data)

        except Exception as e:
            error_handler.handle_error(e, {"step": "pyswrd_search"})
            raise ProcessingError(
                f"pyswrd search failed for {query_file} against {ref_faa}: {e}",
                step="pyswrd_search",
                input_file=str(query_file),
            ) from e

        # Sort and filter top hits per query, then write to output file.
        #
        # Both orderings below are deterministic on purpose. pyswrd searches in
        # parallel, so ``hits_by_query`` is populated in thread-completion order
        # and equal-scoring hits arrive in an arbitrary order. Left as-is, two
        # runs of the same input emit the same hit set in a different order,
        # which changes which subjects survive the top-N cut, and from there the
        # alignment, the gene tree, the nearest neighbour, and finally the
        # per-rank taxonomy vote. This was measured: two identical runs differed
        # in the genus and species columns. Sorting by query id, and breaking
        # score/e-value ties on target id, pins all of it. Matches
        # ``_rank_unique_subjects`` in the species-tree neighbour search.
        total_hits_written = _write_blast_output(
            hits_by_query, output_file, top_per_query
        )

        # Data validation warnings
        if suspicious_patterns:
            print("⚠️  WARNING: Suspicious patterns detected in BLAST results:")
            for pattern in suspicious_patterns[:5]:  # Show first 5 warnings
                print(f"   - {pattern}")
            if len(suspicious_patterns) > 5:
                print(f"   ... and {len(suspicious_patterns) - 5} more warnings")

        if total_hits_processed > 0:
            perfect_ratio = perfect_hits / total_hits_processed
            if perfect_ratio > 0.9:
                print(
                    f"⚠️  WARNING: {perfect_hits}/{total_hits_processed} ({perfect_ratio:.1%}) hits show perfect identity!"
                )
                print(
                    "   This may indicate incorrect alignment statistics or identical sequences."
                )

        # Warn if no hits were found despite having sequences
        if total_hits_processed == 0 and len(query_seqs) > 0 and len(ref_seqs) > 0:
            print(
                "⚠️  WARNING: No hits found despite having query and reference sequences."
            )
            print("   Check e-value threshold or sequence compatibility.")

        if total_hits_written > 0:
            print(
                f"pyswrd search completed. Found {total_hits_written} hits (top {top_per_query} per query). Results written to: {output_file}"
            )
        else:
            print("pyswrd search completed. No significant hits found.")

    except ProcessingError:
        raise
    except Exception as e:
        raise ProcessingError(
            f"BLAST search failed for {query_file} against {ref_faa}: {e}",
            step="blast_search",
            input_file=str(query_file),
        ) from e


def parse_blastp(
    blastpout: str, max_hits_per_query: int = 10, min_identity: float = 30.0
) -> list[str]:
    """
    Parse the BLAST output and return the best hits sorted by score.

    The output file is already sorted by score (descending) from run_blastp.
    We now select top hits based on score and identity percentage.

    Args:
        blastpout (str): Path to the BLAST output file.
        max_hits_per_query (int): Maximum hits to keep per query sequence.
        min_identity (float): Minimum identity percentage to keep a hit.

    Returns:
        List[str]: List of best hit subject IDs sorted by score.
    """
    hits_by_query: dict[str, list[tuple[str, float, float]]] = defaultdict(list)
    for record in iter_blastp_records(blastpout):
        if record.identity < min_identity:
            continue
        if len(hits_by_query[record.query_id]) < max_hits_per_query:
            hits_by_query[record.query_id].append(
                (record.subject_id, record.score, record.identity)
            )

    all_hits = [hit for hits in hits_by_query.values() for hit in hits]
    all_hits.sort(key=lambda hit: -hit[1])

    best_hits = []
    seen_subjects = set()
    for subject, _, _ in all_hits:
        if subject not in seen_subjects:
            best_hits.append(subject)
            seen_subjects.add(subject)
    return best_hits


def reduce_ref(
    queryfaa: str,
    reffaa: str,
    reduced_reffaa_merged: str,
    refdb: str,
    blastpout: str,
    top_per_query: int = 100,
) -> None:
    """
    Reduce the reference sequences based on the best hits from diamond blastp.

    Args:
        queryfaa (str): Path to the query FASTA file.
        reffaa (str): Path to the reference FASTA file.
        reduced_reffaa_merged (str): Path to the output file for reduced reference sequences.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
        top_per_query (int): Maximum number of hits to keep per query (default: 100).
    """
    run_blastp(queryfaa, refdb, blastpout, top_per_query=top_per_query)
    besthits_set = set(parse_blastp(blastpout, max_hits_per_query=top_per_query))

    seenids: set[str] = set()
    outseqs: list[SeqRecord] = []

    for seq_record in SeqIO.parse(queryfaa, "fasta"):
        if seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.add(seq_record.id)

    for seq_record in SeqIO.parse(reffaa, "fasta"):
        if (
            seq_record.id in besthits_set or seq_record.description in besthits_set
        ) and seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.add(seq_record.id)

    SeqIO.write(outseqs, reduced_reffaa_merged, "fasta")


@click.command()
@click.option(
    "--queryfaa",
    "-q",
    type=click.Path(exists=True),
    required=True,
    help="Query sequences from hmmsearch for a single marker",
)
@click.option(
    "--reffaa",
    "-r",
    type=click.Path(exists=True),
    required=True,
    help="Reference sequences for a single marker",
)
@click.option(
    "--refdb",
    "-d",
    type=click.Path(exists=True),
    required=True,
    help="Diamond formatted database for a single marker",
)
@click.option(
    "--blastpout",
    "-b",
    type=click.Path(),
    required=True,
    help="Output file name for diamond blastp results",
)
@click.option(
    "--reduced_reffaa_merged",
    "-o",
    type=click.Path(),
    required=True,
    help="Output file name for reduced reference sequences",
)
def main(
    queryfaa: str, reffaa: str, refdb: str, blastpout: str, reduced_reffaa_merged: str
) -> None:
    """
    Main function to reduce the reference sequences based on the best hits from diamond blastp.

    Args:
        queryfaa (str): Path to the query FASTA file.
        reffaa (str): Path to the reference FASTA file.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
        reduced_reffaa_merged (str): Path to the output file for reduced reference sequences.
    """
    if os.path.getsize(queryfaa) > 0:
        reduce_ref(queryfaa, reffaa, reduced_reffaa_merged, refdb, blastpout)
    else:
        print(f"{queryfaa} is empty, no hits, omitting blast")


if __name__ == "__main__":
    main()
