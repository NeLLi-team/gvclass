"""Protein search with pyswrd and per-protein reference selection for marker trees."""

import logging
import math
import re
import tempfile
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Protocol, TypedDict

import pyswrd
from Bio.SeqIO.FastaIO import SimpleFastaParser

from src.utils.common import validate_file_path
from src.utils.error_handling import ProcessingError

logger = logging.getLogger(__name__)

#: pyswrd hits above this E-value are dropped.
MAX_EVALUE = 1e-5
#: Hits below this percent identity are not used as tree references.
MIN_REFERENCE_IDENTITY = 30.0

# One run of insertions or one run of deletions in an opal alignment string.
_GAP_RUN = re.compile(r"I+|D+")


class FastaRecord(NamedTuple):
    """One FASTA entry; ``id`` is the first word of ``title``."""

    id: str
    title: str
    sequence: str


class _SearchResult(Protocol):
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    alignment: str

    def identity(self) -> float:
        """Return the aligned-residue identity fraction."""
        ...


class _SearchHit(Protocol):
    query_index: int
    target_index: int
    evalue: float
    score: float
    result: _SearchResult


class BlastHit(TypedDict):
    """One pyswrd hit in the column order of a BLAST m8 line."""

    query_id: str
    target_id: str
    identity: float
    aln_len: int
    mismatches: int
    gap_opens: int
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


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read every record of a FASTA file into memory, in file order."""
    with Path(path).open(encoding="utf-8") as handle:
        return [
            FastaRecord(_first_word(title), title, sequence)
            for title, sequence in SimpleFastaParser(handle)
        ]


def run_blastp(
    queryfaa: str,
    reffaa: str,
    blastpout: str,
    threads: int = 8,
    top_per_query: int = 100,
) -> None:
    """Search the proteins of ``queryfaa`` against ``reffaa`` and write an m8 file.

    Args:
        queryfaa: Query protein FASTA.
        reffaa: Reference protein FASTA.
        blastpout: Output m8 file; an empty file means no hits.
        threads: Number of pyswrd threads.
        top_per_query: Maximum number of hits to keep per query protein.
    """
    query_file = validate_file_path(queryfaa, must_exist=True)
    reference_file = validate_file_path(reffaa, must_exist=True)
    output_file = validate_file_path(blastpout, must_exist=False)
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("Threads must be a positive integer")

    search_proteins(
        read_fasta(query_file),
        read_fasta(reference_file),
        output_file,
        threads=threads,
        top_per_query=top_per_query,
    )


def search_proteins(
    queries: Sequence[FastaRecord],
    references: Sequence[FastaRecord],
    m8_path: Path,
    *,
    threads: int,
    top_per_query: int,
) -> dict[str, list[BlastHit]]:
    """Search with pyswrd, write the m8 file and return the hits of each query.

    Each query keeps at most ``top_per_query`` hits with an E-value of
    ``MAX_EVALUE`` or less, ordered by score (high first), E-value and subject
    ID. The order is fixed on purpose: pyswrd emits equal-scoring hits in thread
    order, and an unfixed order changes which subjects survive a cut, then the
    tree and the taxonomy vote.

    The m8 file is published atomically. A search without hits leaves an empty
    file.

    Raises:
        ProcessingError: If pyswrd or the hit conversion fails; no file is written.
    """
    if not queries or not references:
        logger.warning(
            "No search possible: %d query and %d reference sequences; writing %s empty",
            len(queries),
            len(references),
            m8_path,
        )
        _write_m8({}, m8_path)
        return {}

    logger.info(
        "Running pyswrd search: %d query vs %d reference sequences",
        len(queries),
        len(references),
    )
    try:
        hits_by_query, pattern_warnings = _collect_hits(
            queries, references, threads=threads, top_per_query=top_per_query
        )
    except Exception as exc:  # noqa: BLE001  # the native pyswrd search raises arbitrary types
        raise ProcessingError(
            f"pyswrd search failed for {m8_path}: {exc}",
            step="pyswrd_search",
        ) from exc

    _warn_about_odd_hits(hits_by_query, pattern_warnings)
    ranked = {
        query_id: sorted(hits, key=_hit_order)[:top_per_query]
        for query_id, hits in hits_by_query.items()
    }
    _write_m8(ranked, m8_path)
    logger.info(
        "pyswrd search wrote %d hits (top %d per query) to %s",
        sum(len(hits) for hits in ranked.values()),
        top_per_query,
        m8_path,
    )
    return ranked


def select_references(
    hits_by_query: Mapping[str, Sequence[BlastHit]],
    max_refs: int,
    min_identity: float = MIN_REFERENCE_IDENTITY,
) -> list[str]:
    """Choose up to ``max_refs`` reference IDs, giving every query protein its own.

    A marker can have several query proteins, for example the copies of two
    viruses in one bin. Pooling all hits by score would let one long, strong
    copy fill every place. Instead the proteins are taken in ID order: rank 1 of
    each protein, then rank 2 of each, and so on, skipping subjects already
    chosen, until ``max_refs`` subjects are chosen or the hits run out. Within a
    protein, hits rank by score (high first), E-value and subject ID. With one
    protein the result is that protein's ranking cut at ``max_refs``.

    Hits below ``min_identity`` percent are ignored. The result does not depend
    on the order of the input.
    """
    rankings = [
        sorted(
            (hit for hit in hits_by_query[query_id] if hit["identity"] >= min_identity),
            key=_hit_order,
        )
        for query_id in sorted(hits_by_query)
    ]
    selected: dict[str, None] = {}
    for rank in range(max((len(ranking) for ranking in rankings), default=0)):
        for ranking in rankings:
            if rank >= len(ranking) or ranking[rank]["target_id"] in selected:
                continue
            selected[ranking[rank]["target_id"]] = None
            if len(selected) == max_refs:
                return list(selected)
    return list(selected)


def iter_blastp_records(blastpout: str | Path) -> Iterator[BlastM8Record]:
    """Yield validated records, accepting a zero-byte file as no hits."""
    blast_path = Path(blastpout)
    try:
        with blast_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                yield _parse_blastp_record(line, blast_path, line_number)
    except OSError as exc:
        raise OSError(f"Failed to read BLAST m8 file {blast_path}: {exc}") from exc


def _first_word(title: str) -> str:
    words = title.split(None, 1)
    return words[0] if words else ""


def _hit_order(hit: BlastHit) -> tuple[float, float, str]:
    return (-hit["score"], hit["evalue"], hit["target_id"])


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


def _collect_hits(
    queries: Sequence[FastaRecord],
    references: Sequence[FastaRecord],
    *,
    threads: int,
    top_per_query: int,
) -> tuple[dict[str, list[BlastHit]], list[str]]:
    """Run pyswrd and group the E-value-filtered hits by query ID."""
    hits_by_query: dict[str, list[BlastHit]] = defaultdict(list)
    pattern_warnings: list[str] = []
    search = pyswrd.search(
        [record.sequence for record in queries],
        [record.sequence for record in references],
        scorer_name="BLOSUM62",
        score_threshold=0,
        # pyswrd defaults max_alignments to 10, which would cap the hits per
        # query below top_per_query before our own cut.
        max_alignments=top_per_query,
        threads=threads,
    )
    for hit in search:
        if hit.evalue > MAX_EVALUE:
            continue
        query = queries[hit.query_index]
        blast_hit, hit_warnings = _parse_search_hit(
            hit, query, references[hit.target_index]
        )
        hits_by_query[query.id].append(blast_hit)
        pattern_warnings.extend(hit_warnings)
    return hits_by_query, pattern_warnings


def _parse_search_hit(
    hit: _SearchHit,
    query: FastaRecord,
    reference: FastaRecord,
) -> tuple[BlastHit, list[str]]:
    """Convert one pyswrd hit into an m8 hit and a list of suspicious patterns.

    The columns follow BLAST tabular output: the alignment length counts gap
    columns, mismatches count substituted residues and gap openings count runs
    of inserted or deleted residues. The identity is the exception: it is
    pyswrd's own value, matches divided by aligned residue pairs, so gap columns
    do not lower it. BLAST divides by the alignment length. The contamination
    model was trained on pyswrd's value, so it stays.
    """
    result = hit.result
    identity = result.identity() * 100.0
    if not math.isfinite(identity) or not 0.0 <= identity <= 100.0:
        raise ValueError(
            f"Invalid identity {identity} for {query.id} vs {reference.id}"
        )

    alignment = result.alignment
    blast_hit: BlastHit = {
        "query_id": query.id,
        "target_id": reference.id,
        "identity": identity,
        "aln_len": len(alignment),
        "mismatches": alignment.count("X"),
        "gap_opens": len(_GAP_RUN.findall(alignment)),
        "q_start": result.query_start + 1,
        "q_end": result.query_end + 1,
        "s_start": result.target_start + 1,
        "s_end": result.target_end + 1,
        "evalue": hit.evalue,
        "score": hit.score,
    }

    warnings = []
    if _is_perfect(blast_hit) and blast_hit["aln_len"] < 10:
        warnings.append(
            f"Perfect identity with very short alignment ({blast_hit['aln_len']}bp) "
            f"for {query.id} vs {reference.id}"
        )
    if blast_hit["aln_len"] > len(query.sequence) * 1.5:
        warnings.append(
            f"Alignment length ({blast_hit['aln_len']}) exceeds query length "
            f"for {query.id}"
        )
    return blast_hit, warnings


def _is_perfect(hit: BlastHit) -> bool:
    return hit["identity"] == 100.0 and hit["gap_opens"] == 0


def _warn_about_odd_hits(
    hits_by_query: Mapping[str, Sequence[BlastHit]],
    pattern_warnings: Sequence[str],
) -> None:
    if pattern_warnings:
        extra = len(pattern_warnings) - 5
        logger.warning(
            "Suspicious patterns detected in BLAST results: %s%s",
            "; ".join(pattern_warnings[:5]),
            f"; ... and {extra} more" if extra > 0 else "",
        )

    hits = [hit for query_hits in hits_by_query.values() for hit in query_hits]
    if not hits:
        logger.warning(
            "No hits found despite having query and reference sequences. "
            "Check the E-value threshold or sequence compatibility."
        )
        return
    perfect = sum(_is_perfect(hit) for hit in hits)
    if perfect / len(hits) > 0.9:
        logger.warning(
            "%d/%d (%.1f%%) hits show perfect identity! This may indicate "
            "incorrect alignment statistics or identical sequences.",
            perfect,
            len(hits),
            100 * perfect / len(hits),
        )


def _write_m8(
    hits_by_query: Mapping[str, Sequence[BlastHit]], output_file: Path
) -> None:
    """Atomically publish BLAST m8 lines, queries in ID order, hits as given."""
    temp_path: Path | None = None
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
                for hit in hits_by_query[query_id]:
                    out.write(
                        f"{hit['query_id']}\t{hit['target_id']}\t"
                        f"{hit['identity']:.1f}\t{hit['aln_len']}\t"
                        f"{hit['mismatches']}\t{hit['gap_opens']}\t"
                        f"{hit['q_start']}\t{hit['q_end']}\t"
                        f"{hit['s_start']}\t{hit['s_end']}\t"
                        f"{hit['evalue']:.2e}\t{hit['score']:.1f}\n"
                    )
        temp_path.replace(output_file)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
