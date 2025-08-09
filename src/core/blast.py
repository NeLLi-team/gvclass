import os
from collections import defaultdict
from typing import List, Dict
import click
from Bio import SeqIO
import pyswrd
import logging

from src.utils.common import validate_file_path
from src.utils.error_handling import ErrorHandler

# Set up logging
logger = logging.getLogger(__name__)
error_handler = ErrorHandler("blast")


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
            open(output_file, "a").close()
            return

        if not ref_seqs:
            logger.warning(f"No sequences found in reference file: {ref_faa}")
            print(f"No sequences found in reference file: {ref_faa}")
            open(output_file, "a").close()
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
        search_completed = False
        try:
            for hit in pyswrd.search(
                queries,
                targets,
                scorer_name="BLOSUM62",
                score_threshold=0,  # We'll filter by e-value later
                threads=threads,
            ):
                search_completed = True

                # Apply e-value filter (1e-5)
                if hit.evalue > 1e-5:
                    continue

                try:
                    # Get query and target info with bounds checking
                    query_idx = hit.query_index
                    target_idx = hit.target_index

                    if query_idx >= len(query_seqs):
                        logger.error(
                            f"Query index {query_idx} out of bounds (max: {len(query_seqs)-1})"
                        )
                        continue
                    if target_idx >= len(ref_seqs):
                        logger.error(
                            f"Target index {target_idx} out of bounds (max: {len(ref_seqs)-1})"
                        )
                        continue

                    # Get sequence IDs
                    query_id = query_seqs[query_idx].id
                    target_id = ref_seqs[target_idx].id

                    # Get alignment details from the result object with validation
                    result = hit.result
                    if not hasattr(result, "identity"):
                        logger.warning(
                            f"Hit missing identity method for {query_id} vs {target_id}"
                        )
                        continue

                    identity_pct = result.identity() * 100.0  # Convert to percentage

                    # Validate identity percentage
                    if not (0.0 <= identity_pct <= 100.0):
                        logger.warning(
                            f"Invalid identity {identity_pct} for {query_id} vs {target_id}"
                        )
                        continue

                    # Get alignment positions (pyswrd uses 0-based indexing)
                    q_start = result.query_start + 1  # Convert to 1-based
                    q_end = result.query_end + 1
                    s_start = result.target_start + 1  # Convert to 1-based
                    s_end = result.target_end + 1

                    # Calculate alignment length
                    aln_len = result.query_end - result.query_start + 1

                    # Better gap and mismatch calculation
                    gaps = 0
                    mismatches = aln_len - int(aln_len * result.identity())

                    # Try to get more accurate stats from alignment strings if available
                    if hasattr(hit, "query_aln") and hasattr(hit, "target_aln"):
                        # Count actual gaps from alignment strings
                        query_aln = hit.query_aln
                        target_aln = hit.target_aln
                        gaps = query_aln.count("-") + target_aln.count("-")
                        # Calculate mismatches excluding gaps
                        mismatches = sum(
                            1
                            for a, b in zip(query_aln, target_aln)
                            if a != b and a != "-" and b != "-"
                        )
                    elif hasattr(result, "alignment"):
                        # Fallback to CIGAR-like notation if available
                        gaps = result.alignment.count("I") + result.alignment.count("D")

                    # Validate calculated values
                    assert (
                        0.0 <= identity_pct <= 100.0
                    ), f"Invalid identity: {identity_pct}"
                    assert aln_len > 0, f"Invalid alignment length: {aln_len}"
                    assert mismatches >= 0, f"Invalid mismatch count: {mismatches}"

                    # Track statistics for validation
                    total_hits_processed += 1
                    if identity_pct == 100.0 and gaps == 0:
                        perfect_hits += 1

                    # Check for suspicious patterns
                    if identity_pct == 100.0 and aln_len < 10:
                        suspicious_patterns.append(
                            f"Perfect identity with very short alignment ({aln_len}bp) for {query_id} vs {target_id}"
                        )
                    if aln_len > len(queries[query_idx]) * 1.5:
                        suspicious_patterns.append(
                            f"Alignment length ({aln_len}) exceeds query length for {query_id}"
                        )

                    hits_by_query[query_id].append(
                        {
                            "query_id": query_id,
                            "target_id": target_id,
                            "identity": identity_pct,
                            "aln_len": aln_len,
                            "mismatches": mismatches,
                            "gaps": gaps,
                            "q_start": q_start,
                            "q_end": q_end,
                            "s_start": s_start,
                            "s_end": s_end,
                            "evalue": hit.evalue,
                            "score": hit.score,
                        }
                    )

                except AssertionError as e:
                    logger.error(f"Validation error for {query_id} vs {target_id}: {e}")
                    continue
                except Exception as e:
                    logger.error(
                        f"Unexpected error processing hit {query_id} vs {target_id}: {e}"
                    )
                    continue

        except Exception as e:
            logger.error(f"Error during pyswrd search: {e}")
            error_handler.handle_error(e, {"step": "pyswrd_search"})
            if not search_completed:
                print(f"pyswrd search failed: {e}")
                print("Creating empty output file to allow pipeline to continue.")
                open(output_file, "a").close()
                return

        # Sort and filter top hits per query, then write to output file
        with open(output_file, "w") as out:
            total_hits_written = 0
            for query_id, query_hits in hits_by_query.items():
                # Sort this query's hits by score (descending) and e-value (ascending)
                query_hits.sort(key=lambda x: (-x["score"], x["evalue"]))

                # Take only top N hits for this query
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

    except Exception as e:
        print(f"Error running pyswrd: {str(e)}")
        print("Creating empty output file to allow pipeline to continue.")
        open(output_file, "a").close()


def parse_blastp(
    blastpout: str, max_hits_per_query: int = 10, min_identity: float = 30.0
) -> List[str]:
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
    queryids_hits_dict: Dict[str, List[tuple]] = defaultdict(list)
    seen_subjects: set = set()

    # Validation tracking
    line_count = 0
    invalid_lines = 0

    try:
        # Parse the sorted BLAST output
        with open(blastpout, "r") as infile:
            for line in infile:
                line_count += 1
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    invalid_lines += 1
                    continue

                try:
                    queryname = parts[0]
                    subjectname = parts[1]
                    identity = float(parts[2])
                    score = float(parts[11])

                    # Validate parsed values
                    if not (0.0 <= identity <= 100.0):
                        print(
                            f"Warning: Invalid identity {identity} in line {line_count}"
                        )
                        invalid_lines += 1
                        continue

                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse line {line_count}: {e}")
                    invalid_lines += 1
                    continue

                # Skip hits below identity threshold
                if identity < min_identity:
                    continue

                # Store hit with score for later sorting
                if len(queryids_hits_dict[queryname]) < max_hits_per_query:
                    queryids_hits_dict[queryname].append((subjectname, score, identity))
                    seen_subjects.add(subjectname)

        # Collect all hits, maintaining score order
        all_hits = []
        for query, hits in queryids_hits_dict.items():
            # Hits are already in score order from the file
            for subject, score, identity in hits:
                all_hits.append((subject, score, identity))

        # Sort all hits by score (descending) to ensure best hits are selected first
        all_hits.sort(key=lambda x: -x[1])

        # Return unique subject IDs in score order
        besthits = []
        seen = set()
        for subject, _, _ in all_hits:
            if subject not in seen:
                besthits.append(subject)
                seen.add(subject)

        # Report validation issues if found
        if invalid_lines > 0:
            print(
                f"⚠️  WARNING: {invalid_lines}/{line_count} lines had parsing issues in BLAST output"
            )
            if invalid_lines > line_count * 0.5:
                print("   More than 50% of lines failed - check file format!")

    except Exception as e:
        print(f"Error parsing BLAST output: {e}")
        return []

    return besthits


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
    besthits = parse_blastp(blastpout)

    seenids: List[str] = []
    outseqs: List[SeqIO.SeqRecord] = []

    for seq_record in SeqIO.parse(queryfaa, "fasta"):
        if seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)

    for seq_record in SeqIO.parse(reffaa, "fasta"):
        if seq_record.description in besthits and seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)

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
