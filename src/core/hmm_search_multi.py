#!/usr/bin/env python
"""
Multi-file HMM search module for the pipeline.
Processes multiple HMM files separately and combines results.
"""
import os
import tempfile
from typing import List, Tuple
from pathlib import Path
import pyhmmer.plan7
import pyhmmer.easel

from src.utils.error_handling import error_handler


def clean_hmm_file(hmm_file: str) -> str:
    """
    Clean HMM file by removing GA/TC/NC lines that cause pyhmmer to crash.

    Args:
        hmm_file: Path to HMM file

    Returns:
        Path to cleaned HMM file (same as input if no cleaning needed)
    """
    # Check if cleaning is needed
    needs_cleaning = False
    with open(hmm_file, "r") as f:
        for line in f:
            if line.startswith(("GA ", "TC ", "NC ")):
                needs_cleaning = True
                break

    if not needs_cleaning:
        return hmm_file

    error_handler.log_info(
        f"Cleaning {os.path.basename(hmm_file)} (removing GA/TC/NC lines)"
    )

    # Create cleaned version
    with tempfile.NamedTemporaryFile(mode="w", suffix=".hmm", delete=False) as tmp:
        with open(hmm_file, "r") as inf:
            for line in inf:
                # Skip problematic lines
                if not line.startswith(("GA ", "TC ", "NC ")):
                    tmp.write(line)
        return tmp.name


def run_pyhmmer_search_multi(
    hmm_files: List[str], query_file: str, output_file: str, threads: int = 4
) -> None:
    """
    Run HMM search using multiple HMM files separately.

    Args:
        hmm_files: List of HMM model files
        query_file: Path to query sequences
        output_file: Path to output file
        threads: Number of threads
    """
    error_handler.log_info(
        f"Running multi-file HMM search with {len(hmm_files)} model files"
    )

    # Load query sequences once
    error_handler.log_info(f"Loading query sequences from {query_file}")
    with pyhmmer.easel.SequenceFile(query_file, digital=True) as f:
        sequences = list(f)
    error_handler.log_info(f"Loaded {len(sequences)} sequences")

    all_results = []
    total_models = 0

    # Process each HMM file separately
    for hmm_file in hmm_files:
        if not os.path.exists(hmm_file):
            error_handler.log_warning(f"HMM file not found: {hmm_file}")
            continue

        base_name = Path(hmm_file).stem
        error_handler.log_info(f"Processing {base_name}.hmm")

        try:
            # Clean HMM file if needed
            cleaned_hmm = clean_hmm_file(hmm_file)

            # Load HMM models
            with pyhmmer.plan7.HMMFile(cleaned_hmm) as f:
                hmms = list(f)

            total_models += len(hmms)
            error_handler.log_info(f"  Loaded {len(hmms)} models from {base_name}")

            # Run search
            results = pyhmmer.hmmsearch(
                hmms,
                sequences,
                cpus=threads,
                E=10.0,  # E-value threshold
                domE=10.0,  # Domain E-value threshold
            )

            # Store results with source file info
            for result in results:
                all_results.append((base_name, result))

            # Clean up temp file if created
            if cleaned_hmm != hmm_file and os.path.exists(cleaned_hmm):
                os.unlink(cleaned_hmm)

        except Exception as e:
            error_handler.log_error(f"Error processing {base_name}: {e}")
            continue

    error_handler.log_info(f"Total models processed: {total_models}")

    # Write combined results
    write_combined_results(all_results, output_file, hmm_files, sequences)


def write_combined_results(
    all_results: List[Tuple[str, any]],
    output_file: str,
    hmm_files: List[str],
    sequences: List,
) -> None:
    """
    Write combined results from multiple HMM searches.

    Args:
        all_results: List of (source_file, results) tuples
        output_file: Output file path
        hmm_files: List of HMM files used
        sequences: Query sequences
    """
    with open(output_file, "w") as out:
        # Write header
        out.write(
            "# target name\taccession\ttlen\tquery name\taccession\tqlen\t"
            "E-value\tscore\tbias\t#\tof\tc-Evalue\ti-Evalue\tscore\tbias\t"
            "from\tto\tfrom\tto\tfrom\tto\tacc\n"
        )

        total_hits = 0

        # Process results from each HMM file
        for source_file, top_hits in all_results:
            # Get HMM info from the results
            if hasattr(top_hits, "query"):
                hmm = top_hits.query
                hmm_name = (
                    hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                )
                hmm_accession = (
                    hmm.accession.decode()
                    if hmm.accession and isinstance(hmm.accession, bytes)
                    else (hmm.accession or "-")
                )
                hmm_length = hmm.M if hasattr(hmm, "M") else 0
            else:
                hmm_name = source_file
                hmm_accession = "-"
                hmm_length = 0

            # Process hits
            for hit in top_hits:
                target_name = (
                    hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                )
                target_accession = "-"
                target_length = len(sequences[0]) if sequences else 1000

                # Process domains
                for domain_idx, domain in enumerate(hit.domains):
                    # Get alignment info
                    if domain.alignment:
                        ali = domain.alignment
                        hmm_from = ali.hmm_from
                        hmm_to = ali.hmm_to
                        target_from = ali.target_from
                        target_to = ali.target_to
                        target_length = ali.target_length
                    else:
                        hmm_from = 1
                        hmm_to = hmm_length
                        target_from = domain.env_from
                        target_to = domain.env_to

                    # Write output line
                    out.write(f"{target_name}\t{target_accession}\t{target_length}\t")
                    out.write(f"{hmm_name}\t{hmm_accession}\t{hmm_length}\t")
                    out.write(f"{hit.evalue:.2e}\t{hit.score:.1f}\t{hit.bias:.1f}\t")
                    out.write(f"{domain_idx+1}\t{len(hit.domains)}\t")
                    out.write(f"{domain.c_evalue:.2e}\t{domain.i_evalue:.2e}\t")
                    out.write(f"{domain.score:.1f}\t{domain.bias:.1f}\t")
                    # Use 1-based coordinates
                    out.write(f"{hmm_from+1}\t{hmm_to}\t")
                    out.write(f"{target_from+1}\t{target_to}\t")
                    out.write(f"{domain.env_from+1}\t{domain.env_to}\t")
                    out.write("0.90\n")  # Default accuracy

                    total_hits += 1

    error_handler.log_info(f"Wrote {total_hits} hits to {output_file}")
