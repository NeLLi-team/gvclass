#!/usr/bin/env python
"""
Multi-file HMM search module for the pipeline.
Processes multiple HMM files separately and combines results.
"""
import os
import tempfile
from typing import List, Tuple, Any, Iterable, TextIO
from pathlib import Path
import pyhmmer.plan7
import pyhmmer.easel

from src.utils.error_handling import error_handler
from src.core.hmm_search import DEFAULT_EVALUE_THRESHOLD


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


def _decode_name(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode()
    return value


def _decode_accession(value: Any) -> Any:
    if value and isinstance(value, bytes):
        return value.decode()
    return value or "-"


def _load_query_sequences(query_file: str) -> List[Any]:
    error_handler.log_info(f"Loading query sequences from {query_file}")
    with pyhmmer.easel.SequenceFile(query_file, digital=True) as f:
        sequences = list(f)
    error_handler.log_info(f"Loaded {len(sequences)} sequences")
    return sequences


def _load_hmms(cleaned_hmm: str) -> List[Any]:
    with pyhmmer.plan7.HMMFile(cleaned_hmm) as f:
        return list(f)


def _has_model_cutoffs(hmm: Any) -> bool:
    return bool(
        hasattr(hmm, "cutoffs")
        and hmm.cutoffs
        and (
            hmm.cutoffs.gathering is not None
            or hmm.cutoffs.trusted is not None
            or hmm.cutoffs.noise is not None
        )
    )


def _run_search_for_hmms(
    hmms: List[Any], sequences: List[Any], threads: int, sensitive_mode: bool, base_name: str
) -> Iterable[Any]:
    all_have_cutoffs = all(_has_model_cutoffs(hmm) for hmm in hmms)
    if sensitive_mode:
        error_handler.log_info(
            f"  Sensitive mode enabled for {base_name}: using E-value {DEFAULT_EVALUE_THRESHOLD}"
        )
        return pyhmmer.hmmsearch(
            hmms,
            sequences,
            cpus=threads,
            E=DEFAULT_EVALUE_THRESHOLD,
            domE=DEFAULT_EVALUE_THRESHOLD,
        )

    if all_have_cutoffs:
        error_handler.log_info(f"  Using GA/TC/NC bit-score cutoffs for {base_name}")
        return pyhmmer.hmmsearch(
            hmms,
            sequences,
            cpus=threads,
            bit_cutoffs="gathering",
        )

    error_handler.log_info(
        f"  No GA cutoffs, using E-value {DEFAULT_EVALUE_THRESHOLD}"
    )
    return pyhmmer.hmmsearch(
        hmms,
        sequences,
        cpus=threads,
        E=DEFAULT_EVALUE_THRESHOLD,
        domE=DEFAULT_EVALUE_THRESHOLD,
    )


def _process_single_hmm_file(
    hmm_file: str,
    sequences: List[Any],
    threads: int,
    sensitive_mode: bool,
) -> Tuple[str, str, List[Any], Iterable[Any]]:
    base_name = Path(hmm_file).stem
    error_handler.log_info(f"Processing {base_name}.hmm")

    cleaned_hmm = clean_hmm_file(hmm_file)
    hmms = _load_hmms(cleaned_hmm)
    error_handler.log_info(f"  Loaded {len(hmms)} models from {base_name}")
    results = _run_search_for_hmms(hmms, sequences, threads, sensitive_mode, base_name)
    return base_name, cleaned_hmm, hmms, results


def _append_results(
    all_results: List[Tuple[str, Any]], base_name: str, results: Iterable[Any]
) -> None:
    for result in results:
        all_results.append((base_name, result))


def _cleanup_cleaned_hmm(cleaned_hmm: str, hmm_file: str) -> None:
    if cleaned_hmm != hmm_file and os.path.exists(cleaned_hmm):
        os.unlink(cleaned_hmm)


def run_pyhmmer_search_multi(
    hmm_files: List[str],
    query_file: str,
    output_file: str,
    threads: int = 4,
    sensitive_mode: bool = False,
) -> None:
    """
    Run HMM search using multiple HMM files separately.

    Args:
        hmm_files: List of HMM model files
        query_file: Path to query sequences
        output_file: Path to output file
        threads: Number of threads
        sensitive_mode: If True, force E-value filtering (1e-5) and skip GA cutoffs
    """
    error_handler.log_info(
        f"Running multi-file HMM search with {len(hmm_files)} model files"
    )

    sequences = _load_query_sequences(query_file)
    all_results = []
    total_models = 0

    for hmm_file in hmm_files:
        if not os.path.exists(hmm_file):
            error_handler.log_warning(f"HMM file not found: {hmm_file}")
            continue

        try:
            base_name, cleaned_hmm, hmms, results = _process_single_hmm_file(
                hmm_file, sequences, threads, sensitive_mode
            )
            total_models += len(hmms)
            _append_results(all_results, base_name, results)
            _cleanup_cleaned_hmm(cleaned_hmm, hmm_file)
        except Exception as e:
            base_name = Path(hmm_file).stem
            error_handler.log_error(f"Error processing {base_name}: {e}")
            continue

    error_handler.log_info(f"Total models processed: {total_models}")

    # Write combined results
    write_combined_results(all_results, output_file, hmm_files, sequences)


def _write_domtbl_header(out: TextIO) -> None:
    out.write(
        "# target name\taccession\ttlen\tquery name\taccession\tqlen\t"
        "E-value\tscore\tbias\t#\tof\tc-Evalue\ti-Evalue\tscore\tbias\t"
        "from\tto\tfrom\tto\tfrom\tto\tacc\n"
    )


def _get_hmm_info(top_hits: Any, source_file: str) -> Tuple[str, str, int]:
    if hasattr(top_hits, "query"):
        hmm = top_hits.query
        hmm_name = _decode_name(hmm.name)
        hmm_accession = _decode_accession(hmm.accession)
        hmm_length = hmm.M if hasattr(hmm, "M") else 0
        return hmm_name, hmm_accession, hmm_length
    return source_file, "-", 0


def _domain_coordinates(
    domain: Any, hmm_length: int, target_length: int
) -> Tuple[int, int, int, int, int]:
    if domain.alignment:
        ali = domain.alignment
        return (
            ali.hmm_from,
            ali.hmm_to,
            ali.target_from,
            ali.target_to,
            ali.target_length,
        )
    return 1, hmm_length, domain.env_from, domain.env_to, target_length


def _write_domain_line(
    out: TextIO,
    target_name: str,
    target_accession: str,
    target_length: int,
    hmm_name: str,
    hmm_accession: str,
    hmm_length: int,
    hit: Any,
    domain: Any,
    domain_idx: int,
    hmm_from: int,
    hmm_to: int,
    target_from: int,
    target_to: int,
) -> None:
    out.write(f"{target_name}\t{target_accession}\t{target_length}\t")
    out.write(f"{hmm_name}\t{hmm_accession}\t{hmm_length}\t")
    out.write(f"{hit.evalue:.2e}\t{hit.score:.1f}\t{hit.bias:.1f}\t")
    out.write(f"{domain_idx+1}\t{len(hit.domains)}\t")
    out.write(f"{domain.c_evalue:.2e}\t{domain.i_evalue:.2e}\t")
    out.write(f"{domain.score:.1f}\t{domain.bias:.1f}\t")
    out.write(f"{hmm_from + 1}\t{hmm_to + 1}\t")
    out.write(f"{target_from + 1}\t{target_to + 1}\t")
    out.write(f"{domain.env_from + 1}\t{domain.env_to + 1}\t")
    out.write("0.90\n")


def _write_top_hits_results(
    out: TextIO, top_hits: Any, source_file: str, sequences: List[Any]
) -> int:
    hmm_name, hmm_accession, hmm_length = _get_hmm_info(top_hits, source_file)
    written_hits = 0
    for hit in top_hits:
        target_name = _decode_name(hit.name)
        target_accession = "-"
        target_length = len(sequences[0]) if sequences else 1000

        for domain_idx, domain in enumerate(hit.domains):
            hmm_from, hmm_to, target_from, target_to, target_length = _domain_coordinates(
                domain, hmm_length, target_length
            )
            _write_domain_line(
                out,
                target_name,
                target_accession,
                target_length,
                hmm_name,
                hmm_accession,
                hmm_length,
                hit,
                domain,
                domain_idx,
                hmm_from,
                hmm_to,
                target_from,
                target_to,
            )
            written_hits += 1
    return written_hits


def write_combined_results(
    all_results: List[Tuple[str, Any]],
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
    _ = hmm_files
    with open(output_file, "w") as out:
        _write_domtbl_header(out)
        total_hits = 0
        for source_file, top_hits in all_results:
            total_hits += _write_top_hits_results(out, top_hits, source_file, sequences)

    error_handler.log_info(f"Wrote {total_hits} hits to {output_file}")
