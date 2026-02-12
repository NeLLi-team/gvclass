from typing import List, Tuple, Dict, Optional, Any, TextIO, Iterable
from collections import defaultdict
import pandas as pd
import pyhmmer
from pathlib import Path

from src.utils.error_handling import error_handler, ProcessingError

# Default E-value threshold when no GA/TC/NC cutoffs are defined in HMM models
DEFAULT_EVALUE_THRESHOLD = 1e-5


def load_cutoffs_simple(cutoffs_file: str) -> Dict[str, float]:
    """
    Load cutoff values for each HMM model (score only).

    Args:
        cutoffs_file: Path to cutoffs file

    Returns:
        Dictionary mapping model names to cutoff scores
    """
    cutoffs = {}
    try:
        with open(cutoffs_file, "r") as f:
            for line in f:
                parsed = _parse_cutoff_line(line)
                if parsed is None:
                    continue
                model_name, cutoff_score = parsed
                cutoffs[model_name] = cutoff_score
    except Exception as e:
        error_handler.log_error(f"Error loading cutoffs: {e}")
    return cutoffs


def _parse_cutoff_line(line: str) -> Optional[Tuple[str, float]]:
    stripped = line.strip()
    if not stripped or line.startswith("#") or line.startswith("Model"):
        return None

    parts = stripped.split("\t")
    if len(parts) < 2:
        return None

    try:
        return parts[0], float(parts[1])
    except ValueError:
        return None


def extract_cutoffs_from_hmm(
    hmm_file: str,
) -> Tuple[Dict[str, Tuple[float, float]], bool]:
    """
    Extract cutoff values from HMM model headers.

    Args:
        hmm_file: Path to HMM models file

    Returns:
        Tuple of:
        - Dictionary mapping model names to (seq_cutoff, dom_cutoff) tuples
        - Boolean indicating if any GA/TC/NC cutoffs were found
    """
    cutoffs: Dict[str, Tuple[float, float]] = {}
    has_cutoffs = False
    try:
        with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
            for hmm in hmm_handle:
                hmm_name = _decode_hmm_name(hmm.name)
                cutoff = _extract_model_cutoff(hmm)
                if cutoff is None:
                    continue
                cutoffs[hmm_name] = cutoff
                has_cutoffs = True
    except Exception as e:
        error_handler.log_error(f"Error extracting cutoffs from HMM: {e}")
    return cutoffs, has_cutoffs


def _decode_hmm_name(value: bytes | str) -> str:
    return value.decode() if isinstance(value, bytes) else value


def _extract_model_cutoff(hmm: Any) -> Optional[Tuple[float, float]]:
    if not hasattr(hmm, "cutoffs") or not hmm.cutoffs:
        return None
    if hmm.cutoffs.gathering is not None:
        return hmm.cutoffs.gathering
    if hmm.cutoffs.trusted is not None:
        return hmm.cutoffs.trusted
    if hmm.cutoffs.noise is not None:
        return hmm.cutoffs.noise
    return None


def _log_filtering_search_inputs(
    query_file: str, hmm_file: str, output_file: str
) -> None:
    error_handler.log_info("Starting HMM search with filtering:")
    error_handler.log_info(
        f"  Query file: {query_file}, exists: {Path(query_file).exists()}"
    )
    error_handler.log_info(f"  HMM file: {hmm_file}, exists: {Path(hmm_file).exists()}")
    error_handler.log_info(f"  Output file: {output_file}")


def _resolve_filter_cutoffs(
    hmm_file: str, sensitive_mode: bool
) -> Tuple[Dict[str, Tuple[float, float]], bool]:
    if sensitive_mode:
        error_handler.log_info(
            f"Sensitive mode enabled - using E-value threshold {DEFAULT_EVALUE_THRESHOLD} and skipping GA/TC/NC cutoffs"
        )
        return {}, False

    cutoffs, has_ga_cutoffs = extract_cutoffs_from_hmm(hmm_file)
    if has_ga_cutoffs:
        error_handler.log_info(
            f"Found {len(cutoffs)} GA/TC/NC cutoffs - will apply during search"
        )
    else:
        error_handler.log_info(
            f"No GA cutoffs found - using E-value threshold {DEFAULT_EVALUE_THRESHOLD}"
        )
    return cutoffs, has_ga_cutoffs


def _should_keep_filtered_hit(
    model_name: str,
    seq_score: float,
    dom_score: float,
    has_ga_cutoffs: bool,
    cutoffs: Dict[str, Tuple[float, float]],
) -> bool:
    if has_ga_cutoffs and model_name in cutoffs:
        seq_cutoff, dom_cutoff = cutoffs[model_name]
        return seq_score >= seq_cutoff and dom_score >= dom_cutoff
    return True


def _collect_filtered_search_results(
    output_file: str,
    cutoffs: Dict[str, Tuple[float, float]],
    has_ga_cutoffs: bool,
) -> Tuple[Dict[str, List[Tuple[str, float]]], List[str], int, int, int]:
    model_hits: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    filtered_lines: List[str] = []
    total_lines = 0
    hit_lines = 0
    passed_cutoffs = 0

    with open(output_file, "r") as f:
        for line in f:
            total_lines += 1
            if line.startswith("#"):
                filtered_lines.append(line)
                continue

            parts = line.strip().split()
            if len(parts) < 15:
                continue

            hit_lines += 1
            query_name = parts[0]
            model_name = parts[3]
            seq_score = float(parts[7])
            dom_score = float(parts[13])

            if not _should_keep_filtered_hit(
                model_name, seq_score, dom_score, has_ga_cutoffs, cutoffs
            ):
                continue

            filtered_lines.append(line)
            model_hits[model_name].append((query_name, seq_score))
            passed_cutoffs += 1

    return model_hits, filtered_lines, total_lines, hit_lines, passed_cutoffs


def _write_filtered_output_file(filtered_output: str, filtered_lines: List[str]) -> None:
    with open(filtered_output, "w") as f:
        f.writelines(filtered_lines)


def _write_model_hit_counts(
    counts_file: str, model_hits: Dict[str, List[Tuple[str, float]]]
) -> None:
    with open(counts_file, "w") as f:
        for model, hits in sorted(model_hits.items()):
            f.write(f"{model}\t{len(hits)}\n")


def _write_model_hit_scores(
    score_file: str, model_hits: Dict[str, List[Tuple[str, float]]]
) -> None:
    with open(score_file, "w") as f:
        total_score = sum(
            max(score for _, score in hits) for hits in model_hits.values() if hits
        )
        f.write(f"Total score: {total_score}\n")
        for model, hits in sorted(model_hits.items()):
            if hits:
                max_score = max(score for _, score in hits)
                f.write(f"{model}\t{max_score}\n")


def run_pyhmmer_search_with_filtering(
    query_file: str,
    hmm_file: str,
    output_file: str,
    filtered_output: str,
    cutoffs_file: str,  # Kept for compatibility but not used
    counts_file: str,
    score_file: str,
    threads: int = 4,
    sensitive_mode: bool = False,
) -> None:
    """Run pyhmmer, apply post-filtering, and emit filtered/count/score outputs."""
    _ = cutoffs_file
    _log_filtering_search_inputs(query_file, hmm_file, output_file)
    cutoffs, has_ga_cutoffs = _resolve_filter_cutoffs(hmm_file, sensitive_mode)

    run_pyhmmer_search(
        hmm_file,
        query_file,
        output_file,
        threads,
        sensitive_mode=sensitive_mode,
    )

    output_size = (
        Path(output_file).stat().st_size if Path(output_file).exists() else "NOT FOUND"
    )
    error_handler.log_info(f"HMM search output file size: {output_size}")

    model_hits, filtered_lines, total_lines, hit_lines, passed_cutoffs = (
        _collect_filtered_search_results(output_file, cutoffs, has_ga_cutoffs)
    )
    error_handler.log_info(
        f"Parsed {total_lines} lines, {hit_lines} hits, {passed_cutoffs} passed cutoffs"
    )

    _write_filtered_output_file(filtered_output, filtered_lines)
    error_handler.log_info(f"Wrote {len(filtered_lines)} lines to filtered output")
    error_handler.log_info(f"Found hits for {len(model_hits)} unique models")

    _write_model_hit_counts(counts_file, model_hits)
    _write_model_hit_scores(score_file, model_hits)


def _decode_name(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode()
    return value


def _decode_accession(value: Any) -> Any:
    if value and isinstance(value, bytes):
        return value.decode()
    return value or "-"


def _has_model_cutoffs(model: Any) -> bool:
    return bool(
        hasattr(model, "cutoffs")
        and model.cutoffs
        and (
            model.cutoffs.gathering is not None
            or model.cutoffs.trusted is not None
            or model.cutoffs.noise is not None
        )
    )


def _load_hmm_profiles(hmm_file: str) -> List[Any]:
    with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
        return list(hmm_handle)


def _load_query_sequences(query_file: str) -> List[Any]:
    with pyhmmer.easel.SequenceFile(query_file, digital=True) as seq_handle:
        return list(seq_handle)


def _run_hmmsearch_with_strategy(
    profiles: List[Any],
    sequences: List[Any],
    threads: int,
    evalue_threshold: float,
    sensitive_mode: bool,
) -> Iterable[Any]:
    all_have_cutoffs = all(_has_model_cutoffs(hmm) for hmm in profiles)
    if sensitive_mode:
        error_handler.log_info(
            f"Sensitive mode enabled, using E-value threshold: {evalue_threshold}"
        )
        return pyhmmer.hmmsearch(
            profiles,
            sequences,
            cpus=threads,
            E=evalue_threshold,
            domE=evalue_threshold,
        )

    if all_have_cutoffs:
        error_handler.log_info("Using GA/TC/NC bit-score cutoffs from HMM models")
        return pyhmmer.hmmsearch(
            profiles,
            sequences,
            cpus=threads,
            bit_cutoffs="gathering",
        )

    error_handler.log_info(
        f"Not all HMMs have GA cutoffs, using E-value threshold: {evalue_threshold}"
    )
    return pyhmmer.hmmsearch(
        profiles,
        sequences,
        cpus=threads,
        E=evalue_threshold,
        domE=evalue_threshold,
    )


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


def _write_domtbl_header(out_handle: TextIO) -> None:
    out_handle.write(
        "# target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value\tscore\tbias\t#\tof\tc-Evalue\ti-Evalue\tscore\tbias\tfrom\tto\tfrom\tto\tfrom\tto\tacc\n"
    )


def _write_domain_result_line(
    out_handle: TextIO,
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
    out_handle.write(f"{target_name}\t{target_accession}\t{target_length}\t")
    out_handle.write(f"{hmm_name}\t{hmm_accession}\t{hmm_length}\t")
    out_handle.write(f"{hit.evalue:.2e}\t{hit.score:.1f}\t{hit.bias:.1f}\t")
    out_handle.write(f"{domain_idx+1}\t{len(hit.domains)}\t")
    out_handle.write(f"{domain.c_evalue:.2e}\t{domain.i_evalue:.2e}\t")
    out_handle.write(f"{domain.score:.1f}\t{domain.bias:.1f}\t")

    out_handle.write(f"{hmm_from + 1}\t{hmm_to + 1}\t")
    out_handle.write(f"{target_from + 1}\t{target_to + 1}\t")
    out_handle.write(f"{domain.env_from + 1}\t{domain.env_to + 1}\t")
    out_handle.write("0.90\n")


def _write_pyhmmer_results(
    output_file: str, profiles: List[Any], results: Iterable[Any]
) -> None:
    with open(output_file, "w") as out_handle:
        _write_domtbl_header(out_handle)
        for _, (hmm, top_hits) in enumerate(zip(profiles, results)):
            hmm_name = _decode_name(hmm.name)
            hmm_accession = _decode_accession(hmm.accession)
            hmm_length = hmm.M

            for hit in top_hits:
                target_name = _decode_name(hit.name)
                target_accession = "-"
                target_length = 1000

                for domain_idx, domain in enumerate(hit.domains):
                    hmm_from, hmm_to, target_from, target_to, target_length = (
                        _domain_coordinates(domain, hmm_length, target_length)
                    )
                    _write_domain_result_line(
                        out_handle,
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


def run_pyhmmer_search(
    hmm_file: str,
    query_file: str,
    output_file: str,
    threads: int = 4,
    evalue_threshold: Optional[float] = None,
    sensitive_mode: bool = False,
) -> None:
    """
    Run HMM search using pyhmmer instead of external hmmsearch command.

    Args:
        hmm_file: Path to HMM model file
        query_file: Path to query sequence file
        output_file: Path to output file
        threads: Number of threads to use
        evalue_threshold: E-value threshold (uses DEFAULT_EVALUE_THRESHOLD if None)
        sensitive_mode: If True, force E-value filtering and skip GA cutoffs
    """
    evalue_threshold = (
        DEFAULT_EVALUE_THRESHOLD if evalue_threshold is None else evalue_threshold
    )
    try:
        error_handler.log_info(
            f"Running pyhmmer search on {query_file} with {hmm_file}"
        )
        profiles = _load_hmm_profiles(hmm_file)
        sequences = _load_query_sequences(query_file)
        results = _run_hmmsearch_with_strategy(
            profiles, sequences, threads, evalue_threshold, sensitive_mode
        )
        _write_pyhmmer_results(output_file, profiles, results)
        error_handler.log_info(
            f"pyhmmer search completed. Results written to {output_file}"
        )
    except Exception as e:
        raise ProcessingError(
            f"pyhmmer search failed: {str(e)}",
            step="pyhmmer_search",
            input_file=query_file,
        ) from e


def generate_counts_from_output(output_file: str, counts_file: str) -> None:
    """
    Generate counts file from HMM search output.

    Args:
        output_file: Path to HMM search output
        counts_file: Path to write counts file
    """
    from collections import defaultdict

    counts = defaultdict(int)

    try:
        with open(output_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    model_name = parts[3]  # Query name is the model
                    counts[model_name] += 1

        # Write counts
        with open(counts_file, "w") as f:
            for model, count in sorted(counts.items()):
                f.write(f"{model}\t{count}\n")

        error_handler.log_info(f"Generated counts file: {counts_file}")

    except Exception as e:
        error_handler.log_error(f"Error generating counts: {e}")


def get_models(modelsin: str) -> List[str]:
    """
    Generate a list of all model names from the combined HMM model file.

    Args:
        modelsin (str): Path to the combined HMM model file.

    Returns:
        List[str]: List of model names.
    """
    with open(modelsin, "r") as f:
        models = [line.split()[1].strip() for line in f if line.startswith("NAME")]
        # only get stats for general models not order level modes ("OG")
        models_general = [x for x in models if not x.startswith("OG")]

    return models_general


HitKey = Tuple[str, str]
HitRecord = Tuple[str, float, int, str]
GenomeHits = Dict[str, Dict[HitKey, HitRecord]]


def _parse_hmmout_hit(line: str) -> Tuple[str, str, str, float, int]:
    parts = line.split()
    protein_id = parts[0]
    genomeid = protein_id.split("|")[0]
    model = parts[3]
    score = float(parts[7])
    length = int(parts[2])
    return protein_id, genomeid, model, score, length


def _passes_process_cutoff(
    model: str, score: float, length: int, cutoffs: Dict[str, Tuple[float, float]]
) -> bool:
    if model not in cutoffs:
        return False
    score_cutoff, length_cutoff = cutoffs[model]
    return score >= score_cutoff and length >= length_cutoff


def _store_og_hit(
    protein_hits_order: GenomeHits,
    genomeid: str,
    key: HitKey,
    record: HitRecord,
    outfile: TextIO,
) -> None:
    if key not in protein_hits_order[genomeid]:
        protein_hits_order[genomeid][key] = record
        outfile.write(record[3])
        return

    if record[1] > protein_hits_order[genomeid][key][1]:
        protein_hits_order[genomeid][key] = record
        outfile.write(record[3])


def _find_existing_non_og_hit(
    protein_dict: Dict[HitKey, HitRecord], protein_id: str
) -> Optional[HitKey]:
    for existing_key in list(protein_dict.keys()):
        if existing_key[0] == protein_id:
            return existing_key
    return None


def _store_non_og_hit(
    protein_hits: GenomeHits, genomeid: str, key: HitKey, record: HitRecord
) -> None:
    existing_hit = _find_existing_non_og_hit(protein_hits[genomeid], key[0])
    if existing_hit is None:
        protein_hits[genomeid][key] = record
        return

    if record[1] > protein_hits[genomeid][existing_hit][1]:
        del protein_hits[genomeid][existing_hit]
        protein_hits[genomeid][key] = record


def _process_hmmout_hit_line(
    line: str,
    cutoffs: Dict[str, Tuple[float, float]],
    protein_hits: GenomeHits,
    protein_hits_order: GenomeHits,
    outfile: TextIO,
) -> None:
    protein_id, genomeid, model, score, length = _parse_hmmout_hit(line)
    if not _passes_process_cutoff(model, score, length, cutoffs):
        return

    key: HitKey = (protein_id, model)
    record: HitRecord = (model, score, length, line)
    if model.startswith("OG"):
        _store_og_hit(protein_hits_order, genomeid, key, record, outfile)
        return

    _store_non_og_hit(protein_hits, genomeid, key, record)


def _write_non_og_hits(outfile: TextIO, protein_hits: GenomeHits) -> None:
    for _, protein_dict in protein_hits.items():
        for _, (_, _, _, line) in protein_dict.items():
            outfile.write(f"{line}\n")


def _filter_hmmout_to_hits(
    hmmout: str,
    filtered_hmmout: str,
    cutoffs: Dict[str, Tuple[float, float]],
    protein_hits: GenomeHits,
    protein_hits_order: GenomeHits,
) -> None:
    with open(hmmout, "r") as infile, open(filtered_hmmout, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            _process_hmmout_hit_line(
                line, cutoffs, protein_hits, protein_hits_order, outfile
            )
        _write_non_og_hits(outfile, protein_hits)


def _count_non_og_hits(
    protein_hits: GenomeHits,
    counts: Dict[str, Dict[str, int]],
    scores: Dict[str, Dict[str, float]],
) -> int:
    total_non_og_markers = 0
    for genomeid, protein_dict in protein_hits.items():
        for (_, model), (_, score, _, _) in protein_dict.items():
            counts[genomeid][model] += 1
            scores[genomeid][model] = max(scores[genomeid][model], score)
            total_non_og_markers += 1
    return total_non_og_markers


def _count_og_hits(
    protein_hits_order: GenomeHits,
    counts: Dict[str, Dict[str, int]],
    scores: Dict[str, Dict[str, float]],
) -> Tuple[int, Dict[str, set]]:
    total_og_markers = 0
    proteins_seen: Dict[str, set] = defaultdict(set)
    for genomeid, protein_dict in protein_hits_order.items():
        for (protein_id, model), (_, score, _, _) in protein_dict.items():
            counts[genomeid][model] += 1
            scores[genomeid][model] = max(scores[genomeid][model], score)
            total_og_markers += 1
            proteins_seen[protein_id].add(model)
    return total_og_markers, proteins_seen


def _count_proteins_with_multiple_og_markers(proteins_seen: Dict[str, set]) -> int:
    proteins_with_multiple_markers = 0
    for _, model_set in proteins_seen.items():
        if len(model_set) > 1:
            proteins_with_multiple_markers += 1
    return proteins_with_multiple_markers


def _print_process_hmmout_stats(
    total_og_markers: int,
    total_non_og_markers: int,
    proteins_with_multiple_markers: int,
) -> None:
    if proteins_with_multiple_markers > 0:
        print(
            f"Note: {proteins_with_multiple_markers} proteins matched multiple OG markers (expected behavior)"
        )
    if total_og_markers == 0 and total_non_og_markers == 0:
        print("⚠️  WARNING: No markers detected in HMM search results!")
    elif total_og_markers > 0:
        print(
            f"HMM search found {total_og_markers} OG markers and {total_non_og_markers} other markers"
        )


def _initialize_missing_models(
    models: List[str],
    counts: Dict[str, Dict[str, int]],
    scores: Dict[str, Dict[str, float]],
) -> None:
    for genomeid in counts.keys():
        for model in models:
            counts[genomeid].setdefault(model, 0)
            scores[genomeid].setdefault(model, 0.0)


def process_hmmout(
    models: List[str],
    hmmout: str,
    filtered_hmmout: str,
    cutoffs: Dict[str, Tuple[float, float]],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process HMM output to generate count and score matrices including all models.

    Args:
        models (List[str]): List of model names.
        hmmout (str): Path to the HMM output file.
        filtered_hmmout (str): Path to the filtered HMM output file.
        cutoffs (Dict[str, Tuple[float, float]]): Dictionary mapping model names to (score_cutoff, length_cutoff) tuples.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Tuple containing count and score matrices.
    """
    scores: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    protein_hits: GenomeHits = defaultdict(dict)
    protein_hits_order: GenomeHits = defaultdict(dict)

    _filter_hmmout_to_hits(
        hmmout, filtered_hmmout, cutoffs, protein_hits, protein_hits_order
    )
    total_non_og_markers = _count_non_og_hits(protein_hits, counts, scores)
    total_og_markers, proteins_seen = _count_og_hits(
        protein_hits_order, counts, scores
    )
    proteins_with_multiple_markers = _count_proteins_with_multiple_og_markers(
        proteins_seen
    )
    _print_process_hmmout_stats(
        total_og_markers, total_non_og_markers, proteins_with_multiple_markers
    )
    _initialize_missing_models(models, counts, scores)

    count_df = pd.DataFrame(counts).T.fillna(0).astype(int)
    score_df = pd.DataFrame(scores).T.fillna(0)
    return count_df, score_df


def load_cutoffs(cutoff_file: str) -> Dict[str, Tuple[float, float]]:
    """
    Load cutoff values from a file.

    Args:
        cutoff_file (str): Path to the cutoff file.

    Returns:
        Dict[str, Tuple[float, float]]: Dictionary mapping model names to (score_cutoff, length_cutoff) tuples.
    """
    cutoffs = {}
    try:
        with open(cutoff_file, "r") as file:
            next(file)  # Skip header
            for line in file:
                model, score_cutoff, length_cutoff = line.strip().split("\t")
                cutoffs[model] = (float(score_cutoff), float(length_cutoff))
    except Exception as e:
        print(f"Error reading {cutoff_file}: {e}")
    return cutoffs


def main() -> None:
    """Compatibility entry point delegated to src.core.hmm_search_cli."""
    from src.core.hmm_search_cli import main as cli_main

    cli_main()


if __name__ == "__main__":
    main()
