"""Internal helpers for per-query GVClass processing."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import shutil
import tarfile
import time
import traceback

from prefect import get_run_logger

from src.core.genetic_code_optimizer import GeneticCodeOptimizer
from src.core.hmm_search import run_pyhmmer_search_with_filtering
from src.core.marker_extraction import extract_marker_hits, get_marker_database_path
from src.core.marker_processing import MarkerProcessor
from src.core.reformat import calculate_stats, reformat_sequences
from src.core.summarize_full import FullSummarizer
from src.core.tree_analysis import TreeAnalyzer
from src.pipeline.summary_writer import write_individual_summary_file
from src.utils.error_handling import ProcessingError

from Bio import SeqIO

HMM_MODEL_FILES = [
    "busco69.hmm",
    "genomad_hmm_v1.7_selection20.hmm",
    "GVOG9.hmm",
    "mirus.hmm",
    "MRYA.hmm",
    "OGv2_order.hmm",
    "UNI56.hmm",
]

@dataclass
class PreparedQueryInput:
    """Input files and metadata needed by downstream query steps."""

    protein_file: Path
    reformatted_file: Path
    best_code: Optional[int]

def process_single_query(
    query_file: Path,
    output_base: Path,
    database_path: Path,
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    sensitive_mode: bool = False,
    threads: int = 4,
) -> Dict[str, Any]:
    """Run the complete pipeline for one query."""
    logger = get_run_logger()
    query_name, query_output_dir = _initialize_query(query_file, output_base)
    logger.info(
        f"Processing query: {query_name} with {threads} threads (sensitive_mode={sensitive_mode})"
    )

    prepared_input, summary_data = _run_query_stages(
        query_file=query_file,
        query_name=query_name,
        query_output_dir=query_output_dir,
        database_path=database_path,
        genetic_codes=genetic_codes,
        tree_method=tree_method,
        mode_fast=mode_fast,
        sensitive_mode=sensitive_mode,
        threads=threads,
        logger=logger,
    )

    _create_output_dirs(query_output_dir)
    _copy_final_sequence_outputs(
        query_file,
        query_name,
        query_output_dir,
        prepared_input.reformatted_file,
        prepared_input.protein_file,
        prepared_input.best_code,
        logger,
    )
    logger.info(f"Writing summary file to {query_output_dir / f'{query_name}.summary.tab'}")
    write_individual_summary_file(query_output_dir, query_name, summary_data, logger)
    _post_process_query(query_name, query_output_dir, output_base, logger)
    return _build_query_result(query_name, summary_data, prepared_input.best_code, output_base)

def _run_query_stages(
    query_file: Path,
    query_name: str,
    query_output_dir: Path,
    database_path: Path,
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    sensitive_mode: bool,
    threads: int,
    logger,
) -> Tuple[PreparedQueryInput, Dict[str, Any]]:
    prepared_input = _prepare_query_input(
        query_file, query_output_dir, database_path, genetic_codes, threads, logger
    )
    models_out_filtered = _run_hmm_and_blast(
        query_name, query_output_dir, prepared_input.protein_file, database_path, threads, sensitive_mode, logger
    )
    marker_results = _process_markers(
        query_name, models_out_filtered, prepared_input.protein_file, query_output_dir, database_path, tree_method, mode_fast, threads, logger
    )
    tree_nn_results = _analyze_marker_trees(
        query_name, query_output_dir, database_path, marker_results, logger
    )
    summary_data = _generate_summary_data(
        query_name, query_output_dir, database_path, tree_nn_results, mode_fast, logger
    )
    return prepared_input, summary_data

def _initialize_query(query_file: Path, output_base: Path) -> Tuple[str, Path]:
    query_name = query_file.stem
    query_output_dir = output_base / query_name
    query_output_dir.mkdir(parents=True, exist_ok=True)
    return query_name, query_output_dir

def _prepare_query_input(
    query_file: Path,
    query_output_dir: Path,
    database_path: Path,
    genetic_codes: List[int],
    threads: int,
    logger,
) -> PreparedQueryInput:
    if query_file.suffix == ".fna":
        return _prepare_nucleotide_query(
            query_file, query_output_dir, database_path, genetic_codes, threads, logger
        )
    return _prepare_protein_query(query_file, query_output_dir, logger)

def _prepare_nucleotide_query(
    query_file: Path,
    query_output_dir: Path,
    database_path: Path,
    genetic_codes: List[int],
    threads: int,
    logger,
) -> PreparedQueryInput:
    query_name = query_file.stem
    logger.info(f"Reformatting nucleotide sequences: {query_name}")
    reformatted_file = query_output_dir / f"{query_name}_reformatted.fna"
    records = list(SeqIO.parse(str(query_file), "fasta"))
    _write_reformatted_sequences(records, query_file, reformatted_file)
    _write_nucleotide_stats_tsv(records, query_file, query_output_dir, query_name)
    best_code, outputs, all_code_results = _optimize_genetic_code(
        query_name,
        reformatted_file,
        query_output_dir,
        database_path,
        genetic_codes,
        threads,
        logger,
    )
    _write_genecalling_tables(query_output_dir, query_name, best_code, all_code_results)
    return PreparedQueryInput(
        protein_file=Path(outputs["faa"]),
        reformatted_file=reformatted_file,
        best_code=best_code,
    )

def _prepare_protein_query(
    query_file: Path,
    query_output_dir: Path,
    logger,
) -> PreparedQueryInput:
    query_name = query_file.stem
    logger.info(f"Using protein sequences directly: {query_name}")
    reformatted_file = query_output_dir / "query_faa" / f"{query_name}.faa"
    reformatted_file.parent.mkdir(exist_ok=True, parents=True)
    records = list(SeqIO.parse(str(query_file), "fasta"))
    reformatted_records = reformat_sequences(records, str(query_file))
    SeqIO.write(reformatted_records, str(reformatted_file), "fasta")
    _write_protein_stats(query_output_dir, query_name, reformatted_records)
    return PreparedQueryInput(
        protein_file=reformatted_file,
        reformatted_file=reformatted_file,
        best_code=None,
    )

def _write_reformatted_sequences(records, query_file: Path, reformatted_file: Path) -> None:
    reformatted_records = reformat_sequences(records, str(query_file))
    SeqIO.write(reformatted_records, str(reformatted_file), "fasta")

def _write_nucleotide_stats_tsv(
    records, query_file: Path, query_output_dir: Path, query_name: str
) -> None:
    stats_file = query_output_dir / "stats" / f"{query_name}_stats.tsv"
    stats_file.parent.mkdir(exist_ok=True)
    stats_df = calculate_stats(records, str(query_file))
    stats_df.to_csv(str(stats_file), sep="\t", index=False)

def _optimize_genetic_code(
    query_name: str,
    reformatted_file: Path,
    query_output_dir: Path,
    database_path: Path,
    genetic_codes: List[int],
    threads: int,
    logger,
) -> Tuple[int, Dict[str, Any], Dict[int, Dict[str, Any]]]:
    logger.info(f"Optimizing genetic code for: {query_name}")
    optimizer = GeneticCodeOptimizer(database_path=database_path, threads=threads)
    best_code, outputs, _best_metrics, all_code_results = optimizer.optimize_genetic_code(
        query_file=reformatted_file,
        output_dir=query_output_dir,
        genetic_codes=genetic_codes,
    )
    return best_code, outputs, all_code_results

def _write_genecalling_tables(
    query_output_dir: Path,
    query_name: str,
    best_code: int,
    all_code_results: Dict[int, Dict[str, Any]],
) -> None:
    _write_genecalling_comparison_table(query_output_dir, query_name, all_code_results)
    _write_selected_genecode_table(query_output_dir, query_name, best_code, all_code_results)

def _write_genecalling_comparison_table(
    query_output_dir: Path,
    query_name: str,
    all_code_results: Dict[int, Dict[str, Any]],
) -> None:
    comparison_file = query_output_dir / "stats" / f"{query_name}.genecalling.tab"
    with open(comparison_file, "w") as handle:
        handle.write(
            "ttable\tGenes\tBases\tCodingDensity\tcomplete_bestHits\tavg_bestHit_score\n"
        )
        for row in _iter_sorted_genecalling_rows(all_code_results):
            handle.write(
                f"{row['code']}\t{row['genes']}\t{row['bases']}\t"
                f"{row['coding_density']:.3f}\t{row['complete_hits']}\t"
                f"{row['avg_best_hit_score']:.1f}\n"
            )

def _iter_sorted_genecalling_rows(
    all_code_results: Dict[int, Dict[str, Any]]
) -> List[Dict[str, Any]]:
    rows = []
    for code, metrics in all_code_results.items():
        rows.append(
            {
                "code": code,
                "genes": metrics.get("gene_count", 0),
                "bases": metrics.get("total_bases", 0),
                "coding_density": metrics.get("coding_density", 0.0),
                "complete_hits": metrics.get("complete_hits", 0),
                "avg_best_hit_score": metrics.get("avg_score", 0.0),
            }
        )
    return sorted(rows, key=lambda item: (-item["complete_hits"], -item["avg_best_hit_score"]))

def _write_selected_genecode_table(
    query_output_dir: Path,
    query_name: str,
    best_code: int,
    all_code_results: Dict[int, Dict[str, Any]],
) -> None:
    stats_file = query_output_dir / "stats" / f"{query_name}.stats.tab"
    ttable_str = "codemeta" if best_code == 0 else str(best_code)
    with open(stats_file, "w") as handle:
        handle.write(f"ttable\t{ttable_str}\n")
        if best_code not in all_code_results:
            return
        selected = all_code_results[best_code]
        handle.write(f"genes\t{selected.get('gene_count', 0)}\n")
        handle.write(f"bases\t{selected.get('total_bases', 0)}\n")
        handle.write(f"coding_density\t{selected.get('coding_density', 0.0):.3f}\n")

def _write_protein_stats(query_output_dir: Path, query_name: str, reformatted_records) -> None:
    stats_dir = query_output_dir / "stats"
    stats_dir.mkdir(exist_ok=True)
    protein_count = len(reformatted_records)
    total_aa_length = sum(len(record.seq) for record in reformatted_records)
    estimated_genome_length = int(total_aa_length * 3 * 1.1)
    _write_protein_stats_tab(stats_dir / f"{query_name}.stats.tab", protein_count, estimated_genome_length)
    _write_protein_stats_tsv(stats_dir / f"{query_name}_stats.tsv", protein_count, estimated_genome_length)

def _write_protein_stats_tab(stats_file: Path, protein_count: int, genome_length: int) -> None:
    with open(stats_file, "w") as handle:
        handle.write("ttable\tno_fna\n")
        handle.write(f"genes\t{protein_count}\n")
        handle.write(f"bases\t{genome_length}\n")
        handle.write("coding_density\t0.0\n")

def _write_protein_stats_tsv(stats_file: Path, protein_count: int, genome_length: int) -> None:
    with open(stats_file, "w") as handle:
        handle.write("contigs\tLENbp\tGCperc\tgenecount\tCODINGperc\tttable\n")
        handle.write(f"0\t{genome_length}\t0.0\t{protein_count}\t0.0\tno_fna\n")

def _run_hmm_and_blast(
    query_name: str,
    query_output_dir: Path,
    protein_file: Path,
    database_path: Path,
    threads: int,
    sensitive_mode: bool,
    logger,
) -> Path:
    logger.info(f"Running HMM search: {query_name}")
    hmmout_dir = query_output_dir / "hmmout"
    hmmout_dir.mkdir(exist_ok=True)
    hmm_files = _resolve_hmm_files(database_path)
    models_out, models_out_filtered, models_counts, models_score = _get_hmm_output_paths(hmmout_dir)
    _run_hmm_search(
        protein_file,
        hmm_files,
        models_out,
        models_out_filtered,
        models_counts,
        models_score,
        threads,
        sensitive_mode,
    )
    from src.core.marker_extraction import parse_hmm_output

    hmm_results = parse_hmm_output(str(models_out_filtered))
    _run_blast_search(query_name, query_output_dir, database_path, protein_file, hmm_results, threads, logger)
    return models_out_filtered

def _resolve_hmm_files(database_path: Path) -> List[str]:
    model_dir = database_path / "models"
    hmm_files = [str(model_dir / model_name) for model_name in HMM_MODEL_FILES]
    existing_hmm_files = [hmm_file for hmm_file in hmm_files if Path(hmm_file).exists()]
    if existing_hmm_files:
        return existing_hmm_files
    combined_hmm = model_dir / "combined.hmm"
    if combined_hmm.exists():
        get_run_logger().info("Using fallback combined.hmm file")
        return [str(combined_hmm)]
    raise ProcessingError("No HMM model files found", step="hmm_search")

def _get_hmm_output_paths(hmmout_dir: Path) -> Tuple[Path, Path, Path, Path]:
    models_out = hmmout_dir / "models.out"
    models_out_filtered = hmmout_dir / "models.out.filtered"
    models_counts = hmmout_dir / "models.counts"
    models_score = hmmout_dir / "models.score"
    return models_out, models_out_filtered, models_counts, models_score

def _run_hmm_search(
    protein_file: Path,
    hmm_files: List[str],
    models_out: Path,
    models_out_filtered: Path,
    models_counts: Path,
    models_score: Path,
    threads: int,
    sensitive_mode: bool,
) -> None:
    if len(hmm_files) > 1:
        from src.core.hmm_search import generate_counts_from_output
        from src.core.hmm_search_multi import run_pyhmmer_search_multi

        run_pyhmmer_search_multi(
            hmm_files=hmm_files,
            query_file=str(protein_file),
            output_file=str(models_out),
            threads=threads,
            sensitive_mode=sensitive_mode,
        )
        shutil.copy2(str(models_out), str(models_out_filtered))
        generate_counts_from_output(str(models_out_filtered), str(models_counts))
        return

    run_pyhmmer_search_with_filtering(
        query_file=str(protein_file),
        hmm_file=hmm_files[0],
        output_file=str(models_out),
        filtered_output=str(models_out_filtered),
        cutoffs_file=str(None),
        counts_file=str(models_counts),
        score_file=str(models_score),
        threads=threads,
        sensitive_mode=sensitive_mode,
    )

def _run_blast_search(
    query_name: str,
    query_output_dir: Path,
    database_path: Path,
    protein_file: Path,
    hmm_results: Any,
    threads: int,
    logger,
) -> None:
    logger.info(f"Running BLAST search: {query_name}")
    blast_dir = query_output_dir / "blastp_out"
    blast_dir.mkdir(exist_ok=True)
    if not hmm_results:
        return

    from src.core.blast import run_blastp

    for marker in hmm_results:
        ref_faa = database_path / "database" / "alignment" / f"{marker}.ref.faa"
        if not ref_faa.exists():
            continue
        blast_out = blast_dir / f"{query_name}.{marker}.blastpout"
        run_blastp(
            queryfaa=str(protein_file),
            refdb=str(ref_faa),
            blastpout=str(blast_out),
            threads=threads,
        )

def _process_markers(
    query_name: str,
    models_out_filtered: Path,
    protein_file: Path,
    query_output_dir: Path,
    database_path: Path,
    tree_method: str,
    mode_fast: bool,
    threads: int,
    logger,
) -> Dict[str, Dict[str, Any]]:
    logger.info(f"Extracting and processing markers: {query_name}")
    marker_files = _extract_marker_files(models_out_filtered, protein_file, query_output_dir, mode_fast)
    if not marker_files:
        return {}
    max_parallel, threads_per_marker = _calculate_marker_threading(marker_files, threads)
    logger.info(
        f"Processing {len(marker_files)} markers: {max_parallel} parallel × {threads_per_marker} threads each"
    )
    return _run_marker_tasks(
        marker_files,
        database_path,
        query_output_dir,
        tree_method,
        max_parallel,
        threads_per_marker,
        logger,
    )

def _extract_marker_files(
    models_out_filtered: Path, protein_file: Path, query_output_dir: Path, mode_fast: bool
) -> List[Tuple[str, Path]]:
    query_hits_dir = query_output_dir / "query_hits_faa"
    query_hits_dir.mkdir(exist_ok=True)
    if not models_out_filtered.exists():
        return []
    return extract_marker_hits(
        hmm_output=models_out_filtered,
        query_faa=protein_file,
        output_dir=query_hits_dir,
        mode_fast=mode_fast,
    )

def _calculate_marker_threading(
    marker_files: List[Tuple[str, Path]], threads: int
) -> Tuple[int, int]:
    num_markers = len(marker_files)
    if threads >= num_markers * 2:
        return num_markers, min(threads // num_markers, 4)
    return min(threads, num_markers), 1

def _run_marker_tasks(
    marker_files: List[Tuple[str, Path]],
    database_path: Path,
    query_output_dir: Path,
    tree_method: str,
    max_parallel: int,
    threads_per_marker: int,
    logger,
) -> Dict[str, Dict[str, Any]]:
    marker_results: Dict[str, Dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = [
            executor.submit(
                _process_single_marker,
                marker_data,
                database_path,
                query_output_dir,
                tree_method,
                threads_per_marker,
                logger,
            )
            for marker_data in marker_files
        ]
        for future in as_completed(futures):
            try:
                marker, result = future.result()
                if marker and result:
                    marker_results[marker] = result
                    logger.info(f"Completed marker: {marker}")
            except Exception as exc:
                logger.error(f"Error processing marker: {exc}")
    return marker_results

def _process_single_marker(
    marker_data: Tuple[str, Path],
    database_path: Path,
    query_output_dir: Path,
    tree_method: str,
    threads_per_marker: int,
    logger,
) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    marker, marker_faa = marker_data
    logger.debug(f"Processing marker {marker}, file: {marker_faa}")
    if not marker_faa.exists():
        logger.warning(f"Marker file does not exist: {marker_faa}")
        return None, None
    if not _marker_database_exists(marker, database_path, logger):
        return None, None

    processor = MarkerProcessor(marker, database_path, query_output_dir)
    result = processor.process_marker(
        marker_faa,
        max_blast_hits=100,
        tree_method=tree_method,
        threads=threads_per_marker,
    )
    logger.debug(f"Marker {marker} processing result: {result is not None}")
    if result and result.get("tree"):
        logger.info(f"Successfully processed marker {marker} with tree file: {result['tree']}")
        return marker, result
    logger.warning(f"Marker {marker} processing did not produce a tree file")
    return None, None

def _marker_database_exists(marker: str, database_path: Path, logger) -> bool:
    try:
        marker_db = get_marker_database_path(marker, database_path)
        if marker_db.exists():
            return True
        logger.warning(f"Marker database does not exist: {marker_db}")
        return False
    except FileNotFoundError:
        logger.warning(f"Marker database not found for {marker}")
        return False

def _analyze_marker_trees(
    query_name: str,
    query_output_dir: Path,
    database_path: Path,
    marker_results: Dict[str, Dict[str, Any]],
    logger,
) -> Dict[str, Any]:
    logger.info(f"Marker results: {len(marker_results)} markers processed")
    if not marker_results:
        return {}
    tree_dir = query_output_dir / "queryrefs_genetrees"
    logger.info(f"Looking for tree files in {tree_dir}")
    logger.info(f"Tree directory exists: {tree_dir.exists()}")
    if not tree_dir.exists():
        logger.warning(f"Tree directory does not exist: {tree_dir}")
        return {}
    tree_files = list(tree_dir.glob("*.treefile"))
    logger.info(f"Found {len(tree_files)} tree files")
    if not tree_files:
        logger.warning("Tree directory exists but no .treefile files found")
        return {}
    _log_tree_files(tree_files, logger)
    return _run_tree_analysis(query_name, query_output_dir, database_path, tree_dir, logger)

def _log_tree_files(tree_files: List[Path], logger) -> None:
    for index, tree_file in enumerate(tree_files[:5]):
        logger.info(f"  Tree file {index + 1}: {tree_file.name}")

def _run_tree_analysis(
    query_name: str, query_output_dir: Path, database_path: Path, tree_dir: Path, logger
) -> Dict[str, Any]:
    labels_file = database_path / "gvclassFeb26_labels.tsv"
    logger.info(f"Labels file: {labels_file}, exists: {labels_file.exists()}")
    stats_dir = query_output_dir / "stats"
    stats_dir.mkdir(exist_ok=True, parents=True)
    tree_nn_file = stats_dir / f"{query_name}.tree_nn"
    try:
        analyzer = TreeAnalyzer(labels_file)
        logger.info(f"TreeAnalyzer initialized with labels file: {labels_file}")
        logger.info(f"Created stats directory: {stats_dir}")
        logger.info(f"Will write tree_nn results to {tree_nn_file}")
        tree_nn_results = analyzer.process_marker_trees(tree_dir, query_name, tree_nn_file)
        logger.info(f"Tree analysis complete. Results: {len(tree_nn_results)} markers")
    except Exception as exc:
        logger.error(f"Tree analysis failed: {exc}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        tree_nn_results = {}
    _log_tree_nn_output(tree_nn_file, logger)
    return tree_nn_results

def _log_tree_nn_output(tree_nn_file: Path, logger) -> None:
    if tree_nn_file.exists():
        logger.info(f"Tree_nn file created successfully: {tree_nn_file}")
        logger.info(f"Tree_nn file size: {tree_nn_file.stat().st_size} bytes")
        return
    logger.error(f"Tree_nn file was NOT created: {tree_nn_file}")

def _generate_summary_data(
    query_name: str,
    query_output_dir: Path,
    database_path: Path,
    tree_nn_results: Dict[str, Any],
    mode_fast: bool,
    logger,
) -> Dict[str, Any]:
    logger.info(f"Generating full summary: {query_name}")
    try:
        summarizer = FullSummarizer(database_path)
        summary_data = summarizer.summarize_query_full(
            query_id=query_name,
            query_output_dir=query_output_dir,
            tree_nn_results=tree_nn_results,
            mode_fast=mode_fast,
        )
        logger.info(f"Summary generated successfully for {query_name}")
        return summary_data
    except Exception as exc:
        logger.error(f"Failed to generate summary for {query_name}: {exc}")
        logger.error(traceback.format_exc())
        return _fallback_summary_data(query_name)

def _fallback_summary_data(query_name: str) -> Dict[str, Any]:
    return {
        "query": query_name,
        "length": "",
        "genecount": "",
        "coding_density": "",
        "ttable": "",
        "classification": "Error during classification",
    }

def _create_output_dirs(query_output_dir: Path) -> None:
    for dir_name in ["query_faa", "query_fna", "query_gff"]:
        (query_output_dir / dir_name).mkdir(exist_ok=True, parents=True)

def _copy_final_sequence_outputs(
    query_file: Path,
    query_name: str,
    query_output_dir: Path,
    reformatted_file: Path,
    protein_file: Path,
    best_code: Optional[int],
    logger,
) -> None:
    if query_file.suffix != ".fna" or best_code is None:
        return
    query_faa_dir = query_output_dir / "query_faa"
    query_fna_dir = query_output_dir / "query_fna"
    query_gff_dir = query_output_dir / "query_gff"
    if protein_file.exists():
        shutil.copy2(protein_file, query_faa_dir / f"{query_name}.faa")
    if reformatted_file.exists():
        shutil.copy2(reformatted_file, query_fna_dir / f"{query_name}.fna")
    gff_file = _find_query_gff_file(query_output_dir, query_name)
    if gff_file:
        shutil.copy2(gff_file, query_gff_dir / f"{query_name}.gff")
    else:
        logger.warning(f"No GFF file found for {query_name}")

def _find_query_gff_file(query_output_dir: Path, query_name: str) -> Optional[Path]:
    gff_candidates = [
        query_output_dir
        / "gene_calling"
        / f"{query_name}_reformatted"
        / f"{query_name}_reformatted.gff",
        query_output_dir / "gene_calling" / f"{query_name}.gff",
    ]
    for gff_candidate in gff_candidates:
        if gff_candidate.exists():
            return gff_candidate
    for gff_file in query_output_dir.glob("**/*.gff"):
        return gff_file
    return None

def _post_process_query(
    query_name: str, query_output_dir: Path, output_base: Path, logger
) -> None:
    logger.info(f"Post-processing results for query {query_name}")
    logger.debug(f"Query output directory: {query_output_dir}")
    logger.debug(f"Output base directory: {output_base}")
    post_process_start = time.time()
    try:
        _copy_query_summary_tab(query_name, query_output_dir, output_base, logger)
        _log_tree_nn_file(query_name, query_output_dir, logger)
        _archive_query_output(query_name, query_output_dir, output_base, logger)
        _remove_reformatted_file(query_name, output_base, logger)
        shutil.rmtree(query_output_dir)
        logger.info(f"✓ Removed original directory: {query_output_dir}")
        elapsed = time.time() - post_process_start
        logger.info(f"✓ Post-processing completed in {elapsed:.2f} seconds")
    except Exception as exc:
        elapsed = time.time() - post_process_start
        logger.error(
            f"❌ Post-processing failed for {query_name} after {elapsed:.2f} seconds: {exc}"
        )
        logger.error(f"Traceback: {traceback.format_exc()}")
        _log_partial_outputs(query_name, query_output_dir, output_base, logger)
        raise RuntimeError(f"Post-processing failed for {query_name}: {exc}") from exc

def _copy_query_summary_tab(
    query_name: str, query_output_dir: Path, output_base: Path, logger
) -> None:
    summary_tab_file = query_output_dir / f"{query_name}.summary.tab"
    logger.debug(f"Looking for summary file: {summary_tab_file}")
    if summary_tab_file.exists():
        file_size = summary_tab_file.stat().st_size
        logger.debug(f"Found summary file, size: {file_size} bytes")
        destination = output_base / f"{query_name}.summary.tab"
        shutil.copy2(summary_tab_file, destination)
        logger.info(f"✓ Copied summary.tab to {destination} ({file_size} bytes)")
        return
    summary_tab_file_stats = query_output_dir / "stats" / f"{query_name}.summary.tab"
    if summary_tab_file_stats.exists():
        destination = output_base / f"{query_name}.summary.tab"
        shutil.copy2(summary_tab_file_stats, destination)
        logger.info(f"Copied summary.tab from stats dir to {destination}")
        return
    logger.warning(
        f"Summary tab file not found in {summary_tab_file} or {summary_tab_file_stats}"
    )

def _log_tree_nn_file(query_name: str, query_output_dir: Path, logger) -> None:
    tree_nn_file = query_output_dir / "stats" / f"{query_name}.tree_nn"
    if tree_nn_file.exists():
        file_size = tree_nn_file.stat().st_size
        logger.debug(f"Found tree_nn file in stats directory: {file_size} bytes")

def _archive_query_output(
    query_name: str, query_output_dir: Path, output_base: Path, logger
) -> None:
    tar_file = output_base / f"{query_name}.tar.gz"
    logger.debug(f"Creating archive: {tar_file}")
    file_count = sum(1 for file_path in query_output_dir.rglob("*") if file_path.is_file())
    logger.debug(f"Archiving {file_count} files from {query_output_dir}")
    with tarfile.open(tar_file, "w:gz") as tar_handle:
        tar_handle.add(query_output_dir, arcname=query_name)
    archive_size = tar_file.stat().st_size
    logger.info(f"✓ Created archive: {tar_file} ({archive_size:,} bytes, {file_count} files)")

def _remove_reformatted_file(query_name: str, output_base: Path, logger) -> None:
    reformatted_file = output_base / f"{query_name}_reformatted.fna"
    if reformatted_file.exists():
        reformatted_file.unlink()
        logger.info(f"✓ Removed reformatted file: {reformatted_file}")

def _log_partial_outputs(
    query_name: str, query_output_dir: Path, output_base: Path, logger
) -> None:
    logger.debug("Checking partial outputs:")
    if (output_base / f"{query_name}.summary.tab").exists():
        logger.debug("  - Summary file was copied")
    if (output_base / f"{query_name}.tar.gz").exists():
        logger.debug("  - Archive was created")
    if query_output_dir.exists():
        logger.debug("  - Original directory still exists")

def _build_query_result(
    query_name: str,
    summary_data: Dict[str, Any],
    best_code: Optional[int],
    output_base: Path,
) -> Dict[str, Any]:
    return {
        "query": query_name,
        "status": "complete",
        "summary_data": summary_data,
        "genetic_code": best_code if best_code is not None else "N/A",
        "output_dir": str(output_base),
    }
