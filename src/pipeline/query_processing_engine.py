"""Internal helpers for per-query GVClass processing."""

from __future__ import annotations

import csv
import logging
import os
import shutil
import tarfile
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO

from src.core.genetic_code_optimizer import (
    META_CODE,
    GeneCallingResult,
    optimize_genetic_code,
)
from src.core.hmm_search import COMBINED_HMM, run_marker_search
from src.core.marker_extraction import extract_marker_hits, get_marker_database_path
from src.core.marker_processing import MarkerProcessor
from src.core.reformat import calculate_stats, reformat_sequences
from src.core.summarize_full import FullSummarizer
from src.core.tree_analysis import TreeAnalyzer
from src.pipeline.progress_events import emit_progress_event
from src.pipeline.run_config import RunConfig
from src.pipeline.summary_writer import (
    write_individual_final_summary_file,
    write_individual_summary_file,
)
from src.utils.error_handling import ProcessingError
from src.utils.input_validation import InputValidator
from src.utils.resource_store import ResourceStore

#: Genetic codes evaluated for every nucleotide query. This is deliberately not
#: configurable: GVClass always runs the full optimization and picks the
#: best-scoring code per query. Code 0 is pyrodigal meta mode.
GENETIC_CODES = [0, 1, 4, 6, 11, 15, 29, 106, 129]

#: Folder in the output directory that holds the working files of the queries that
#: run now: ``<output>/.gvclass_work/<query>/``. A query folder is archived and
#: removed when its query ends.
WORK_DIRNAME = ".gvclass_work"

_THREAD_LOCAL = threading.local()


def _pipeline_logger():
    return logging.getLogger("gvclass_runner")


@dataclass
class PreparedQueryInput:
    """Input files and metadata needed by downstream query steps."""

    protein_file: Path
    reformatted_file: Path
    best_code: Optional[int]
    is_nucleotide: bool = False


def process_single_query(
    query_file: Path, config: RunConfig, threads: int
) -> Dict[str, Any]:
    """Run the complete pipeline for one query and publish its archive.

    The archive ``<output>/<query>.tar.gz`` exists when this returns. The caller
    writes the completion record.
    """
    logger = _pipeline_logger()
    output_base = config.output_dir
    query_name, query_output_dir = _initialize_query(query_file, output_base)
    _emit_query_progress(query_name, 0, "started")
    logger.info(
        "Processing query: %s with %d threads (sensitive_mode=%s)",
        query_name,
        threads,
        config.sensitive_mode,
    )

    prepared_input, summary_data = _run_query_stages(
        query_file, query_name, query_output_dir, config, threads, logger
    )

    # Species-tree per-query hook (the default product): runs while the per-query
    # gene trees still exist (before _post_process_query archives + deletes the
    # dir) and merges this query's tree-placement columns into summary_data before
    # the per-query summary is written below. Failures must reach the query collector.
    if config.species_tree:
        _emit_query_progress(query_name, 93, "species_tree")
        from src.core.species_tree.orchestration import (  # noqa: PLC0415  # Load optional species-tree code only when requested.
            run_per_query_species_tree,
        )

        run_per_query_species_tree(
            query_name,
            query_output_dir,
            config.database_path,
            summary_data,
            output_base,
            threads,
            config.species_tree_trim,
            config.tree_method,
        )
    _emit_query_progress(query_name, 95, "writing_outputs")

    _write_contamination_candidates_file(
        query_output_dir, query_name, summary_data, logger
    )
    logger.info(
        "Writing summary file to %s", query_output_dir / f"{query_name}.summary.tab"
    )
    write_individual_summary_file(query_output_dir, query_name, summary_data)
    write_individual_final_summary_file(query_output_dir, query_name, summary_data)
    _emit_query_progress(query_name, 97, "post_processing")
    _post_process_query(query_name, query_output_dir, output_base, logger)
    _emit_query_progress(query_name, 100, "complete")
    return _build_query_result(
        query_name, summary_data, prepared_input.best_code, output_base
    )


def _run_query_stages(
    query_file: Path,
    query_name: str,
    query_output_dir: Path,
    config: RunConfig,
    threads: int,
    logger,
) -> Tuple[PreparedQueryInput, Dict[str, Any]]:
    prepared_input = _prepare_query_input(
        query_file, query_output_dir, config.database_path, threads, logger
    )
    _create_output_dirs(query_output_dir)
    _copy_final_sequence_outputs(
        query_file,
        query_name,
        query_output_dir,
        prepared_input.reformatted_file,
        prepared_input.protein_file,
        prepared_input.best_code,
        prepared_input.is_nucleotide,
        logger,
    )
    _emit_query_progress(query_name, 15, "input_prepared")
    models_out_filtered = _run_hmm_and_blast(
        query_name,
        query_output_dir,
        prepared_input.protein_file,
        config,
        threads,
        logger,
    )
    _emit_query_progress(query_name, 35, "marker_search_complete")
    marker_results = _process_markers(
        query_name,
        models_out_filtered,
        prepared_input.protein_file,
        query_output_dir,
        config,
        threads,
        logger,
    )
    _emit_query_progress(query_name, 78, "marker_trees_complete")
    tree_nn_results = _analyze_marker_trees(
        query_name, query_output_dir, config.database_path, marker_results, logger
    )
    _emit_query_progress(query_name, 86, "tree_analysis_complete")
    summary_data = _generate_summary_data(
        query_name, query_output_dir, config, tree_nn_results, logger
    )
    _emit_query_progress(query_name, 92, "summary_complete")
    return prepared_input, summary_data


def _emit_query_progress(
    query_name: str,
    progress: float,
    stage: str,
    **extra: Any,
) -> None:
    emit_progress_event(
        {
            "event": "query_progress",
            "query": query_name,
            "progress": max(0, min(100, progress)),
            "stage": stage,
            **extra,
        }
    )


def _initialize_query(query_file: Path, output_base: Path) -> tuple[str, Path]:
    """Return the query name and its empty private work directory.

    The directory is ``<output>/.gvclass_work/<query>/``. A directory that an earlier
    attempt left is deleted, so old search results cannot affect the new scores.
    """
    query_name = query_file.stem
    query_output_dir = output_base / WORK_DIRNAME / query_name
    if query_file.resolve().is_relative_to(query_output_dir.resolve()):
        raise ProcessingError(
            f"The working directory {query_output_dir} for query {query_name} "
            f"holds the input file {query_file}. Choose an output directory "
            "outside the input directory.",
            step="initialize_query",
            input_file=str(query_file),
        )
    if query_output_dir.exists():
        shutil.rmtree(query_output_dir)
    query_output_dir.mkdir(parents=True)
    return query_name, query_output_dir


def _prepare_query_input(
    query_file: Path,
    query_output_dir: Path,
    database_path: Path,
    threads: int,
    logger,
) -> PreparedQueryInput:
    if is_nucleotide_input(query_file):
        return _prepare_nucleotide_query(
            query_file, query_output_dir, database_path, threads, logger
        )
    return _prepare_protein_query(query_file, query_output_dir, logger)


def is_nucleotide_input(query_file: Path) -> bool:
    """Decide nucleotide vs protein by alphabet, not just the suffix.

    ``.fna`` is always nucleotide and ``.faa`` always protein. For the
    extension-agnostic ``.fasta``/``.fas``/``.fa``, which the validator admits as
    either, the alphabet comes from the first record, so a DNA ``.fasta`` is
    gene-called instead of being fed raw to the protein path.
    """
    suffix = query_file.suffix.lower()
    if suffix == ".fna":
        return True
    if suffix in (".fasta", ".fas", ".fa"):
        return InputValidator.infer_fasta_type(query_file) == "fna"
    return False


def _prepare_nucleotide_query(
    query_file: Path,
    query_output_dir: Path,
    database_path: Path,
    threads: int,
    logger,
) -> PreparedQueryInput:
    query_name = query_file.stem
    logger.info(f"Reformatting nucleotide sequences: {query_name}")
    reformatted_file = query_output_dir / f"{query_name}_reformatted.fna"
    records = list(SeqIO.parse(str(query_file), "fasta"))
    _write_reformatted_sequences(records, query_file, reformatted_file)
    _write_nucleotide_stats_tsv(records, query_file, query_output_dir, query_name)
    gene_calling = _optimize_genetic_code(
        query_name,
        reformatted_file,
        query_output_dir,
        database_path,
        threads,
        logger,
    )
    _write_genecalling_tables(query_output_dir, query_name, gene_calling)
    return PreparedQueryInput(
        protein_file=gene_calling.faa,
        reformatted_file=reformatted_file,
        best_code=gene_calling.best_code,
        is_nucleotide=True,
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
        is_nucleotide=False,
    )


def _write_reformatted_sequences(
    records, query_file: Path, reformatted_file: Path
) -> None:
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
    threads: int,
    logger,
) -> GeneCallingResult:
    logger.info(f"Optimizing genetic code for: {query_name}")
    return optimize_genetic_code(
        reformatted_file,
        query_output_dir,
        database_path / COMBINED_HMM,
        GENETIC_CODES,
        threads=threads,
    )


def _write_genecalling_tables(
    query_output_dir: Path,
    query_name: str,
    gene_calling: GeneCallingResult,
) -> None:
    _write_genecalling_comparison_table(query_output_dir, query_name, gene_calling)
    _write_selected_genecode_table(query_output_dir, query_name, gene_calling)


def _write_genecalling_comparison_table(
    query_output_dir: Path,
    query_name: str,
    gene_calling: GeneCallingResult,
) -> None:
    comparison_file = query_output_dir / "stats" / f"{query_name}.genecalling.tab"
    ranked = sorted(
        gene_calling.evaluations.values(),
        key=lambda item: (-item.complete_hits, -item.avg_best_hit_score),
    )
    with open(comparison_file, "w") as handle:
        handle.write(
            "ttable\tGenes\tBases\tCodingDensity\tcomplete_bestHits\tavg_bestHit_score\n"
        )
        handle.writelines(
            f"{item.code}\t{item.gene_count}\t{gene_calling.assembly.length_bp}\t"
            f"{item.coding_density:.3f}\t{item.complete_hits}\t"
            f"{item.avg_best_hit_score:.1f}\n"
            for item in ranked
        )


def _write_selected_genecode_table(
    query_output_dir: Path,
    query_name: str,
    gene_calling: GeneCallingResult,
) -> None:
    stats_file = query_output_dir / "stats" / f"{query_name}.stats.tab"
    selected = gene_calling.evaluations[gene_calling.best_code]
    ttable_str = (
        "codemeta" if gene_calling.best_code == META_CODE else str(selected.code)
    )
    with open(stats_file, "w") as handle:
        handle.write(f"ttable\t{ttable_str}\n")
        handle.write(f"genes\t{selected.gene_count}\n")
        handle.write(f"bases\t{gene_calling.assembly.length_bp}\n")
        handle.write(f"coding_density\t{selected.coding_density:.3f}\n")


def _write_protein_stats(
    query_output_dir: Path, query_name: str, reformatted_records
) -> None:
    stats_dir = query_output_dir / "stats"
    stats_dir.mkdir(exist_ok=True)
    protein_count = len(reformatted_records)
    total_aa_length = sum(len(record.seq) for record in reformatted_records)
    estimated_genome_length = int(total_aa_length * 3 * 1.1)
    _write_protein_stats_tab(
        stats_dir / f"{query_name}.stats.tab", protein_count, estimated_genome_length
    )
    _write_protein_stats_tsv(
        stats_dir / f"{query_name}_stats.tsv", protein_count, estimated_genome_length
    )


def _write_protein_stats_tab(
    stats_file: Path, protein_count: int, genome_length: int
) -> None:
    with open(stats_file, "w") as handle:
        handle.write("ttable\tno_fna\n")
        handle.write(f"genes\t{protein_count}\n")
        handle.write(f"bases\t{genome_length}\n")
        handle.write("coding_density\t0.0\n")


def _write_protein_stats_tsv(
    stats_file: Path, protein_count: int, genome_length: int
) -> None:
    with open(stats_file, "w") as handle:
        handle.write("contigs\tLENbp\tGCperc\tgenecount\tCODINGperc\tttable\n")
        handle.write(f"0\t{genome_length}\t0.0\t{protein_count}\t0.0\tno_fna\n")


def _run_hmm_and_blast(
    query_name: str,
    query_output_dir: Path,
    protein_file: Path,
    config: RunConfig,
    threads: int,
    logger,
) -> Path:
    """Run the per-query marker HMM search and return ``models.out.filtered``.

    The reference BLAST pass belongs to :class:`MarkerProcessor`.
    """
    logger.info(f"Running HMM search: {query_name}")
    return run_marker_search(
        protein_file,
        config.database_path / COMBINED_HMM,
        query_output_dir / "hmmout",
        threads=threads,
        sensitive=config.sensitive_mode,
    )


# Removed in v1.4.3: a former ``_run_blast_search`` helper ran a per-marker
# BLAST pass that emitted ``<query>.<marker>.blastpout`` files. Those files
# were never read — contamination scoring reads ``<marker>.m8`` files
# produced by :class:`MarkerProcessor` (see ``src/core/marker_processing.py``)
# which also feeds the tree-building pipeline. Removing the duplicate pass
# halves the per-query BLAST cost.


def _process_markers(
    query_name: str,
    models_out_filtered: Path,
    protein_file: Path,
    query_output_dir: Path,
    config: RunConfig,
    threads: int,
    logger,
) -> Dict[str, Dict[str, Any]]:
    logger.info(f"Extracting and processing markers: {query_name}")
    _emit_query_progress(query_name, 40, "extracting_markers")
    marker_files = _extract_marker_files(
        models_out_filtered, protein_file, query_output_dir, config.mode_fast
    )
    if not marker_files:
        _emit_query_progress(query_name, 75, "no_markers_found")
        return {}
    max_parallel, threads_per_marker = _calculate_marker_threading(
        marker_files, threads
    )
    logger.info(
        f"Processing {len(marker_files)} markers: {max_parallel} parallel × {threads_per_marker} threads each"
    )
    _emit_query_progress(
        query_name,
        45,
        "building_marker_trees",
        completed_markers=0,
        total_markers=len(marker_files),
    )
    return _run_marker_tasks(
        query_name,
        marker_files,
        query_output_dir,
        config,
        max_parallel,
        threads_per_marker,
        logger,
    )


def _extract_marker_files(
    models_out_filtered: Path,
    protein_file: Path,
    query_output_dir: Path,
    mode_fast: bool,
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
    query_name: str,
    marker_files: List[Tuple[str, Path]],
    query_output_dir: Path,
    config: RunConfig,
    max_parallel: int,
    threads_per_marker: int,
    logger,
) -> Dict[str, Dict[str, Any]]:
    marker_results: Dict[str, Dict[str, Any]] = {}
    total_markers = len(marker_files)
    completed_markers = 0
    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = [
            executor.submit(
                _process_single_marker,
                marker_data,
                query_output_dir,
                config,
                threads_per_marker,
                logger,
            )
            for marker_data in marker_files
        ]
        for future in as_completed(futures):
            completed_markers += 1
            marker, result = future.result()
            if marker and result:
                marker_results[marker] = result
                logger.info("Completed marker: %s", marker)
            _emit_query_progress(
                query_name,
                45 + (30 * completed_markers / total_markers),
                "building_marker_trees",
                marker=marker,
                completed_markers=completed_markers,
                total_markers=total_markers,
            )
    return marker_results


def _process_single_marker(
    marker_data: Tuple[str, Path],
    query_output_dir: Path,
    config: RunConfig,
    threads_per_marker: int,
    logger,
) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    marker, marker_faa = marker_data
    logger.debug(f"Processing marker {marker}, file: {marker_faa}")
    if not marker_faa.exists():
        raise FileNotFoundError(f"Marker input does not exist: {marker_faa}")
    if not _marker_database_exists(marker, config.database_path, logger):
        raise FileNotFoundError(f"Reference proteins are missing for marker {marker}")

    processor = MarkerProcessor(marker, config.database_path, query_output_dir)
    result = processor.process_marker(
        marker_faa,
        max_blast_hits=100,
        tree_method=config.tree_method,
        threads=threads_per_marker,
    )
    logger.debug(f"Marker {marker} processing result: {result is not None}")
    if result and result.get("tree"):
        logger.info(
            f"Successfully processed marker {marker} with tree file: {result['tree']}"
        )
        return marker, result
    if result and result.get("tree_skipped"):
        return None, None
    raise ProcessingError(
        f"Marker {marker} produced neither a tree nor an explicit skip",
        step="marker_tree",
    )


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
        raise FileNotFoundError(f"Generated tree directory is missing: {tree_dir}")
    tree_files = [Path(result["tree"]) for result in marker_results.values()]
    for tree_file in tree_files:
        if not tree_file.is_file():
            raise FileNotFoundError(f"Generated marker tree is missing: {tree_file}")
    logger.info("Analyzing %d trees from the current attempt", len(tree_files))
    _log_tree_files(tree_files, logger)
    return _run_tree_analysis(
        query_name, query_output_dir, database_path, tree_dir, tree_files, logger
    )


def _log_tree_files(tree_files: List[Path], logger) -> None:
    for index, tree_file in enumerate(tree_files[:5]):
        logger.info(f"  Tree file {index + 1}: {tree_file.name}")


def _run_tree_analysis(
    query_name: str,
    query_output_dir: Path,
    database_path: Path,
    tree_dir: Path,
    tree_files: list[Path],
    logger: logging.Logger,
) -> Dict[str, Any]:
    labels_file = ResourceStore(database_path).label_path("labels.tsv")
    logger.info(f"Labels file: {labels_file}, exists: {labels_file.exists()}")
    stats_dir = query_output_dir / "stats"
    stats_dir.mkdir(exist_ok=True, parents=True)
    tree_nn_file = stats_dir / f"{query_name}.tree_nn"
    try:
        analyzer = TreeAnalyzer(labels_file)
        logger.info(f"TreeAnalyzer initialized with labels file: {labels_file}")
        logger.info(f"Created stats directory: {stats_dir}")
        logger.info(f"Will write tree_nn results to {tree_nn_file}")
        tree_nn_results = analyzer.process_marker_trees(
            tree_dir, query_name, tree_nn_file, tree_paths=tree_files
        )
        logger.info(f"Tree analysis complete. Results: {len(tree_nn_results)} markers")
    except Exception as exc:
        logger.error(f"Tree analysis failed: {exc}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        raise
    _log_tree_nn_output(tree_nn_file, logger)
    return tree_nn_results


def _log_tree_nn_output(tree_nn_file: Path, logger) -> None:
    if tree_nn_file.exists():
        logger.info(f"Tree_nn file created successfully: {tree_nn_file}")
        logger.info(f"Tree_nn file size: {tree_nn_file.stat().st_size} bytes")
        return
    logger.error(f"Tree_nn file was NOT created: {tree_nn_file}")


def _get_thread_summarizer(config: RunConfig) -> FullSummarizer:
    cache: dict[tuple[Any, str, str, bool], FullSummarizer] = getattr(
        _THREAD_LOCAL, "summarizers", None
    )
    if cache is None:
        cache = {}
        _THREAD_LOCAL.summarizers = cache
    key = (
        FullSummarizer,
        str(config.database_path.resolve()),
        config.completeness_mode,
        config.sensitive_mode,
    )
    if key not in cache:
        cache[key] = FullSummarizer(
            config.database_path,
            completeness_mode=config.completeness_mode,
            sensitive_mode=config.sensitive_mode,
        )
    return cache[key]


def _generate_summary_data(
    query_name: str,
    query_output_dir: Path,
    config: RunConfig,
    tree_nn_results: Dict[str, Any],
    logger,
) -> Dict[str, Any]:
    logger.info(f"Generating full summary: {query_name}")
    try:
        summary_data = _get_thread_summarizer(config).summarize_query_full(
            query_id=query_name,
            query_output_dir=query_output_dir,
            tree_nn_results=tree_nn_results,
            mode_fast=config.mode_fast,
        )
        logger.info(f"Summary generated successfully for {query_name}")
        return summary_data
    except Exception as exc:
        logger.error(f"Failed to generate summary for {query_name}: {exc}")
        logger.error(traceback.format_exc())
        raise


def _create_output_dirs(query_output_dir: Path) -> None:
    for dir_name in ["query_faa", "query_fna", "query_gff"]:
        (query_output_dir / dir_name).mkdir(exist_ok=True, parents=True)


def _write_contamination_candidates_file(
    query_output_dir: Path,
    query_name: str,
    summary_data: Dict[str, Any],
    logger,
) -> Optional[Path]:
    """Write review evidence independently of the contamination score gate."""
    suppressed_type = summary_data.get("contamination_type", "clean") in {
        "clean",
        "uncertain",
        "uncertain_sensitive_mode",
    }
    candidates = (
        [] if suppressed_type else summary_data.get("_contamination_candidates", [])
    )
    reviews = summary_data.get("_contamination_reviews", [])
    stats_dir = query_output_dir / "stats"
    fieldnames = [
        "query",
        "contamination_type",
        "estimated_contamination",
        "contig_id",
        "candidate_type",
        "reason",
        "recommended_action",
        "length_bp",
        "length_basis",
        "cellular_marker_count",
        "phage_marker_count",
        "viral_marker_count",
        "ambiguous_eve_hit_count",
        "nonviral_fraction",
        "foreign_viral_fraction",
    ]
    candidate_file = None
    for suffix, rows in (
        ("contamination_review", reviews),
        ("contamination_candidates", candidates),
    ):
        path = stats_dir / f"{query_name}.{suffix}.tsv"
        if not rows:
            if path.exists():
                path.unlink()
                logger.info("Removed stale contamination evidence file: %s", path)
            continue
        stats_dir.mkdir(exist_ok=True, parents=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in rows:
                row = {field: record.get(field, "") for field in fieldnames}
                row.update(
                    query=query_name,
                    contamination_type=summary_data.get("contamination_type", ""),
                    estimated_contamination=summary_data.get(
                        "estimated_contamination", ""
                    ),
                    recommended_action=record.get("recommended_action", "review"),
                )
                writer.writerow(row)
        logger.info("Contamination evidence file written: %s", path)
        if suffix == "contamination_candidates":
            candidate_file = path
    return candidate_file


def _copy_final_sequence_outputs(
    query_file: Path,
    query_name: str,
    query_output_dir: Path,
    reformatted_file: Path,
    protein_file: Path,
    best_code: Optional[int],
    is_nucleotide: bool,
    logger,
) -> None:
    # Gate on the inferred alphabet (not the literal .fna suffix) so a DNA
    # .fasta/.fas that was gene-called still gets its <name>.fna/.faa/.gff
    # written; protein inputs (.faa or AA .fasta) have no gene-calling output.
    if not is_nucleotide or best_code is None:
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
    gff_file = (
        query_output_dir
        / "gene_calling"
        / f"{query_name}_reformatted"
        / f"{query_name}_reformatted.gff"
    )
    return gff_file if gff_file.is_file() else None


def _post_process_query(
    query_name: str, query_output_dir: Path, output_base: Path, logger
) -> None:
    """Archive the work directory of a query and remove it.

    A failure leaves the work directory for inspection and any archive of an
    earlier run in place: the new archive replaces it only once it is complete.
    """
    logger.info(f"Post-processing results for query {query_name}")
    logger.debug(f"Query output directory: {query_output_dir}")
    logger.debug(f"Output base directory: {output_base}")
    post_process_start = time.time()
    try:
        _log_tree_nn_file(query_name, query_output_dir, logger)
        _remove_reformatted_file(query_name, query_output_dir, logger)
        _archive_query_output(query_name, query_output_dir, output_base, logger)
        shutil.rmtree(query_output_dir)
        logger.info("Removed original directory: %s", query_output_dir)
        elapsed = time.time() - post_process_start
        logger.info("Post-processing completed in %.2f seconds", elapsed)
    except Exception as exc:
        elapsed = time.time() - post_process_start
        logger.error(
            "Post-processing failed for %s after %.2f seconds: %s",
            query_name,
            elapsed,
            exc,
        )
        logger.error(f"Traceback: {traceback.format_exc()}")
        raise RuntimeError(f"Post-processing failed for {query_name}: {exc}") from exc


def _log_tree_nn_file(query_name: str, query_output_dir: Path, logger) -> None:
    tree_nn_file = query_output_dir / "stats" / f"{query_name}.tree_nn"
    if tree_nn_file.exists():
        file_size = tree_nn_file.stat().st_size
        logger.debug(f"Found tree_nn file in stats directory: {file_size} bytes")


def _archive_query_output(
    query_name: str, query_output_dir: Path, output_base: Path, logger
) -> None:
    """Write the per-query tarball atomically.

    The archive lands at ``<query_name>.tar.gz`` only after: (a) the tarball
    has been written to a sibling ``.tar.gz.part``, (b) the file handle has
    been fsynced, (c) ``tarfile.is_tarfile`` confirms the archive is well
    formed, (d) ``os.replace`` atomically swaps the final name into place,
    and (e) the parent directory is fsynced. A crash or SIGKILL at any
    point leaves either the previous final archive untouched or a
    ``.tar.gz.part`` that resume correctly ignores — never a truncated
    ``<query_name>.tar.gz`` that resume would wrongly accept.
    """
    tar_file = output_base / f"{query_name}.tar.gz"
    tar_tmp = output_base / f"{query_name}.tar.gz.part"
    if tar_tmp.exists():
        # Orphan from a prior crash; safe to discard before we open a fresh
        # writer at the same path.
        tar_tmp.unlink()
    logger.debug(f"Creating archive: {tar_file} (via {tar_tmp})")
    file_count = sum(
        1 for file_path in query_output_dir.rglob("*") if file_path.is_file()
    )
    logger.debug(f"Archiving {file_count} files from {query_output_dir}")
    with tarfile.open(tar_tmp, "w:gz") as tar_handle:
        tar_handle.add(query_output_dir, arcname=query_name)
    # fsync the data file before any rename so the payload is durable.
    with open(tar_tmp, "rb") as fh:
        os.fsync(fh.fileno())
    if not tarfile.is_tarfile(tar_tmp):
        tar_tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"Archive verification failed for {tar_tmp}; refusing to publish {tar_file}"
        )
    os.replace(tar_tmp, tar_file)
    # fsync the parent directory so the rename itself is durable on POSIX.
    try:
        dir_fd = os.open(str(tar_file.parent), os.O_RDONLY)
    except OSError:
        dir_fd = None
    if dir_fd is not None:
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    archive_size = tar_file.stat().st_size
    logger.info(
        "Created archive: %s (%s bytes, %s files)",
        tar_file,
        f"{archive_size:,}",
        file_count,
    )


def _remove_reformatted_file(
    query_name: str, query_output_dir: Path, logger: logging.Logger
) -> None:
    reformatted_file = query_output_dir / f"{query_name}_reformatted.fna"
    if reformatted_file.exists():
        reformatted_file.unlink()
        logger.info("Removed reformatted file: %s", reformatted_file)


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
