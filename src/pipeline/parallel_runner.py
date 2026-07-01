"""GVClass pipeline execution engine."""

from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import logging
import traceback
from concurrent.futures import Future, ThreadPoolExecutor, as_completed

from src.pipeline.query_processor import run_query_processing
from src.pipeline.progress_events import emit_progress_event
from src.pipeline.run_status import (
    RunStatusRecorder,
    query_completion_state_in_run_status,
)
from src.pipeline.summary_writer import (
    archive_final_summary_extended_files,
    read_individual_summary_results,
    write_final_summary_extended_files,
    write_final_summary_files,
)
from src.utils import InputValidator
from src.utils.database_manager import DatabaseManager


def validate_and_setup_task(
    query_dir: str,
    output_dir: str,
    database_path: Optional[str] = None,
    allow_short: bool = False,
) -> Dict[str, Any]:
    """Validate inputs and setup directories."""
    logger = logging.getLogger("gvclass_runner")
    logger.info(f"Validating inputs - Query: {query_dir}, Output: {output_dir}")
    query_path = InputValidator.validate_query_directory(
        query_dir, allow_short=allow_short
    )
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    db_path = DatabaseManager.setup_database(database_path)
    logger.info(f"Using database at: {db_path}")

    # InputValidator already accepts .fna / .faa / .fasta / .fas via
    # alphabet inference, so enumerate the same set here — otherwise a
    # directory of .fasta inputs would validate but yield zero queries
    # (Codex-audit finding).
    query_files = []
    for ext in (".fna", ".faa", ".fasta", ".fas"):
        query_files.extend(query_path.glob(f"*{ext}"))

    logger.info(f"Found {len(query_files)} query files")
    return {
        "query_path": query_path,
        "output_path": output_path,
        "database_path": db_path,
        "query_files": query_files,
    }


def process_query_task(
    query_file: Path,
    output_base: Path,
    database_path: Path,
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str = "legacy",
    sensitive_mode: bool = False,
    threads: int = 4,
    species_tree: bool = False,
    species_tree_trim: str = "witchi",
    run_status: Optional[RunStatusRecorder] = None,
) -> Dict[str, Any]:
    """
    Process a single query through the entire pipeline.
    This is the main unit of parallelization.
    """
    return run_query_processing(
        query_file=query_file,
        output_base=output_base,
        database_path=database_path,
        genetic_codes=genetic_codes,
        tree_method=tree_method,
        mode_fast=mode_fast,
        completeness_mode=completeness_mode,
        sensitive_mode=sensitive_mode,
        threads=threads,
        species_tree=species_tree,
        species_tree_trim=species_tree_trim,
        run_status=run_status,
    )


def create_final_summary_task(results: List[Dict[str, Any]], output_dir: Path) -> Path:
    """Create the final summary file combining all results."""
    logger = logging.getLogger("gvclass_runner")
    logger.info(f"Creating final summary for {len(results)} queries")
    summary_tsv = write_final_summary_files(results, output_dir)
    summary_csv = output_dir / "gvclass_summary.csv"
    extended_tsv = write_final_summary_extended_files(results, output_dir)
    extended_archive = archive_final_summary_extended_files(output_dir)
    logger.info(f"Summary written to: {summary_tsv}")
    logger.info(f"CSV summary written to: {summary_csv}")
    logger.info(f"Extended diagnostics archived in: {extended_archive or extended_tsv}")
    logger.info("All query processing complete")
    logger.info("Post-processing complete")
    return summary_tsv


def calculate_optimal_workers(
    n_queries: int,
    total_threads: int,
    min_threads_per_worker: int = 4,
    max_workers: Optional[int] = None,
) -> Tuple[int, int]:
    """
    Calculate optimal number of workers and threads per worker.

    Strategy:
    - Each worker needs at least min_threads_per_worker threads
    - Maximize parallelism while maintaining thread efficiency
    - For few queries with many threads, use fewer workers with more threads
    - For many queries with limited threads, use more workers with fewer threads
    """
    if max_workers is None:
        max_workers = n_queries

    max_possible_workers = min(
        n_queries,
        total_threads // min_threads_per_worker,
        max_workers,
    )

    if n_queries <= 4 and total_threads >= n_queries * 8:
        n_workers = n_queries
        threads_per_worker = total_threads // n_workers
    else:
        target_threads_per_worker = 4
        n_workers = total_threads // target_threads_per_worker
        n_workers = min(n_workers, max_possible_workers)
        if n_queries > 100:
            n_workers = min(n_workers, 20)
        n_workers = max(1, n_workers)
        threads_per_worker = total_threads // n_workers

    return max(1, n_workers), max(min_threads_per_worker, threads_per_worker)


def _query_is_resume_complete(output_path: Path, query_name: str, logger) -> bool:
    """Return True when ``query_name`` is safe to skip under ``--resume``.

    Preferred gate: the consolidated ``run_status.json`` manifest. For
    backward compatibility, v1.4.3+ ``<query_name>.SUCCESS`` files are still
    accepted, followed by the pre-1.4.3 ``<query_name>.summary.tab`` +
    ``<query_name>.tar.gz`` fallback.
    """
    import tarfile

    manifest_state = query_completion_state_in_run_status(output_path, query_name)
    if manifest_state is not None:
        return manifest_state

    sentinel = output_path / f"{query_name}.SUCCESS"
    if sentinel.exists():
        logger.info(
            f"Resume mode: accepting legacy SUCCESS sentinel for {query_name}. "
            "New runs write run_status.json instead."
        )
        return True

    summary_file = output_path / f"{query_name}.summary.tab"
    tar_file = output_path / f"{query_name}.tar.gz"
    if not (summary_file.exists() and tar_file.exists()):
        return False
    try:
        if not tarfile.is_tarfile(tar_file):
            logger.warning(
                f"Legacy resume fallback: tar file is malformed, re-processing {query_name}"
            )
            return False
    except (OSError, tarfile.TarError) as exc:
        logger.warning(
            f"Legacy resume fallback: cannot verify {tar_file} ({exc}); re-processing {query_name}"
        )
        return False
    logger.info(
        f"Legacy resume fallback: accepting pre-v1.4.3 outputs for {query_name} "
        "(not recorded in run_status.json). Re-run without --resume to upgrade."
    )
    return True


def _apply_resume_filter(config: Dict[str, Any], resume: bool, logger) -> None:
    config["resume_skipped_query_names"] = []
    if not resume:
        return

    output_path = config["output_path"]
    filtered_queries = []
    skipped_count = 0
    for query_file in list(config["query_files"]):
        query_name = query_file.stem
        if _query_is_resume_complete(output_path, query_name, logger):
            logger.info(f"Skipping completed query: {query_name}")
            config["resume_skipped_query_names"].append(query_name)
            skipped_count += 1
        else:
            filtered_queries.append(query_file)

    config["query_files"] = filtered_queries
    logger.info(
        f"Resume mode: skipped {skipped_count} completed queries, {len(filtered_queries)} remaining"
    )


def _resolve_worker_distribution(
    n_queries: int,
    total_threads: int,
    max_workers: Optional[int],
    threads_per_worker: Optional[int],
) -> Tuple[int, int]:
    if threads_per_worker is None:
        return calculate_optimal_workers(
            n_queries=n_queries,
            total_threads=total_threads,
            max_workers=max_workers,
        )

    if threads_per_worker <= 0:
        raise ValueError("threads_per_worker must be a positive integer")

    max_by_threads = max(1, total_threads // threads_per_worker)
    n_workers = min(n_queries, max_by_threads, max_workers or n_queries)
    return n_workers, threads_per_worker


def _submit_query_jobs(
    executor: ThreadPoolExecutor,
    query_files: List[Path],
    config: Dict[str, Any],
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive_mode: bool,
    threads_per_worker: int,
    species_tree: bool = False,
    species_tree_trim: str = "witchi",
    run_status: Optional[RunStatusRecorder] = None,
) -> Dict[Future, Path]:
    future_to_query: Dict[Future, Path] = {}
    for query_file in query_files:
        future = executor.submit(
            process_query_task,
            query_file=query_file,
            output_base=config["output_path"],
            database_path=config["database_path"],
            genetic_codes=genetic_codes,
            tree_method=tree_method,
            mode_fast=mode_fast,
            completeness_mode=completeness_mode,
            sensitive_mode=sensitive_mode,
            threads=threads_per_worker,
            species_tree=species_tree,
            species_tree_trim=species_tree_trim,
            run_status=run_status,
        )
        future_to_query[future] = query_file
    return future_to_query


def _collect_query_results(
    future_to_query: Dict[Future, Path],
    logger,
    run_status: Optional[RunStatusRecorder] = None,
) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for future in as_completed(future_to_query):
        query_file = future_to_query[future]
        try:
            result = future.result()
            results.append(result)
            logger.info(f"Completed: {result['query']} - {result['status']}")
        except Exception as exc:
            logger.error(f"Query {query_file.stem} failed: {exc}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            if run_status is not None:
                run_status.record_query_failed(query_file.stem, str(exc))
            emit_progress_event(
                {
                    "event": "query_failed",
                    "query": query_file.stem,
                    "progress": 100,
                    "stage": "failed",
                }
            )
            results.append(
                {"query": query_file.stem, "status": "failed", "error": str(exc)}
            )
    return results


def _process_queries_in_parallel(
    config: Dict[str, Any],
    n_workers: int,
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive_mode: bool,
    threads_per_worker: int,
    logger,
    species_tree: bool = False,
    species_tree_trim: str = "witchi",
    run_status: Optional[RunStatusRecorder] = None,
) -> Tuple[List[Dict[str, Any]], int]:
    query_files = list(config["query_files"])
    total_queries = len(query_files)
    logger.info(
        f"Processing {total_queries} queries with {n_workers} workers ({threads_per_worker} threads each)"
    )

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        future_to_query = _submit_query_jobs(
            executor,
            query_files,
            config,
            genetic_codes,
            tree_method,
            mode_fast,
            completeness_mode,
            sensitive_mode,
            threads_per_worker,
            species_tree=species_tree,
            species_tree_trim=species_tree_trim,
            run_status=run_status,
        )
        results = _collect_query_results(future_to_query, logger, run_status)

    return results, total_queries


def _count_completed(results: List[Dict[str, Any]]) -> int:
    return len([result for result in results if result["status"] == "complete"])


def _combine_with_resume_skipped_results(
    config: Dict[str, Any], results: List[Dict[str, Any]], logger
) -> List[Dict[str, Any]]:
    skipped_names = config.get("resume_skipped_query_names", [])
    if not skipped_names:
        return results

    skipped_results = read_individual_summary_results(
        config["output_path"], skipped_names
    )
    combined_results = {result["query"]: result for result in skipped_results}
    combined_results.update({result["query"]: result for result in results})

    missing_names = sorted(set(skipped_names) - set(combined_results))
    if missing_names:
        logger.warning(
            "Resume mode: skipped queries missing readable summary files: %s",
            ", ".join(missing_names),
        )
        # A query that resume declared complete but whose
        # loose summary cannot be re-read must still appear in the run, as a
        # failed row, rather than silently disappearing from the final summary.
        for name in missing_names:
            combined_results[name] = {"query": name, "status": "missing_summary"}

    return list(combined_results.values())


def gvclass_flow(
    query_dir: str,
    output_dir: str,
    database_path: Optional[str] = None,
    total_threads: int = 16,
    max_workers: Optional[int] = None,
    threads_per_worker: Optional[int] = None,
    tree_method: str = "veryfasttree",
    iqtree_mode: str = "fast",
    mode_fast: bool = True,
    completeness_mode: str = "legacy",
    sensitive_mode: bool = False,
    genetic_codes: List[int] = [0, 1, 4, 6, 11, 15, 29, 106, 129],
    cluster_type: str = "local",
    cluster_config: Optional[Dict[str, Any]] = None,
    resume: bool = False,
    allow_short: bool = False,
    species_tree: bool = False,
    species_tree_combined: bool = False,
    species_tree_trim: str = "witchi",
):
    """Main GVClass pipeline."""
    logger = logging.getLogger("gvclass_runner")
    if species_tree_combined:
        species_tree = True  # --species-tree-combined implies --species-tree
    from src.core.alignment import configure_iqtree

    configure_iqtree(iqtree_mode)
    logger.info("Validating inputs and setting up directories")
    config = validate_and_setup_task(
        query_dir, output_dir, database_path, allow_short=allow_short
    )
    all_query_names = [query_file.stem for query_file in config["query_files"]]
    _apply_resume_filter(config, resume, logger)

    n_queries = len(config["query_files"])
    logger.info(f"Found {n_queries} queries to process")

    run_status = RunStatusRecorder(
        output_dir=config["output_path"],
        database_path=config["database_path"],
        query_names=all_query_names,
        resume_skipped_query_names=config.get("resume_skipped_query_names", []),
        settings={
            "tree_method": tree_method,
            "iqtree_mode": iqtree_mode,
            "mode_fast": mode_fast,
            "completeness_mode": completeness_mode,
            "sensitive_mode": sensitive_mode,
            "total_threads": total_threads,
            "max_workers": max_workers,
            "threads_per_worker": threads_per_worker,
            "resume": resume,
            "species_tree": species_tree,
            "species_tree_combined": species_tree_combined,
            "species_tree_trim": species_tree_trim,
        },
    )
    run_status.start()

    # When --resume is set and every query is already complete in the manifest
    # (or a valid legacy marker), we have nothing to process. Creating a
    # ThreadPoolExecutor with max_workers=0 would raise ValueError, so short-
    # circuit and emit the summary from whatever is already on disk.
    if n_queries == 0:
        logger.info(
            "No queries remaining after resume filter; skipping parallel processing."
        )
        if species_tree_combined:
            logger.info(
                "All queries resume-skipped; the combined species tree is built only "
                "from this run's queries, so it is not regenerated. Re-run without "
                "--resume to rebuild out/species_tree/combined.*"
            )
        all_results = _combine_with_resume_skipped_results(config, [], logger)
        summary_file = create_final_summary_task(all_results, config["output_path"])
        run_status.finish("complete")
        logger.info(f"Pipeline completed! Results in: {config['output_path']}")
        logger.info(f"Summary written to: {summary_file}")
        _ = cluster_type
        _ = cluster_config
        return summary_file

    n_workers, threads_per_worker = _resolve_worker_distribution(
        n_queries, total_threads, max_workers, threads_per_worker
    )
    logger.info(
        f"Parallelization strategy: {n_workers} workers × {threads_per_worker} threads"
    )

    # Clear any previous species-tree scratch + outputs and record this run's
    # queries so the per-query hook and the opt-in combined step read exactly this
    # run's sidecars. Wrapped so a scratch/filesystem failure disables the feature
    # rather than aborting the run.
    if species_tree:
        skipped = config.get("resume_skipped_query_names", [])
        if resume and skipped:
            logger.warning(
                "%d resume-skipped queries will have no species-tree placement this run "
                "(per-query trees + species_tree_* columns are produced only for queries "
                "processed this run); re-run those without --resume to place them.",
                len(skipped),
            )
        try:
            from src.core.species_tree.handoff import init_handoff

            init_handoff(
                Path(config["output_path"]),
                [query_file.stem for query_file in config["query_files"]],
                resume=resume,
            )
        except Exception as exc:
            logger.warning(
                "Species-tree handoff init failed (%s); disabling --species-tree for this run",
                exc,
            )
            species_tree = False
            species_tree_combined = False

    logger.info("Starting parallel query processing")
    results, total_queries = _process_queries_in_parallel(
        config,
        n_workers,
        genetic_codes,
        tree_method,
        mode_fast,
        completeness_mode,
        sensitive_mode,
        threads_per_worker,
        logger,
        species_tree=species_tree,
        species_tree_trim=species_tree_trim,
        run_status=run_status,
    )
    completed_count = _count_completed(results)
    logger.info(
        f"All queries complete. Successfully processed {completed_count}/{total_queries} queries"
    )

    # Opt-in combined species tree from this run's sidecars (--species-tree-combined).
    # The per-query hook already owns the summary columns, so this writes only the
    # additional <out>/species_tree/ combined artifacts. Never blocks the summary
    # (internally wrapped). Resume-skipped queries are excluded (they have no
    # sidecar this run) — a documented limitation.
    if species_tree_combined:
        try:
            from src.core.species_tree.orchestration import run_combined_species_tree

            run_combined_species_tree(
                Path(config["output_path"]),
                Path(config["database_path"]),
                results,
                total_threads,
                species_tree_trim,
                tree_method,
            )
        except Exception as exc:
            logger.warning("Combined species-tree step failed: %s", exc)

    logger.info("Creating final summary...")
    all_results = _combine_with_resume_skipped_results(config, results, logger)
    summary_file = create_final_summary_task(all_results, config["output_path"])
    run_status.finish("complete" if completed_count == total_queries else "failed")
    logger.info(f"Pipeline completed! Results in: {config['output_path']}")
    logger.info(f"Summary written to: {summary_file}")
    logger.info(f"Successfully processed {completed_count}/{len(results)} queries")

    _ = cluster_type
    _ = cluster_config
    return summary_file


# NOTE (v1.4.3 cleanup): the former create_local_deployment /
# create_slurm_deployment helpers were removed. They imported
# `prefect.deployments.Deployment`, were never wired to the CLI, and
# were the only source of a live Prefect dependency in this module.
# Renaming this file to `parallel_runner.py` (and the sibling CLI
# wrapper to `gvclass_runner`) is deferred to v1.5.0 because the CLI
# currently spawns `python -m src.bin.gvclass_runner` as a subprocess.
# If/when we truly adopt Prefect (@flow / @task / retry / persistence)
# file a follow-up issue.
