"""GVClass pipeline execution engine."""

from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import logging
import traceback
from concurrent.futures import Future, ThreadPoolExecutor, as_completed

from src.pipeline.query_processor import run_query_processing
from src.pipeline.summary_writer import write_final_summary_files
from src.utils import InputValidator
from src.utils.database_manager import DatabaseManager


def validate_and_setup_task(
    query_dir: str, output_dir: str, database_path: Optional[str] = None
) -> Dict[str, Any]:
    """Validate inputs and setup directories."""
    logger = logging.getLogger("gvclass_prefect")
    logger.info(f"Validating inputs - Query: {query_dir}, Output: {output_dir}")
    query_path = InputValidator.validate_query_directory(query_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    db_path = DatabaseManager.setup_database(database_path)
    logger.info(f"Using database at: {db_path}")

    query_files = []
    for ext in [".fna", ".faa"]:
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
    sensitive_mode: bool = False,
    threads: int = 4,
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
        sensitive_mode=sensitive_mode,
        threads=threads,
    )


def create_final_summary_task(results: List[Dict[str, Any]], output_dir: Path) -> Path:
    """Create the final summary file combining all results."""
    logger = logging.getLogger("gvclass_prefect")
    logger.info(f"Creating final summary for {len(results)} queries")
    summary_tsv = write_final_summary_files(results, output_dir)
    summary_csv = output_dir / "gvclass_summary.csv"
    logger.info(f"Summary written to: {summary_tsv}")
    logger.info(f"CSV summary written to: {summary_csv}")
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


def _apply_resume_filter(config: Dict[str, Any], resume: bool, logger) -> None:
    if not resume:
        return

    output_path = config["output_path"]
    filtered_queries = []
    skipped_count = 0
    for query_file in list(config["query_files"]):
        query_name = query_file.stem
        summary_file = output_path / f"{query_name}.summary.tab"
        tar_file = output_path / f"{query_name}.tar.gz"
        if summary_file.exists() and tar_file.exists():
            logger.info(f"Skipping completed query: {query_name}")
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
    sensitive_mode: bool,
    threads_per_worker: int,
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
            sensitive_mode=sensitive_mode,
            threads=threads_per_worker,
        )
        future_to_query[future] = query_file
    return future_to_query


def _collect_query_results(future_to_query: Dict[Future, Path], logger) -> List[Dict[str, Any]]:
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
    sensitive_mode: bool,
    threads_per_worker: int,
    logger,
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
            sensitive_mode,
            threads_per_worker,
        )
        results = _collect_query_results(future_to_query, logger)

    return results, total_queries


def _count_completed(results: List[Dict[str, Any]]) -> int:
    return len([result for result in results if result["status"] == "complete"])


def gvclass_flow(
    query_dir: str,
    output_dir: str,
    database_path: Optional[str] = None,
    total_threads: int = 16,
    max_workers: Optional[int] = None,
    threads_per_worker: Optional[int] = None,
    tree_method: str = "fasttree",
    mode_fast: bool = True,
    sensitive_mode: bool = False,
    genetic_codes: List[int] = [0, 1, 4, 6, 11, 15, 29, 106, 129],
    cluster_type: str = "local",
    cluster_config: Optional[Dict[str, Any]] = None,
    resume: bool = False,
):
    """Main GVClass pipeline."""
    logger = logging.getLogger("gvclass_prefect")
    logger.info("Validating inputs and setting up directories")
    config = validate_and_setup_task(query_dir, output_dir, database_path)
    _apply_resume_filter(config, resume, logger)

    n_queries = len(config["query_files"])
    logger.info(f"Found {n_queries} queries to process")
    n_workers, threads_per_worker = _resolve_worker_distribution(
        n_queries, total_threads, max_workers, threads_per_worker
    )
    logger.info(
        f"Parallelization strategy: {n_workers} workers × {threads_per_worker} threads"
    )

    logger.info("Starting parallel query processing")
    results, total_queries = _process_queries_in_parallel(
        config, n_workers, genetic_codes, tree_method, mode_fast, sensitive_mode, threads_per_worker, logger
    )
    completed_count = _count_completed(results)
    logger.info(
        f"All queries complete. Successfully processed {completed_count}/{total_queries} queries"
    )

    logger.info("Creating final summary...")
    summary_file = create_final_summary_task(results, config["output_path"])
    logger.info(f"Pipeline completed! Results in: {config['output_path']}")
    logger.info(f"Summary written to: {summary_file}")
    logger.info(f"Successfully processed {completed_count}/{len(results)} queries")

    _ = cluster_type
    _ = cluster_config
    return summary_file


# Deployment configurations for different environments
def create_local_deployment():
    """Create deployment for local execution."""
    from prefect.deployments import Deployment

    deployment = Deployment.build_from_flow(
        flow=gvclass_flow,
        name="gvclass-local",
        parameters={
            "cluster_type": "local",
            "total_threads": 16,
            "mode_fast": True,
            "sensitive_mode": False,
        },
        tags=["gvclass", "local"],
        description="GVClass pipeline for local execution",
    )

    return deployment



def create_slurm_deployment():
    """Create deployment for SLURM cluster execution."""
    from prefect.deployments import Deployment

    deployment = Deployment.build_from_flow(
        flow=gvclass_flow,
        name="gvclass-slurm",
        parameters={
            "cluster_type": "slurm",
            "total_threads": 128,
            "sensitive_mode": False,
            "cluster_config": {
                "queue": "normal",
                "project": "gvclass",
                "walltime": "08:00:00",
            },
        },
        tags=["gvclass", "hpc", "slurm"],
        description="GVClass pipeline for SLURM cluster execution",
    )

    return deployment
