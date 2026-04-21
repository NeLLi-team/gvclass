#!/usr/bin/env python
"""
GVClass command-line interface.
"""

import click
import sys
from pathlib import Path
from src.utils import setup_logging

from src.pipeline.parallel_runner import gvclass_flow  # noqa: E402


def resolve_paths(querydir: str, output_dir: str | None) -> tuple[Path, Path]:
    query_path = Path(querydir).resolve()
    if not query_path.exists():
        raise FileNotFoundError(f"Query directory {querydir} does not exist")
    if output_dir:
        return query_path, Path(output_dir).resolve()
    return query_path, Path.cwd() / f"{query_path.name}_results"


def build_cluster_config(
    cluster_queue: str | None, cluster_project: str | None, cluster_walltime: str | None
) -> dict:
    config = {}
    if cluster_queue:
        config["queue"] = cluster_queue
    if cluster_project:
        config["project"] = cluster_project
    if cluster_walltime:
        config["walltime"] = cluster_walltime
    return config


def print_run_configuration(
    query_path: Path,
    output_path: Path,
    db_path: str | None,
    threads: int,
    max_workers: int | None,
    threads_per_worker: int | None,
    cluster_type: str,
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive: bool,
) -> None:
    from src.__version__ import __version__ as _gvclass_version

    click.echo(f"GVClass Pipeline v{_gvclass_version}")
    click.echo(f"Query directory: {query_path}")
    click.echo(f"Output directory: {output_path}")
    click.echo(f"Database: {db_path if db_path else 'Will download/use default'}")
    click.echo(f"Total threads: {threads}")
    if max_workers:
        click.echo(f"Max workers: {max_workers}")
    if threads_per_worker:
        click.echo(f"Threads per worker: {threads_per_worker}")
    else:
        click.echo("Worker distribution: Auto-calculated")
    click.echo(f"Cluster type: {cluster_type}")
    click.echo(f"Tree method: {tree_method}")
    click.echo(f"Fast mode: {mode_fast}")
    click.echo(f"Completeness mode: {completeness_mode}")
    click.echo(f"Sensitive mode: {sensitive}")
    click.echo("")


def run_flow(
    query_path: Path,
    output_path: Path,
    db_path: str | None,
    threads: int,
    max_workers: int | None,
    threads_per_worker: int | None,
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive: bool,
    cluster_type: str,
    cluster_config: dict,
    resume: bool,
    allow_short: bool = False,
):
    return gvclass_flow(
        query_dir=str(query_path),
        output_dir=str(output_path),
        database_path=db_path,
        total_threads=threads,
        max_workers=max_workers,
        threads_per_worker=threads_per_worker,
        tree_method=tree_method,
        mode_fast=mode_fast,
        completeness_mode=completeness_mode,
        sensitive_mode=sensitive,
        cluster_type=cluster_type,
        cluster_config=cluster_config if cluster_config else None,
        resume=resume,
        allow_short=allow_short,
    )


def prepare_cli_context(
    querydir: str,
    output_dir: str | None,
    database: str | None,
    cluster_queue: str | None,
    cluster_project: str | None,
    cluster_walltime: str | None,
    threads: int,
    max_workers: int | None,
    threads_per_worker: int | None,
    cluster_type: str,
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive: bool,
) -> tuple[Path, Path, str | None, dict]:
    query_path, output_path = resolve_paths(querydir, output_dir)
    db_path = database if database else None
    cluster_config = build_cluster_config(
        cluster_queue=cluster_queue,
        cluster_project=cluster_project,
        cluster_walltime=cluster_walltime,
    )
    print_run_configuration(
        query_path=query_path,
        output_path=output_path,
        db_path=db_path,
        threads=threads,
        max_workers=max_workers,
        threads_per_worker=threads_per_worker,
        cluster_type=cluster_type,
        tree_method=tree_method,
        mode_fast=mode_fast,
        completeness_mode=completeness_mode,
        sensitive=sensitive,
    )
    return query_path, output_path, db_path, cluster_config


def execute_cli_flow(
    query_path: Path,
    output_path: Path,
    db_path: str | None,
    threads: int,
    max_workers: int | None,
    threads_per_worker: int | None,
    tree_method: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive: bool,
    cluster_type: str,
    cluster_config: dict,
    resume: bool,
    allow_short: bool = False,
):
    result = run_flow(
        query_path=query_path,
        output_path=output_path,
        db_path=db_path,
        threads=threads,
        max_workers=max_workers,
        threads_per_worker=threads_per_worker,
        tree_method=tree_method,
        mode_fast=mode_fast,
        completeness_mode=completeness_mode,
        sensitive=sensitive,
        cluster_type=cluster_type,
        cluster_config=cluster_config,
        resume=resume,
        allow_short=allow_short,
    )
    click.echo("\nPipeline completed successfully!")
    click.echo(f"Results written to: {result}")


@click.command()
@click.argument("querydir", type=click.Path(exists=True))
@click.option("--output-dir", "-o", help="Output directory (default: querydir_results)")
@click.option(
    "--database",
    "-d",
    envvar="GVCLASS_DB",
    default="resources",
    help="Path to GVClass database (or set GVCLASS_DB env var)",
)
@click.option("--threads", "-t", default=16, help="Total threads available")
@click.option(
    "--max-workers", "-j", type=int, help="Maximum parallel workers (auto if not set)"
)
@click.option(
    "--threads-per-worker", type=int, help="Threads per worker (auto if not set)"
)
@click.option(
    "--tree-method",
    default="fasttree",
    type=click.Choice(["fasttree", "iqtree"]),
    help="Tree building method",
)
@click.option(
    "--mode-fast/--extended",
    "-f/-e",
    default=True,
    help="Fast mode (default) or extended mode with all marker trees",
)
@click.option(
    "--completeness-mode",
    default="novelty-aware",
    type=click.Choice(["legacy", "novelty-aware"]),
    show_default=True,
    help="Completeness estimator to surface as the primary estimate",
)
@click.option(
    "--sensitive",
    is_flag=True,
    help="Sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
)
@click.option(
    "--cluster-type",
    default="local",
    type=click.Choice(["local", "slurm", "pbs", "sge"]),
    help="Type of cluster to use",
)
@click.option("--cluster-queue", help="Queue/partition for cluster jobs")
@click.option("--cluster-project", help="Project/account for cluster billing")
@click.option(
    "--cluster-walltime", default="04:00:00", help="Walltime for cluster jobs"
)
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
@click.option(
    "--resume",
    is_flag=True,
    help="Resume from previous run, skipping completed queries",
)
@click.option(
    "--allow-short",
    is_flag=True,
    help=(
        "Accept input FNA files shorter than the 20 kb minimum. Use for "
        "per-contig runs (often paired with --contigs) or for exploratory "
        "runs on short inputs."
    ),
)
def main(querydir, output_dir, database, threads, max_workers, threads_per_worker, tree_method, mode_fast, completeness_mode, sensitive, cluster_type, cluster_queue, cluster_project, cluster_walltime, verbose, resume, allow_short):
    log_level = "DEBUG" if verbose else "INFO"
    logger = setup_logging("gvclass_runner", level=log_level)

    try:
        query_path, output_path, db_path, cluster_config = prepare_cli_context(
            querydir=querydir,
            output_dir=output_dir,
            database=database,
            cluster_queue=cluster_queue,
            cluster_project=cluster_project,
            cluster_walltime=cluster_walltime,
            threads=threads,
            max_workers=max_workers,
            threads_per_worker=threads_per_worker,
            cluster_type=cluster_type,
            tree_method=tree_method,
            mode_fast=mode_fast,
            completeness_mode=completeness_mode,
            sensitive=sensitive,
        )
        execute_cli_flow(query_path, output_path, db_path, threads, max_workers, threads_per_worker, tree_method, mode_fast, completeness_mode, sensitive, cluster_type, cluster_config, resume, allow_short=allow_short)
    except FileNotFoundError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)
    except Exception as e:
        logger.exception("Pipeline failed")
        click.echo(f"Pipeline failed with error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
