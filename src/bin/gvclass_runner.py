#!/usr/bin/env python
"""
GVClass command-line interface.
"""

import sys
from pathlib import Path

import click

from src.pipeline.parallel_runner import gvclass_flow
from src.pipeline.run_config import RunConfig
from src.utils import setup_logging


def resolve_paths(querydir: str, output_dir: str | None) -> tuple[Path, Path]:
    query_path = Path(querydir).resolve()
    if not query_path.exists():
        raise FileNotFoundError(f"Query directory {querydir} does not exist")
    if output_dir:
        return query_path, Path(output_dir).resolve()
    return query_path, Path.cwd() / f"{query_path.name}_results"


def print_run_configuration(config: RunConfig) -> None:
    """Print the resolved settings of the run."""
    from src.__version__ import __version__ as _gvclass_version

    click.echo(f"GVClass Pipeline v{_gvclass_version}")
    click.echo(f"Query directory: {config.query_path}")
    click.echo(f"Output directory: {config.output_dir}")
    click.echo(f"Database: {config.database_path}")
    click.echo(f"Total threads: {config.threads}")
    if config.max_workers:
        click.echo(f"Max workers: {config.max_workers}")
    if config.threads_per_worker:
        click.echo(f"Threads per worker: {config.threads_per_worker}")
    else:
        click.echo("Worker distribution: Auto-calculated")
    click.echo(f"Tree method: {config.tree_method}")
    click.echo(f"IQ-TREE mode: {config.iqtree_mode}")
    click.echo(f"Fast mode: {config.mode_fast}")
    click.echo(f"Completeness mode: {config.completeness_mode}")
    click.echo(f"Sensitive mode: {config.sensitive_mode}")
    if config.min_sequence_length is not None:
        click.echo(f"Minimum FNA length: {config.min_sequence_length} bp")
    click.echo("")


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
    default="veryfasttree",
    type=click.Choice(["veryfasttree", "iqtree", "fasttree"]),
    help="Tree building method (default veryfasttree; iqtree is slower/more accurate)",
)
@click.option(
    "--iqtree-mode",
    default="fast",
    type=click.Choice(["fast", "ufboot"]),
    help="IQ-TREE species-tree search: fast (--fast) or ufboot (ultrafast bootstrap)",
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
    "--sensitive/--no-sensitive",
    default=False,
    help="Sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
)
@click.option("--cluster-type", hidden=True, help="Has no effect.")
@click.option("--cluster-queue", hidden=True, help="Has no effect.")
@click.option("--cluster-project", hidden=True, help="Has no effect.")
@click.option("--cluster-walltime", hidden=True, help="Has no effect.")
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
        "Bypass the nucleotide length floor. Use for exploratory runs on "
        "short inputs."
    ),
)
@click.option(
    "--min-length",
    type=int,
    default=None,
    help=(
        "Minimum total nucleotide length (bp) for each input FNA. "
        "Use 0 to disable the length floor without disabling other validation."
    ),
)
@click.option(
    "--species-tree",
    is_flag=True,
    help=(
        "Build one supermatrix species tree per query (NCLDV GVOG8 by default; "
        "PPV/MIRUS routed by taxonomy) and assign tree-placement taxonomy "
        "(writes out/species_tree/<query>/ and the summary columns)."
    ),
)
@click.option(
    "--species-tree-combined",
    is_flag=True,
    help=(
        "Implies --species-tree; additionally build one combined species tree "
        "over all queries at once (writes out/species_tree/combined.*)."
    ),
)
@click.option(
    "--species-tree-trim",
    type=click.Choice(["witchi", "pytrimal", "none"]),
    default="witchi",
    show_default=True,
    help="Supermatrix column trimming for --species-tree: witchi (rigorous), pytrimal (fast), or none.",
)
@click.option(
    "--screen/--no-screen",
    "contig_screen",
    default=True,
    help="Screen each nucleotide bin for contaminating contigs (default on).",
)
@click.option(
    "--contigs",
    "contigs_mode",
    is_flag=True,
    help="The queries are single contigs from larger inputs; the screen skips them.",
)
def main(  # noqa: PLR0913  # Click supplies one keyword per declared option.
    *,
    querydir: str,
    output_dir: str | None,
    database: str,
    threads: int,
    max_workers: int | None,
    threads_per_worker: int | None,
    tree_method: str,
    iqtree_mode: str,
    mode_fast: bool,
    completeness_mode: str,
    sensitive: bool,
    verbose: bool,
    resume: bool,
    allow_short: bool,
    min_length: int | None,
    species_tree: bool,
    species_tree_combined: bool,
    species_tree_trim: str,
    contig_screen: bool,
    contigs_mode: bool,
    **_ignored_cluster_options: str | None,
) -> None:
    """Run the pipeline with resolved command-line settings."""
    log_level = "DEBUG" if verbose else "INFO"
    logger = setup_logging("gvclass_runner", level=log_level)

    try:
        query_path, output_path = resolve_paths(querydir, output_dir)
        config = RunConfig(
            query_path=query_path,
            output_dir=output_path,
            database_path=Path(database),
            threads=threads,
            max_workers=max_workers,
            threads_per_worker=threads_per_worker,
            tree_method=tree_method,
            iqtree_mode=iqtree_mode,
            mode_fast=mode_fast,
            completeness_mode=completeness_mode,
            sensitive_mode=sensitive,
            resume=resume,
            allow_short=allow_short,
            min_sequence_length=min_length,
            species_tree=species_tree,
            species_tree_combined=species_tree_combined,
            species_tree_trim=species_tree_trim,
            contig_screen=contig_screen,
            contigs_mode=contigs_mode,
            contigs_min_length=0,  # The CLI splits contigs before the runner starts.
        )
        print_run_configuration(config)
        summary_file = gvclass_flow(config)
        click.echo("\nPipeline completed successfully!")
        click.echo(f"Results written to: {summary_file}")
    except FileNotFoundError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)
    except Exception as e:
        logger.exception("Pipeline failed")
        click.echo(f"Pipeline failed with error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
