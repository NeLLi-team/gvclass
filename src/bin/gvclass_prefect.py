#!/usr/bin/env python
"""
GVClass command-line interface with proper Prefect orchestration.
"""

import click
import sys
from pathlib import Path
import logging
import os
import tempfile

# Set default Prefect home early to avoid permission issues
# This will be overridden later with the actual output directory if possible
default_prefect_home = Path(tempfile.mkdtemp(prefix="prefect_", dir="/tmp"))
os.environ["PREFECT_HOME"] = str(default_prefect_home)

from src.utils import setup_logging

# Disable Prefect server requirement for local runs
os.environ.setdefault("PREFECT_API_URL", "")

# Configure logging before importing Prefect
logging.getLogger('prefect.events').setLevel(logging.ERROR)
logging.getLogger('prefect._internal').setLevel(logging.ERROR)
logging.getLogger('prefect').setLevel(logging.WARNING)

# Import after logging configuration to suppress Prefect warnings
from src.pipeline.prefect_flow import gvclass_flow  # noqa: E402


@click.command()
@click.argument('querydir', type=click.Path(exists=True))
@click.option('--output-dir', '-o', help='Output directory (default: querydir_results)')
@click.option('--database', '-d', envvar='GVCLASS_DB', default='resources',
              help='Path to GVClass database (or set GVCLASS_DB env var)')
@click.option('--threads', '-t', default=16, help='Total threads available')
@click.option('--max-workers', '-j', type=int, help='Maximum parallel workers (auto if not set)')
@click.option('--threads-per-worker', type=int, help='Threads per worker (auto if not set)')
@click.option('--tree-method', default='fasttree', type=click.Choice(['fasttree', 'iqtree']),
              help='Tree building method')
@click.option('--mode-fast/--no-mode-fast', default=True,
              help='Skip order-level marker trees when True')
@click.option('--cluster-type', default='local', 
              type=click.Choice(['local', 'slurm', 'pbs', 'sge']),
              help='Type of cluster to use')
@click.option('--cluster-queue', help='Queue/partition for cluster jobs')
@click.option('--cluster-project', help='Project/account for cluster billing')
@click.option('--cluster-walltime', default='04:00:00', help='Walltime for cluster jobs')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--resume', is_flag=True, help='Resume from previous run, skipping completed queries')
def main(querydir, output_dir, database, threads, max_workers, threads_per_worker,
        tree_method, mode_fast, cluster_type, cluster_queue, cluster_project,
        cluster_walltime, verbose, resume):
    """
    Run GVClass pipeline with proper Prefect+Dask orchestration.
    
    This implementation properly uses Prefect's @flow and @task decorators
    with dynamic DaskTaskRunner configuration for true parallel execution.
    
    QUERYDIR: Directory containing query sequences (.fna or .faa files)
    """
    # Setup logging
    log_level = "DEBUG" if verbose else "INFO"
    logger = setup_logging("gvclass_prefect", level=log_level)
    
    # Validate inputs
    query_path = Path(querydir).resolve()
    if not query_path.exists():
        click.echo(f"Error: Query directory {querydir} does not exist", err=True)
        sys.exit(1)
    
    # Set output directory
    if output_dir:
        output_path = Path(output_dir).resolve()
    else:
        output_path = Path.cwd() / f"{query_path.name}_results"
    
    # Database path
    db_path = database if database else None
    
    # Prepare cluster configuration
    cluster_config = {}
    if cluster_queue:
        cluster_config['queue'] = cluster_queue
    if cluster_project:
        cluster_config['project'] = cluster_project
    if cluster_walltime:
        cluster_config['walltime'] = cluster_walltime
    
    click.echo("GVClass Pipeline v1.1.0 (Proper Prefect+Dask)")
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
    click.echo("")
    
    # Set Prefect home to output directory to avoid conflicts in parallel runs
    # First ensure output directory exists and is writable
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Try to set Prefect home in output directory, fall back to temp if not writable
    prefect_home = output_path / ".prefect"
    try:
        prefect_home.mkdir(parents=True, exist_ok=True)
        os.environ["PREFECT_HOME"] = str(prefect_home)
    except PermissionError:
        # Fall back to temp directory
        import tempfile
        temp_prefect = Path(tempfile.mkdtemp(prefix="prefect_", dir="/tmp"))
        os.environ["PREFECT_HOME"] = str(temp_prefect)
        logger.debug(f"Using temporary Prefect home: {temp_prefect}")
    
    try:
        # Run the proper flow
        result = gvclass_flow(
            query_dir=str(query_path),
            output_dir=str(output_path),
            database_path=db_path,
            total_threads=threads,
            max_workers=max_workers,
            threads_per_worker=threads_per_worker,
            tree_method=tree_method,
            mode_fast=mode_fast,
            cluster_type=cluster_type,
            cluster_config=cluster_config if cluster_config else None,
            resume=resume,
        )
        
        click.echo("\nPipeline completed successfully!")
        click.echo(f"Results written to: {result}")
        
    except Exception as e:
        logger.exception("Pipeline failed")
        click.echo(f"Pipeline failed with error: {e}", err=True)
        sys.exit(1)




if __name__ == '__main__':
    main()