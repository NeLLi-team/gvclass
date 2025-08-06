#!/usr/bin/env python
"""
Run GVClass pipeline with clean output and resource monitoring.
"""

import sys
import os
from pathlib import Path
import subprocess
import argparse
from datetime import datetime
import threading
import time
import yaml
import psutil


def count_files_in_directory(directory):
    """Count .fna and .faa files in a directory."""
    path = Path(directory)
    if not path.exists():
        return 0
    fna_files = list(path.glob("*.fna"))
    faa_files = list(path.glob("*.faa"))
    return len(fna_files) + len(faa_files)


def load_config(config_file="config/gvclass_config.yaml"):
    """Load configuration from YAML file."""
    config_path = Path(config_file)

    # Try different locations for config file
    if not config_path.exists():
        # Try in current directory
        config_path = Path.cwd() / config_file

    if not config_path.exists():
        # Try relative to script
        config_path = Path(__file__).parent / config_file

    if not config_path.exists():
        # Try old location
        config_path = Path.cwd() / "gvclass_config.yaml"

    if config_path.exists():
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    # Return default config if file not found
    return {
        "database": {
            "path": "resources",
            "download_url": "https://portal.nersc.gov/cfs/m342/GVClass/gvclass_db_v1.1.tar.gz",
            "expected_size": 701,
        },
        "pipeline": {
            "tree_method": "fasttree",
            "mode_fast": True,
            "threads": 16,
            "output_pattern": "{query_dir}_results",
        },
    }


def check_and_setup_database(db_path, config):
    """Check if database exists, download if needed."""

    db_path = Path(db_path)

    # Check for key database files
    required_checks = [
        db_path / "models" / "combined.hmm",
        db_path / "database",
        db_path / "gvclassJuly25_labels.tsv",
    ]

    if db_path.exists() and any(check.exists() for check in required_checks):
        print(f"✅ Database found at: {db_path}")
        return True

    print(f"\n📦 Database not found at {db_path}")
    print("📥 Setting up database...")

    try:
        # Use DatabaseManager to setup the database
        from src.utils.database_manager import DatabaseManager

        DatabaseManager.setup_database(str(db_path))
        print("✅ Database setup complete!")
        return True

    except Exception as e:
        print(f"\n❌ Failed to setup database: {e}")
        print("   You can try manually with: pixi run setup-db")
        return False


class ResourceMonitor:
    """Monitor progress and resource usage."""

    def __init__(self, output_dir, total_queries, process_pid, initial_completed=0):
        self.output_dir = Path(output_dir)
        self.total_queries = total_queries
        self.process_pid = process_pid
        self.current_query = 0
        self.initial_completed = initial_completed  # For resume mode
        self.current_task = "Starting"
        self.running = True
        self.peak_memory_mb = 0
        self.current_memory_mb = 0

    def monitor(self):
        """Monitor progress and resources."""
        last_output = ""
        try:
            process = psutil.Process(self.process_pid)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            process = None

        while self.running:
            # Count completed queries by looking for summary.tab files in output root
            if self.output_dir.exists():
                total_completed = len(list(self.output_dir.glob("*.summary.tab")))
                # Subtract initial completed to get queries completed in this run
                self.current_query = total_completed - self.initial_completed

            # Monitor memory usage
            if process:
                try:
                    mem_info = process.memory_info()
                    self.current_memory_mb = mem_info.rss / 1024 / 1024
                    self.peak_memory_mb = max(
                        self.peak_memory_mb, self.current_memory_mb
                    )

                    # Include children processes
                    for child in process.children(recursive=True):
                        try:
                            child_mem = child.memory_info().rss / 1024 / 1024
                            self.current_memory_mb += child_mem
                            self.peak_memory_mb = max(
                                self.peak_memory_mb, self.current_memory_mb
                            )
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

            # Display progress
            if self.total_queries > 0:
                percent = int((self.current_query / self.total_queries) * 100)
                output = f"Progress: [{percent:3d}%] Query {self.current_query}/{self.total_queries} | Memory: {self.current_memory_mb:.0f}MB (peak: {self.peak_memory_mb:.0f}MB)"
                if output != last_output:
                    print(f"\r{output:<80}", end="", flush=True)
                    last_output = output

            time.sleep(1)

    def stop(self):
        """Stop monitoring."""
        self.running = False


def main():
    """Run the pipeline with clean output and resource monitoring."""
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Run GVClass pipeline with clean output"
    )
    parser.add_argument(
        "query_dir", nargs="?", help="Query directory containing .fna or .faa files"
    )
    parser.add_argument(
        "-o", "--output-dir", help="Output directory (default: <query_dir>_results)"
    )
    parser.add_argument(
        "-c",
        "--config",
        default="config/gvclass_config.yaml",
        help="Configuration file (default: config/gvclass_config.yaml)",
    )
    parser.add_argument("-d", "--database", help="Database path (overrides config)")
    parser.add_argument(
        "-t", "--threads", type=int, help="Number of threads (overrides config)"
    )
    parser.add_argument(
        "-j", "--max-workers", type=int, help="Maximum parallel workers"
    )
    parser.add_argument("--threads-per-worker", type=int, help="Threads per worker")
    parser.add_argument(
        "--tree-method",
        choices=["fasttree", "iqtree"],
        help="Tree building method (overrides config)",
    )
    parser.add_argument(
        "--mode-fast", action="store_true", help="Enable fast mode (overrides config)"
    )
    parser.add_argument(
        "--no-mode-fast",
        action="store_true",
        help="Disable fast mode (overrides config)",
    )
    parser.add_argument(
        "--cluster-type",
        default="local",
        choices=["local", "slurm", "pbs", "sge"],
        help="Type of cluster to use",
    )
    parser.add_argument("--cluster-queue", help="Queue/partition for cluster jobs")
    parser.add_argument("--cluster-project", help="Project/account for cluster billing")
    parser.add_argument(
        "--cluster-walltime", default="04:00:00", help="Walltime for cluster jobs"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version information and exit",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous run, skipping completed queries",
    )

    args = parser.parse_args()

    # Handle --version flag
    if args.version:
        from src.utils.database_manager import DatabaseManager

        print("GVClass Pipeline")
        print("  Software version: v1.1.0")

        # Load config to get database path
        config = load_config(args.config)
        db_path = Path(args.database if args.database else config["database"]["path"])

        if db_path.exists():
            db_version = DatabaseManager.get_database_version(db_path)
            print(f"  Database version: {db_version}")
        else:
            print("  Database version: not installed")

        sys.exit(0)

    # Check if query_dir is provided for non-version commands
    if not args.query_dir:
        parser.error("query_dir is required unless using --version")

    # Load configuration
    config = load_config(args.config)

    # Apply configuration with command-line overrides
    database = args.database if args.database else config["database"]["path"]
    threads = args.threads if args.threads else config["pipeline"].get("threads", 16)
    tree_method = (
        args.tree_method
        if args.tree_method
        else config["pipeline"].get("tree_method", "fasttree")
    )

    # Handle mode_fast with explicit flags taking priority
    if args.no_mode_fast:
        mode_fast = False
    elif args.mode_fast:
        mode_fast = True
    else:
        mode_fast = config["pipeline"].get("mode_fast", True)

    # Handle output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        query_base = Path(args.query_dir).name
        output_pattern = config["pipeline"].get("output_pattern", "{query_dir}_results")
        output_dir = output_pattern.format(query_dir=query_base)

    # Check and setup database
    if not check_and_setup_database(database, config):
        sys.exit(1)

    # Count input files
    n_queries = count_files_in_directory(args.query_dir)
    n_skipped = 0

    if n_queries == 0:
        print(
            f"❌ ERROR: No .fna or .faa files found in {args.query_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Create output directory if it doesn't exist (needed for resume check)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # If resume mode, count already completed queries
    if args.resume:
        output_path = Path(output_dir)
        # Count queries that have both summary.tab and tar.gz
        # Need to remove the full extension for both file types
        summary_files = set(
            f.name.replace(".summary.tab", "")
            for f in output_path.glob("*.summary.tab")
        )
        tar_files = set(
            f.name.replace(".tar.gz", "") for f in output_path.glob("*.tar.gz")
        )
        n_skipped = len(summary_files & tar_files)  # Intersection

    # Calculate workers and threads per worker
    # If both -j and --threads-per-worker are specified, use them directly
    if args.max_workers and args.threads_per_worker:
        n_workers = args.max_workers
        threads_per_worker = args.threads_per_worker
    # If only -j is specified, use it with remaining threads distributed
    elif args.max_workers:
        n_workers = args.max_workers
        threads_per_worker = max(1, threads // n_workers)
    # If only --threads-per-worker is specified, calculate workers
    elif args.threads_per_worker:
        threads_per_worker = args.threads_per_worker
        n_workers = min(n_queries, max(1, threads // threads_per_worker))
    # If neither is specified, use heuristics
    else:
        if n_queries <= 4 and threads >= n_queries * 2:
            n_workers = n_queries
            threads_per_worker = threads // n_workers
        else:
            # Default to 4 threads per worker as a good balance
            threads_per_worker = 4
            n_workers = min(n_queries, max(1, threads // threads_per_worker))

    # Print clean header with ASCII art and color gradient
    print("\n")
    # ANSI color codes for gradient: purple -> pink -> red -> orange -> yellow
    colors = [
        "\033[95m",  # Purple (G)
        "\033[95m",  # Purple (V)
        "\033[91m",  # Pink (C)
        "\033[31m",  # Red (L)
        "\033[31m",  # Red (A)
        "\033[38;5;208m",  # Orange (S) - 256 color
        "\033[93m",  # Yellow (S)
    ]
    reset = "\033[0m"

    print("    ╔═══════════════════════════════════════════════════════════╗")
    print("    ║                                                           ║")
    print(
        f"    ║  {colors[0]}██████╗{reset} {colors[1]}██╗   ██╗{reset}{colors[2]} ██████╗{reset}{colors[3]}██╗{reset}      {colors[4]}█████╗{reset} {colors[5]}███████╗{reset}{colors[6]}███████╗{reset} ║"
    )
    print(
        f"    ║ {colors[0]}██╔════╝{reset} {colors[1]}██║   ██║{reset}{colors[2]}██╔════╝{reset}{colors[3]}██║{reset}     {colors[4]}██╔══██╗{reset}{colors[5]}██╔════╝{reset}{colors[6]}██╔════╝{reset} ║"
    )
    print(
        f"    ║ {colors[0]}██║  ███╗{reset}{colors[1]}██║   ██║{reset}{colors[2]}██║{reset}     {colors[3]}██║{reset}     {colors[4]}███████║{reset}{colors[5]}███████╗{reset}{colors[6]}███████╗{reset} ║"
    )
    print(
        f"    ║ {colors[0]}██║   ██║{reset}{colors[1]}╚██╗ ██╔╝{reset}{colors[2]}██║{reset}     {colors[3]}██║{reset}     {colors[4]}██╔══██║{reset}{colors[5]}╚════██║{reset}{colors[6]}╚════██║{reset} ║"
    )
    print(
        f"    ║ {colors[0]}╚██████╔╝{reset} {colors[1]}╚████╔╝{reset} {colors[2]}╚██████╗{reset}{colors[3]}███████╗{reset}{colors[4]}██║  ██║{reset}{colors[5]}███████║{reset}{colors[6]}███████║{reset} ║"
    )
    print(
        f"    ║  {colors[0]}╚═════╝{reset}   {colors[1]}╚═══╝{reset}   {colors[2]}╚═════╝{reset}{colors[3]}╚══════╝{reset}{colors[4]}╚═╝  ╚═╝{reset}{colors[5]}╚══════╝{reset}{colors[6]}╚══════╝{reset} ║"
    )
    print("    ║                                                           ║")
    print("    ║         Giant Virus Classification Tool v1.1.0            ║")
    print("    ║                                                           ║")
    print("    ╚═══════════════════════════════════════════════════════════╝")
    print()
    print("=" * 60)
    print("                    Pipeline Configuration")
    print("=" * 60)
    print(f"Config file: {args.config}")
    print(f"Query directory: {args.query_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Database: {database}")
    print(f"Threads: {threads} (Workers: {n_workers} × {threads_per_worker} threads)")
    print(f"Tree method: {tree_method}")
    print(f"Fast mode: {mode_fast}")
    if args.resume:
        print("Resume mode: ENABLED")
        if n_skipped > 0:
            print(f"  - Skipping {n_skipped} completed queries")
            print(f"  - Processing {n_queries - n_skipped} remaining queries")
        else:
            print("  - No completed queries found")
    else:
        print("Resume mode: DISABLED")
    print(f"Total queries: {n_queries}")
    print("=" * 60)

    start_time = datetime.now()

    # Build command
    cmd = [
        sys.executable,
        "-m",
        "src.bin.gvclass_prefect",
        args.query_dir,
        "-o",
        output_dir,
        "-d",
        database,
        "-t",
        str(threads),
        "-j",
        str(n_workers),
        "--threads-per-worker",
        str(threads_per_worker),
        "--tree-method",
        tree_method,
        "--cluster-type",
        args.cluster_type,
    ]

    if mode_fast:
        cmd.append("--mode-fast")
    else:
        cmd.append("--no-mode-fast")

    if args.cluster_queue:
        cmd.extend(["--cluster-queue", args.cluster_queue])
    if args.cluster_project:
        cmd.extend(["--cluster-project", args.cluster_project])
    if args.cluster_walltime:
        cmd.extend(["--cluster-walltime", args.cluster_walltime])
    if args.verbose:
        cmd.append("--verbose")
    if args.resume:
        cmd.append("--resume")

    # Set environment
    env = os.environ.copy()
    env["PREFECT_API_URL"] = ""
    env["PREFECT_LOGGING_LEVEL"] = "ERROR"
    env["PYTHONWARNINGS"] = "ignore"

    # Start process
    if args.verbose:
        process = subprocess.Popen(cmd, env=env)
    else:
        process = subprocess.Popen(
            cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )

    # Start resource monitor with actual queries to process
    queries_to_process = n_queries - n_skipped if args.resume else n_queries
    monitor = ResourceMonitor(
        output_dir, queries_to_process, process.pid, initial_completed=n_skipped
    )
    monitor_thread = threading.Thread(target=monitor.monitor)
    monitor_thread.start()

    try:
        # Wait for process
        return_code = process.wait()

        # Stop monitor
        monitor.stop()
        monitor_thread.join()
        print("\r" + " " * 80 + "\r", end="")  # Clear line

        if return_code == 0:
            elapsed = datetime.now() - start_time
            minutes = int(elapsed.total_seconds() / 60)
            seconds = int(elapsed.total_seconds() % 60)

            print("\n✅ Pipeline completed successfully!")

            # Validate outputs before generating summary
            print("🔍 Validating pipeline outputs...")
            from src.utils.output_validator import validate_pipeline_outputs

            validation_results = validate_pipeline_outputs(
                output_dir=Path(output_dir), verbose=args.verbose
            )

            if not validation_results["success"]:
                print(
                    f"⚠️  Validation found {len(validation_results['issues'])} issue(s):"
                )
                for issue in validation_results["issues"][:5]:  # Show first 5 issues
                    print(f"  - {issue}")
                if len(validation_results["issues"]) > 5:
                    print(f"  ... and {len(validation_results['issues']) - 5} more")
            else:
                print("✅ All outputs validated successfully")

            # Generate combined summary from individual summary.tab files
            print("📊 Generating combined summary...")
            output_path = Path(output_dir)
            summary_files = list(output_path.glob("*.summary.tab"))

            if summary_files:
                # Read all individual summaries and combine them
                combined_data = []
                headers = None

                for summary_file in sorted(summary_files):
                    with open(summary_file, "r") as f:
                        lines = f.readlines()
                        if lines:
                            if headers is None and len(lines) > 0:
                                headers = lines[0].strip()
                            if len(lines) > 1:  # Has data row
                                combined_data.append(lines[1].strip())

                # Write combined summary
                if combined_data and headers:
                    combined_summary = output_path / "gvclass_summary.tsv"
                    with open(combined_summary, "w") as f:
                        f.write(headers + "\n")
                        for data_line in combined_data:
                            f.write(data_line + "\n")
                    print(f"✅ Combined summary written to: {combined_summary}")
                else:
                    print("⚠️  No summary data found to combine")
            else:
                print("⚠️  No individual summary files found")

            print(f"⏱️  Runtime: {minutes}m {seconds}s")
            print(f"💾 Peak memory usage: {monitor.peak_memory_mb:.0f} MB")
            print(f"📊 Results saved to: {output_dir}")

            # Final validation summary
            if validation_results["success"]:
                print(
                    f"✅ Output validation: PASSED ({validation_results['queries_validated']} queries)"
                )
            else:
                print(
                    f"⚠️  Output validation: {validation_results['queries_failed']} queries with issues"
                )
        else:
            print(f"\n❌ Pipeline failed with exit code: {return_code}")
            sys.exit(1)

    except KeyboardInterrupt:
        monitor.stop()
        monitor_thread.join()
        print("\r" + " " * 80 + "\r", end="")
        print("\n⚠️  Pipeline interrupted by user")
        process.terminate()
        sys.exit(1)
    except Exception as e:
        monitor.stop()
        monitor_thread.join()
        print("\r" + " " * 80 + "\r", end="")
        print(f"\n❌ Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
