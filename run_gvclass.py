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
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    # Return default config if file not found
    return {
        'database': {
            'path': 'resources',
            'download_url': 'https://portal.nersc.gov/cfs/m342/GVClass/gvclass_db_v1.1.tar.gz',
            'expected_size': 701
        },
        'pipeline': {
            'tree_method': 'fasttree',
            'mode_fast': True,
            'threads': 16,
            'output_pattern': '{query_dir}_results'
        }
    }


def check_and_setup_database(db_path, config):
    """Check if database exists, download if needed."""
    
    db_path = Path(db_path)
    
    # Check for key database files
    required_checks = [
        db_path / 'models' / 'combined.hmm',
        db_path / 'database',
        db_path / 'ncldvApril24_labels.txt'
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
    def __init__(self, output_dir, total_queries, process_pid):
        self.output_dir = Path(output_dir)
        self.total_queries = total_queries
        self.process_pid = process_pid
        self.current_query = 0
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
            # Count completed queries
            if self.output_dir.exists():
                completed = len(list(self.output_dir.glob("*/stats/*.stats.tab")))
                if completed > self.current_query:
                    self.current_query = completed
            
            # Monitor memory usage
            if process:
                try:
                    mem_info = process.memory_info()
                    self.current_memory_mb = mem_info.rss / 1024 / 1024
                    self.peak_memory_mb = max(self.peak_memory_mb, self.current_memory_mb)
                    
                    # Include children processes
                    for child in process.children(recursive=True):
                        try:
                            child_mem = child.memory_info().rss / 1024 / 1024
                            self.current_memory_mb += child_mem
                            self.peak_memory_mb = max(self.peak_memory_mb, self.current_memory_mb)
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
    parser = argparse.ArgumentParser(description='Run GVClass pipeline with clean output')
    parser.add_argument('query_dir', help='Query directory containing .fna or .faa files')
    parser.add_argument('-o', '--output-dir', help='Output directory (default: <query_dir>_results)')
    parser.add_argument('-c', '--config', default='config/gvclass_config.yaml', 
                       help='Configuration file (default: config/gvclass_config.yaml)')
    parser.add_argument('-d', '--database', help='Database path (overrides config)')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads (overrides config)')
    parser.add_argument('-j', '--max-workers', type=int, help='Maximum parallel workers')
    parser.add_argument('--threads-per-worker', type=int, help='Threads per worker')
    parser.add_argument('--tree-method', choices=['fasttree', 'iqtree'], 
                       help='Tree building method (overrides config)')
    parser.add_argument('--mode-fast', action='store_true', help='Enable fast mode (overrides config)')
    parser.add_argument('--no-mode-fast', action='store_true', help='Disable fast mode (overrides config)')
    parser.add_argument('--cluster-type', default='local', 
                       choices=['local', 'slurm', 'pbs', 'sge'],
                       help='Type of cluster to use')
    parser.add_argument('--cluster-queue', help='Queue/partition for cluster jobs')
    parser.add_argument('--cluster-project', help='Project/account for cluster billing')
    parser.add_argument('--cluster-walltime', default='04:00:00', help='Walltime for cluster jobs')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Apply configuration with command-line overrides
    database = args.database if args.database else config['database']['path']
    threads = args.threads if args.threads else config['pipeline'].get('threads', 16)
    tree_method = args.tree_method if args.tree_method else config['pipeline'].get('tree_method', 'fasttree')
    
    # Handle mode_fast with explicit flags taking priority
    if args.no_mode_fast:
        mode_fast = False
    elif args.mode_fast:
        mode_fast = True
    else:
        mode_fast = config['pipeline'].get('mode_fast', True)
    
    # Handle output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        query_base = Path(args.query_dir).name
        output_pattern = config['pipeline'].get('output_pattern', '{query_dir}_results')
        output_dir = output_pattern.format(query_dir=query_base)
    
    # Check and setup database
    if not check_and_setup_database(database, config):
        sys.exit(1)
    
    # Count input files
    n_queries = count_files_in_directory(args.query_dir)
    
    # Calculate workers if not specified
    if args.max_workers:
        n_workers = args.max_workers
    else:
        if n_queries <= 4 and threads >= n_queries * 2:
            n_workers = n_queries
        else:
            n_workers = min(n_queries, max(1, threads // 4))
    
    if args.threads_per_worker:
        threads_per_worker = args.threads_per_worker
    else:
        threads_per_worker = max(1, threads // n_workers) if n_workers else threads
    
    # Print clean header with ASCII art and color gradient
    print("\n")
    # ANSI color codes for gradient: purple -> pink -> red -> orange -> yellow
    colors = [
        '\033[95m',  # Purple (G)
        '\033[95m',  # Purple (V)
        '\033[91m',  # Light red/pink (C)
        '\033[91m',  # Red (L)
        '\033[38;5;208m',  # Orange (A) - 256 color
        '\033[38;5;208m',  # Orange (S) - 256 color
        '\033[93m'   # Yellow (S)
    ]
    reset = '\033[0m'
    
    print("    ╔═══════════════════════════════════════════════════════════╗")
    print("    ║                                                           ║")
    print(f"    ║  {colors[0]}██████╗{reset} {colors[1]}██╗   ██╗{reset}{colors[2]} ██████╗{reset}{colors[3]}██╗{reset}      {colors[4]}█████╗{reset} {colors[5]}███████╗{reset}{colors[6]}███████╗{reset} ║")
    print(f"    ║ {colors[0]}██╔════╝{reset} {colors[1]}██║   ██║{reset}{colors[2]}██╔════╝{reset}{colors[3]}██║{reset}     {colors[4]}██╔══██╗{reset}{colors[5]}██╔════╝{reset}{colors[6]}██╔════╝{reset} ║")
    print(f"    ║ {colors[0]}██║  ███╗{reset}{colors[1]}██║   ██║{reset}{colors[2]}██║{reset}     {colors[3]}██║{reset}     {colors[4]}███████║{reset}{colors[5]}███████╗{reset}{colors[6]}███████╗{reset} ║")
    print(f"    ║ {colors[0]}██║   ██║{reset}{colors[1]}╚██╗ ██╔╝{reset}{colors[2]}██║{reset}     {colors[3]}██║{reset}     {colors[4]}██╔══██║{reset}{colors[5]}╚════██║{reset}{colors[6]}╚════██║{reset} ║")
    print(f"    ║ {colors[0]}╚██████╔╝{reset} {colors[1]}╚████╔╝{reset} {colors[2]}╚██████╗{reset}{colors[3]}███████╗{reset}{colors[4]}██║  ██║{reset}{colors[5]}███████║{reset}{colors[6]}███████║{reset} ║")
    print(f"    ║  {colors[0]}╚═════╝{reset}   {colors[1]}╚═══╝{reset}   {colors[2]}╚═════╝{reset}{colors[3]}╚══════╝{reset}{colors[4]}╚═╝  ╚═╝{reset}{colors[5]}╚══════╝{reset}{colors[6]}╚══════╝{reset} ║")
    print("    ║                                                           ║")
    print("    ║         Giant Virus Classification Tool v1.1.0            ║")
    print("    ║                                                           ║")
    print("    ╚═══════════════════════════════════════════════════════════╝")
    print()
    print("="*60)
    print("                    Pipeline Configuration")
    print("="*60)
    print(f"Config file: {args.config}")
    print(f"Query directory: {args.query_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Database: {database}")
    print(f"Threads: {threads} (Workers: {n_workers} × {threads_per_worker} threads)")
    print(f"Tree method: {tree_method}")
    print(f"Fast mode: {mode_fast}")
    print(f"Queries to process: {n_queries}")
    print("="*60)
    
    start_time = datetime.now()
    
    # Build command
    cmd = [
        sys.executable, "-m", "src.bin.gvclass_prefect",
        args.query_dir,
        "-o", output_dir,
        "-d", database,
        "-t", str(threads),
        "-j", str(n_workers),
        "--threads-per-worker", str(threads_per_worker),
        "--tree-method", tree_method,
        "--cluster-type", args.cluster_type
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
            cmd,
            env=env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    
    # Start resource monitor
    monitor = ResourceMonitor(output_dir, n_queries, process.pid)
    monitor_thread = threading.Thread(target=monitor.monitor)
    monitor_thread.start()
    
    try:
        # Wait for process
        return_code = process.wait()
        
        # Stop monitor
        monitor.stop()
        monitor_thread.join()
        print("\r" + " "*80 + "\r", end="")  # Clear line
        
        if return_code == 0:
            elapsed = datetime.now() - start_time
            minutes = int(elapsed.total_seconds() / 60)
            seconds = int(elapsed.total_seconds() % 60)
            
            print("\n✅ Pipeline completed successfully!")
            print(f"⏱️  Runtime: {minutes}m {seconds}s")
            print(f"💾 Peak memory usage: {monitor.peak_memory_mb:.0f} MB")
            print(f"📊 Results saved to: {output_dir}")
        else:
            print(f"\n❌ Pipeline failed with exit code: {return_code}")
            sys.exit(1)
            
    except KeyboardInterrupt:
        monitor.stop()
        monitor_thread.join()
        print("\r" + " "*80 + "\r", end="")
        print("\n⚠️  Pipeline interrupted by user")
        process.terminate()
        sys.exit(1)
    except Exception as e:
        monitor.stop()
        monitor_thread.join()
        print("\r" + " "*80 + "\r", end="")
        print(f"\n❌ Pipeline failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()