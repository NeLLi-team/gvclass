#!/usr/bin/env python3
"""
Simple clean output wrapper for Snakemake pipeline
"""
import subprocess
import sys
import os
from pathlib import Path

def main():
    # Parse command line arguments
    args = sys.argv[1:]
    
    # Default values
    threads = "8"
    config_args = []
    has_querydir = False
    querydir = None
    
    # Process arguments
    i = 0
    while i < len(args):
        if args[i] == "-j" and i + 1 < len(args):
            threads = args[i + 1]
            i += 2
        elif args[i] == "--config" and i + 1 < len(args):
            config_args = args[i + 1].split()
            # Check if querydir is specified
            for config_item in config_args:
                if config_item.startswith("querydir="):
                    has_querydir = True
                    querydir = config_item.split("=")[1].strip('"')
            i += 2
        elif args[i] == "--configfile" and i + 1 < len(args):
            # Pass through configfile arguments
            config_args.append(f"--configfile {args[i + 1]}")
            i += 2
        else:
            # Pass through other arguments
            config_args.append(args[i])
            i += 1
    
    # Default querydir if not specified
    if not has_querydir:
        config_args.append("querydir=example")
        querydir = "example"
    
    # If outdir not specified, create default based on querydir
    has_outdir = any("outdir=" in arg for arg in config_args)
    if not has_outdir and querydir:
        outdir = f"{os.path.basename(querydir)}_results"
        config_args.append(f"outdir={outdir}")
    
    print("🚀 Starting GVClass Pipeline v1.1.0")
    print("=" * 50)
    
    # Check for database
    workflow_dir = Path(__file__).parent / "workflow"
    database_path = None
    
    # Check if custom database_path is specified
    for arg in config_args:
        if arg.startswith("database_path="):
            database_path = Path(arg.split("=")[1].strip('"'))
            break
    
    # Default database path if not specified
    if database_path is None:
        database_path = workflow_dir / "resources"
    
    # Check database status
    if database_path.exists() and (database_path / "models" / "combined.hmm").exists():
        print("  ✅ Database found")
    else:
        print("  📦 Database not found. Will be downloaded (~3GB)...")
    
    # Build command
    cmd = ["pixi", "run", "snakemake", "-j", threads]
    
    # Add config arguments
    if config_args:
        cmd.extend(["--config"] + config_args)
    
    # Add database_path to config_args display if not already present
    has_database_path = any("database_path=" in arg for arg in config_args)
    if not has_database_path:
        config_display = config_args + [f"database_path=resources"]
    else:
        config_display = config_args
    
    print(f"  Configuration: {' '.join(config_display)}")
    print(f"  Threads: {threads}")
    print("=" * 50)
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1
    )
    
    # Track progress and state
    errors = []
    database_downloaded = False
    database_verified = False
    last_progress = ""
    
    for line in iter(process.stdout.readline, ''):
        line = line.strip()
        
        # Handle database messages with state tracking
        if "Resources directory not found" in line and not database_downloaded:
            # Clear progress line if exists
            if last_progress:
                print("\r" + " " * len(last_progress), end="\r")
            print("  📦 Downloading database (~3GB)...")
            database_downloaded = True
        elif "Resources downloaded and extracted" in line and database_downloaded:
            print("  ✅ Database downloaded successfully!")
        elif "Resources directory found" in line and not database_verified:
            # Only show if we didn't already show it in pre-flight check
            pass  # Skip this as we already showed it above
        elif "Resources verified" in line and not database_verified:
            if not database_downloaded:  # Only show if we didn't download
                pass  # Skip as we showed it in pre-flight
            database_verified = True
        
        # Show progress updates
        elif "steps" in line and "done" in line:
            # Store progress for potential clearing
            last_progress = f"  Progress: {line}        "
            print(f"\r{last_progress}", end="", flush=True)
        
        # Capture errors
        elif "Error" in line or "failed" in line:
            errors.append(line)
            if "Error in rule" in line or "CalledProcessError" in line:
                # Clear progress line if exists
                if last_progress:
                    print("\r" + " " * len(last_progress), end="\r")
                print(f"  ❌ {line}")
    
    # Wait for completion
    process.wait()
    
    # Final status
    print("\n" + "=" * 50)
    if process.returncode == 0:
        print("✅ Pipeline completed successfully!")
    else:
        print("❌ Pipeline failed with errors")
        if errors:
            print("\nLast errors:")
            for error in errors[-3:]:
                print(f"  • {error}")
    
    return process.returncode

if __name__ == "__main__":
    sys.exit(main())