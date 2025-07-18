#!/usr/bin/env python3
"""
Simple clean output wrapper for Snakemake pipeline
"""
import subprocess
import sys
import os

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
    
    # Build command
    cmd = ["pixi", "run", "snakemake", "-j", threads]
    
    # Add config arguments
    if config_args:
        cmd.extend(["--config"] + config_args)
    
    print(f"  Configuration: {' '.join(config_args)}")
    print(f"  Threads: {threads}")
    print("=" * 50)
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1
    )
    
    # Track progress
    errors = []
    
    for line in iter(process.stdout.readline, ''):
        line = line.strip()
        
        # Show progress updates
        if "steps" in line and "done" in line:
            # Clear previous line and show new progress
            print(f"\r  Progress: {line}", end="        ", flush=True)
        
        # Capture errors
        elif "Error" in line or "failed" in line:
            errors.append(line)
            if "Error in rule" in line or "CalledProcessError" in line:
                print(f"\n  ❌ {line}")
    
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