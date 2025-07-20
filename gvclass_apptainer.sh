#!/bin/bash
set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <query_directory> <num_processes>"
    echo "  query_directory: Path to directory containing query files (must be under current directory)"
    echo "  num_processes: Number of parallel processes (positive integer, max 32)"
    echo ""
    echo "Example:"
    echo "  $0 example 8"
    exit 1
}

# Validate argument count
if [ $# -ne 2 ]; then
    echo "Error: Expected 2 arguments, got $#"
    usage
fi

QUERYDIR="$1"
PROCESSES="$2"

# Validate query directory
if [ -z "$QUERYDIR" ]; then
    echo "Error: Query directory cannot be empty"
    usage
fi

if [ ! -d "$QUERYDIR" ]; then
    echo "Error: Query directory does not exist: $QUERYDIR"
    exit 1
fi

# Validate number of processes
if ! [[ "$PROCESSES" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of processes must be a positive integer"
    exit 1
fi

if [ "$PROCESSES" -gt 32 ]; then
    echo "Error: Number of processes cannot exceed 32"
    exit 1
fi

# Ensure query directory is under current working directory (security check)
QUERYDIR_ABS=$(realpath "$QUERYDIR")
CWD_ABS=$(realpath "$(pwd)")
if [[ ! "$QUERYDIR_ABS" == "$CWD_ABS"* ]]; then
    echo "Error: Query directory must be under current working directory"
    echo "  Current directory: $CWD_ABS"
    echo "  Query directory: $QUERYDIR_ABS"
    exit 1
fi

# Check if query directory contains valid input files
VALID_FILES=$(find "$QUERYDIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.faa" \) | wc -l)
if [ "$VALID_FILES" -eq 0 ]; then
    echo "Error: No valid input files (*.fna or *.faa) found in query directory: $QUERYDIR"
    exit 1
fi

echo "Starting GVClass with Apptainer..."
echo "Query directory: $QUERYDIR"
echo "Number of processes: $PROCESSES"
echo "Valid input files found: $VALID_FILES"

apptainer run --containall --bind "$(pwd):/workdir" --pwd /workdir \
  docker://docker.io/doejgi/gvclass:latest \
  pixi run snakemake --snakefile /gvclass/workflow/Snakefile \
           -j "$PROCESSES" \
           --config querydir="/workdir/$QUERYDIR" \
           database_path="/gvclass/resources"
