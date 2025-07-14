#!/bin/bash
#SBATCH -J gvclass_array
#SBATCH -o logs/gvclass_array_%A_%a.out
#SBATCH -e logs/gvclass_array_%A_%a.err
#SBATCH -c 24
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --account=grp-org-sc-mgs
#SBATCH --qos=jgi_normal

# Default parameters
QUERYDIR=allbins""
OUTDIR="allbins_results"
DATABASE_PATH=""
QUERIES_PER_CHUNK=50
MODE_FAST=true
MAFFTOPTION="linsi"
TREEOPTION="iqtree"
RESUME=false
EXTRA_ARGS=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --querydir)
            QUERYDIR="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --database_path)
            DATABASE_PATH="$2"
            shift 2
            ;;
        --queries-per-chunk)
            QUERIES_PER_CHUNK="$2"
            shift 2
            ;;
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --mode_fast)
            MODE_FAST="$2"
            shift 2
            ;;
        --mafftoption)
            MAFFTOPTION="$2"
            shift 2
            ;;
        --treeoption)
            TREEOPTION="$2"
            shift 2
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        --time)
            SBATCH_TIME="$2"
            shift 2
            ;;
        --mem)
            SBATCH_MEM="$2"
            shift 2
            ;;
        --cpus)
            SBATCH_CPUS="$2"
            shift 2
            ;;
        *)
            EXTRA_ARGS="$EXTRA_ARGS $1"
            shift
            ;;
    esac
done

# Check required parameters
if [ -z "$QUERYDIR" ]; then
    echo "Error: --querydir is required"
    exit 1
fi

# Set default output directory if not specified
if [ -z "$OUTDIR" ]; then
    OUTDIR="${QUERYDIR}_results"
fi

# Get absolute paths
QUERYDIR=$(realpath "$QUERYDIR")
OUTDIR=$(realpath "$OUTDIR")

# Set PROJECT_DIR if not provided
if [ -z "$PROJECT_DIR" ]; then
    # Get the directory containing this script
    SCRIPT_DIR=$(dirname $(realpath "$0"))
    # If we're in a bin directory, go up one level to get project root
    if [[ "$SCRIPT_DIR" == */bin ]]; then
        PROJECT_DIR=$(dirname "$SCRIPT_DIR")
    else
        PROJECT_DIR="$SCRIPT_DIR"
    fi
fi

# Set database path if not specified
if [ -z "$DATABASE_PATH" ]; then
    DATABASE_PATH="${PROJECT_DIR}/resources"
fi

# Create necessary directories
mkdir -p logs
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/chunks"

# Get list of query files
QUERY_FILES=($(find "$QUERYDIR" -name "*.fna" -o -name "*.faa" | sort))
TOTAL_QUERIES=${#QUERY_FILES[@]}

if [ $TOTAL_QUERIES -eq 0 ]; then
    echo "Error: No .fna or .faa files found in $QUERYDIR"
    exit 1
fi

echo "Found $TOTAL_QUERIES query files"

# Calculate number of chunks needed
CHUNKS=$(( (TOTAL_QUERIES + QUERIES_PER_CHUNK - 1) / QUERIES_PER_CHUNK ))
echo "Creating $CHUNKS chunks with up to $QUERIES_PER_CHUNK queries each"

# Calculate chunk size
CHUNK_SIZE=$QUERIES_PER_CHUNK

# Create chunk directories with symlinks
echo "Creating chunk directories..."
for ((i=0; i<$CHUNKS; i++)); do
    CHUNK_DIR="$OUTDIR/chunks/chunk_${i}"
    mkdir -p "$CHUNK_DIR"
    
    # Calculate start and end indices for this chunk
    START=$((i * CHUNK_SIZE))
    END=$((START + CHUNK_SIZE))
    if [ $END -gt $TOTAL_QUERIES ]; then
        END=$TOTAL_QUERIES
    fi
    
    # Create symlinks for files in this chunk
    for ((j=$START; j<$END; j++)); do
        if [ $j -lt $TOTAL_QUERIES ]; then
            ln -sf "${QUERY_FILES[$j]}" "$CHUNK_DIR/"
        fi
    done
    
    # Count files in chunk
    CHUNK_FILES=$(find "$CHUNK_DIR" -name "*.fna" -o -name "*.faa" | wc -l)
    echo "Chunk $i: $CHUNK_FILES files"
    
    # Check if chunk has only .faa files (no .fna files)
    FNA_COUNT=$(find "$CHUNK_DIR" -name "*.fna" | wc -l)
    if [ $FNA_COUNT -eq 0 ] && [ $CHUNK_FILES -gt 0 ]; then
        # Create a dummy .fna file to satisfy Nextflow's requirement
        touch "$CHUNK_DIR/dummy.fna"
        echo "  Added dummy.fna to chunk $i (protein-only input)"
    fi
done

# Submit array job if not already in SLURM job
if [ -z "$SLURM_JOB_ID" ]; then
    # Update SBATCH directives if custom values provided
    SBATCH_OPTS=""
    if [ ! -z "$SBATCH_TIME" ]; then
        SBATCH_OPTS="$SBATCH_OPTS --time=$SBATCH_TIME"
    fi
    if [ ! -z "$SBATCH_MEM" ]; then
        SBATCH_OPTS="$SBATCH_OPTS --mem=$SBATCH_MEM"
    fi
    if [ ! -z "$SBATCH_CPUS" ]; then
        SBATCH_OPTS="$SBATCH_OPTS --cpus-per-task=$SBATCH_CPUS"
    fi
    
    echo "Submitting SLURM array job with $CHUNKS tasks..."
    sbatch --array=0-$((CHUNKS-1)) $SBATCH_OPTS "$0" \
        --querydir "$QUERYDIR" \
        --outdir "$OUTDIR" \
        --database_path "$DATABASE_PATH" \
        --project-dir "$PROJECT_DIR" \
        --queries-per-chunk "$QUERIES_PER_CHUNK" \
        --mode_fast "$MODE_FAST" \
        --mafftoption "$MAFFTOPTION" \
        --treeoption "$TREEOPTION" \
        $( [ "$RESUME" = true ] && echo "--resume" ) \
        $EXTRA_ARGS
    exit 0
fi

# This part runs within the SLURM job
echo "Running chunk $SLURM_ARRAY_TASK_ID on $(hostname)"
echo "Start time: $(date)"

# Ensure PROJECT_DIR is set correctly
# The PROJECT_DIR should have been passed via --export from the submission
echo "Project directory: $PROJECT_DIR"

# Verify the project directory exists and contains expected files
if [ ! -f "${PROJECT_DIR}/main.nf" ]; then
    echo "ERROR: main.nf not found in PROJECT_DIR: $PROJECT_DIR"
    echo "This suggests PROJECT_DIR was not correctly passed from submission."
    exit 1
fi

# Set up environment
# Check if pixi is available (we'll use pixi run for nextflow)
if [ -f "${PROJECT_DIR}/pixi.toml" ] && command -v pixi &> /dev/null; then
    echo "Using pixi environment for Nextflow execution"
else
    # Fall back to modules if available
    if command -v module &> /dev/null 2>&1; then
        echo "Loading modules..."
        # Load Java 11 or higher for Nextflow
        module purge 2>/dev/null || true
        module load java/11 2>/dev/null || module load openjdk/11 2>/dev/null || module load java/17 2>/dev/null || true
        module load nextflow 2>/dev/null || true
    fi
    
    # Verify Nextflow is available
    if ! command -v nextflow &> /dev/null; then
        echo "ERROR: Nextflow not found. Please ensure pixi is installed or nextflow is in PATH."
        exit 1
    fi
fi

# Set up environment
export NXF_EXECUTOR=local
export NXF_WORK="$OUTDIR/work_chunk_${SLURM_ARRAY_TASK_ID}"

# Get chunk directory
CHUNK_DIR="$OUTDIR/chunks/chunk_${SLURM_ARRAY_TASK_ID}"
CHUNK_OUTDIR="$OUTDIR/chunk_${SLURM_ARRAY_TASK_ID}_results"

# Build nextflow command
# Use pixi if available, otherwise use system nextflow
if [ -f "${PROJECT_DIR}/pixi.toml" ] && command -v pixi &> /dev/null; then
    # Change to project directory for pixi
    cd "$PROJECT_DIR"
    NF_CMD="pixi run nextflow run ${PROJECT_DIR}/main.nf"
else
    NF_CMD="nextflow run ${PROJECT_DIR}/main.nf"
fi
NF_CMD="$NF_CMD --querydir $CHUNK_DIR"
NF_CMD="$NF_CMD --outdir $CHUNK_OUTDIR"
NF_CMD="$NF_CMD --database_path $DATABASE_PATH"
NF_CMD="$NF_CMD --mode_fast $MODE_FAST"
NF_CMD="$NF_CMD --mafftoption $MAFFTOPTION"
NF_CMD="$NF_CMD --treeoption $TREEOPTION"
NF_CMD="$NF_CMD -process.cpus=$SLURM_CPUS_PER_TASK"
NF_CMD="$NF_CMD -process.memory=${SLURM_MEM_PER_NODE}MB"

if [ "$RESUME" = true ]; then
    NF_CMD="$NF_CMD -resume"
fi

# Add any extra arguments
NF_CMD="$NF_CMD $EXTRA_ARGS"

echo "Executing: $NF_CMD"

# Run nextflow
$NF_CMD

# Check exit status
if [ $? -eq 0 ]; then
    echo "Chunk $SLURM_ARRAY_TASK_ID completed successfully"
    
    # Move results to main output directory
    echo "Moving results to main output directory..."
    for result_file in "$CHUNK_OUTDIR"/*; do
        if [ -f "$result_file" ]; then
            base_name=$(basename "$result_file")
            # Append chunk ID to avoid conflicts
            new_name="${base_name%.txt}_chunk${SLURM_ARRAY_TASK_ID}.txt"
            cp "$result_file" "$OUTDIR/$new_name"
        fi
    done
else
    echo "Chunk $SLURM_ARRAY_TASK_ID failed with exit code $?"
    exit 1
fi

echo "End time: $(date)"
