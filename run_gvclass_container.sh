#!/bin/bash
# Wrapper script to run GVClass container with proper bind mounts

# Get absolute path of input directory
INPUT_DIR=$(realpath "$1")
shift  # Remove first argument to pass the rest to the pipeline

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Get the directory containing this script (where the SIF should be)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SIF_FILE="$SCRIPT_DIR/gvclass.sif"

# Check if SIF exists
if [ ! -f "$SIF_FILE" ]; then
    echo "Error: Container not found at $SIF_FILE"
    echo "Please build it first with: apptainer build --fakeroot gvclass.sif gvclass-pixi-flexible.def"
    exit 1
fi

# Run the container with proper bind mounts
# Bind the parent directory of input to allow output alongside input
PARENT_DIR=$(dirname "$INPUT_DIR")
INPUT_NAME=$(basename "$INPUT_DIR")

echo "Running GVClass on: $INPUT_DIR"
echo "Container: $SIF_FILE"

# Process arguments to handle output directory
PROCESSED_ARGS=""
OUTPUT_SPECIFIED=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output-dir)
            # Convert relative output path to absolute under /data
            OUTPUT_DIR="$2"
            if [[ "$OUTPUT_DIR" != /* ]]; then
                # Relative path - put it under /data
                OUTPUT_DIR="/data/$OUTPUT_DIR"
            fi
            PROCESSED_ARGS="$PROCESSED_ARGS -o $OUTPUT_DIR"
            OUTPUT_SPECIFIED=true
            shift 2
            ;;
        *)
            PROCESSED_ARGS="$PROCESSED_ARGS $1"
            shift
            ;;
    esac
done

# Add default output if not specified
if [ "$OUTPUT_SPECIFIED" = false ]; then
    PROCESSED_ARGS="-o /data/${INPUT_NAME}_results $PROCESSED_ARGS"
fi

apptainer run \
    -B "$PARENT_DIR:/data" \
    --pwd /opt/gvclass \
    "$SIF_FILE" \
    "/data/$INPUT_NAME" \
    $PROCESSED_ARGS