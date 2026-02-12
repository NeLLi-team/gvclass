#!/bin/bash
set -euo pipefail

# Legacy build script - use containers/build.sh instead
echo "Note: This is a legacy script. Use containers/build.sh instead."
echo ""

# Get to project root
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$PROJECT_ROOT"

echo "Building GVClass containers from $PROJECT_ROOT..."

# Build Docker image
echo "Building Docker image..."
docker build -t gvclass:1.2.1 -t gvclass:latest -f containers/docker/Dockerfile .

# Export Docker image for Apptainer conversion
echo "Exporting Docker image..."
docker save gvclass:1.2.1 -o gvclass_docker.tar

# Build Apptainer/Singularity image from Docker image
if command -v apptainer &> /dev/null; then
    echo "Building Apptainer image..."
    apptainer build gvclass.sif docker-archive://gvclass_docker.tar
elif command -v singularity &> /dev/null; then
    echo "Building Singularity image..."
    singularity build gvclass.sif docker-archive://gvclass_docker.tar
else
    echo "Warning: Neither Apptainer nor Singularity found. Skipping SIF build."
    echo "Docker image built successfully: gvclass:1.2.1"
fi

# Clean up temporary file
rm -f gvclass_docker.tar

echo "Container build complete!"
echo ""
echo "Usage examples:"
echo "Docker:"
echo "  docker run -v /path/to/data:/data -v /path/to/results:/results gvclass:1.2.1 pixi run gvclass /data -o /results"
echo ""
if [ -f gvclass.sif ]; then
    echo "Apptainer/Singularity:"
    echo "  apptainer run -B /path/to/data:/data,/path/to/results:/results gvclass.sif pixi run gvclass /data -o /results"
fi