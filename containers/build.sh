#!/bin/bash
# Build container images for GVClass

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}Building GVClass container images...${NC}"

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Change to project root for Docker context
cd "$PROJECT_ROOT"

# Build Docker image
echo -e "\n${BLUE}Building Docker image...${NC}"
if docker build -t gvclass:1.1.0 -f containers/docker/Dockerfile .; then
    docker tag gvclass:1.1.0 gvclass:latest
    echo -e "${GREEN}✓ Docker image built successfully${NC}"
    docker images | grep gvclass
else
    echo -e "${RED}✗ Docker build failed${NC}"
    exit 1
fi

# Build Apptainer/Singularity image if apptainer is available
if command -v apptainer &> /dev/null; then
    echo -e "\n${BLUE}Building Apptainer/Singularity image...${NC}"
    if apptainer build gvclass.sif docker-daemon://gvclass:1.1.0; then
        echo -e "${GREEN}✓ Apptainer image built successfully${NC}"
        ls -lah gvclass.sif
    else
        echo -e "${RED}✗ Apptainer build failed${NC}"
        exit 1
    fi
elif command -v singularity &> /dev/null; then
    echo -e "\n${BLUE}Building Singularity image...${NC}"
    if singularity build gvclass.sif docker-daemon://gvclass:1.1.0; then
        echo -e "${GREEN}✓ Singularity image built successfully${NC}"
        ls -lah gvclass.sif
    else
        echo -e "${RED}✗ Singularity build failed${NC}"
        exit 1
    fi
else
    echo -e "\n${RED}Neither apptainer nor singularity found. Skipping SIF build.${NC}"
    echo -e "To build Apptainer image later, run:"
    echo -e "  apptainer build gvclass.sif docker-daemon://gvclass:1.1.0"
fi

echo -e "\n${GREEN}Container build complete!${NC}"
echo -e "\nUsage:"
echo -e "  Docker:     docker run -v /path/to/data:/data -v /path/to/results:/results gvclass:1.1.0 pixi run gvclass /data -o /results"
echo -e "  Apptainer:  apptainer run -B /path/to/data:/data,/path/to/results:/results gvclass.sif pixi run gvclass /data -o /results"