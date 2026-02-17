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

# Change to project root for build context
cd "$PROJECT_ROOT"

# Ensure database exists before building container
if [ ! -d "resources/database" ] || [ ! -f "resources/models/combined.hmm" ]; then
    echo -e "${BLUE}Database not found locally. Downloading...${NC}"
    if command -v pixi &> /dev/null; then
        pixi run python -c 'from src.utils.database_manager import DatabaseManager; DatabaseManager.setup_database("resources")'
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Database downloaded successfully${NC}"
        else
            echo -e "${RED}✗ Failed to download database${NC}"
            exit 1
        fi
    else
        echo -e "${RED}Pixi not found. Please install pixi and run 'pixi install' first.${NC}"
        exit 1
    fi
fi

# Setup temporary directory for Apptainer/Singularity builds
# Use /tmp to avoid recursive copy issues when temp dir is inside build context
BUILD_TMPDIR="/tmp/gvclass-build-$$"
mkdir -p "$BUILD_TMPDIR"
trap "rm -rf $BUILD_TMPDIR" EXIT  # Clean up on exit

# Primary method: Build Apptainer/Singularity image directly
if command -v apptainer &> /dev/null; then
    echo -e "\n${BLUE}Building Apptainer image (includes 850MB database)...${NC}"
    if APPTAINER_TMPDIR="$BUILD_TMPDIR" apptainer build --force gvclass.sif containers/apptainer/gvclass.def; then
        echo -e "${GREEN}✓ Apptainer image built successfully${NC}"
        ls -lah gvclass.sif
    else
        echo -e "${RED}✗ Apptainer build failed${NC}"
        exit 1
    fi
elif command -v singularity &> /dev/null; then
    echo -e "\n${BLUE}Building Singularity image (includes 850MB database)...${NC}"
    if SINGULARITY_TMPDIR="$BUILD_TMPDIR" singularity build --force gvclass.sif containers/apptainer/gvclass.def; then
        echo -e "${GREEN}✓ Singularity image built successfully${NC}"
        ls -lah gvclass.sif
    else
        echo -e "${RED}✗ Singularity build failed${NC}"
        exit 1
    fi
else
    echo -e "\n${RED}Neither apptainer nor singularity found.${NC}"
    echo -e "Please install Apptainer or Singularity first:"
    echo -e "  https://apptainer.org/docs/admin/main/installation.html"
    exit 1
fi

# Optional: Build Docker image if docker is available
if command -v docker &> /dev/null; then
    echo -e "\n${BLUE}Optional: Build Docker image? (y/N)${NC}"
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        echo -e "${BLUE}Building Docker image...${NC}"
        if docker build -t gvclass:1.2.2 -f containers/docker/Dockerfile .; then
            docker tag gvclass:1.2.2 gvclass:latest
            echo -e "${GREEN}✓ Docker image built successfully${NC}"
            docker images | grep gvclass
        else
            echo -e "${RED}✗ Docker build failed${NC}"
        fi
    fi
fi

echo -e "\n${GREEN}Container build complete!${NC}"
echo -e "\nUsage:"
echo -e "  Apptainer:  singularity run -B /path/to/data:/data gvclass.sif /data/input_dir -t 32"
echo -e "  Docker:     docker run -v /path/to/data:/data gvclass:1.2.2 /data -t 32"
echo -e "\nNote: The container includes the complete 850MB database."