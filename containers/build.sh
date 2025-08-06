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

# Primary method: Build Apptainer/Singularity image directly
if command -v apptainer &> /dev/null; then
    echo -e "\n${BLUE}Building Apptainer image (includes 850MB database)...${NC}"
    if apptainer build --force gvclass.sif containers/apptainer/gvclass.def; then
        echo -e "${GREEN}✓ Apptainer image built successfully${NC}"
        ls -lah gvclass.sif
    else
        echo -e "${RED}✗ Apptainer build failed${NC}"
        exit 1
    fi
elif command -v singularity &> /dev/null; then
    echo -e "\n${BLUE}Building Singularity image (includes 850MB database)...${NC}"
    if singularity build --force gvclass.sif containers/apptainer/gvclass.def; then
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
        if docker build -t gvclass:1.1.0 -f containers/docker/Dockerfile .; then
            docker tag gvclass:1.1.0 gvclass:latest
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
echo -e "  Docker:     docker run -v /path/to/data:/data gvclass:1.1.0 /data -t 32"
echo -e "\nNote: The container includes the complete 850MB database."