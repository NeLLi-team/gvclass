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

DOCKER_TAG="gvclass:v2.0-dev"
SIF_IMAGE="gvclass_v2.0-dev.sif"

# Setup temporary directory for Apptainer/Singularity builds
# Use /tmp to avoid recursive copy issues when temp dir is inside build context
BUILD_TMPDIR="/tmp/gvclass-build-$$"
mkdir -p "$BUILD_TMPDIR"
trap "rm -rf $BUILD_TMPDIR" EXIT  # Clean up on exit

CONTAINER_RUNTIME=""
RUNTIME_TMPDIR_VAR=""
if command -v apptainer &> /dev/null; then
    CONTAINER_RUNTIME="apptainer"
    RUNTIME_TMPDIR_VAR="APPTAINER_TMPDIR"
elif command -v singularity &> /dev/null; then
    CONTAINER_RUNTIME="singularity"
    RUNTIME_TMPDIR_VAR="SINGULARITY_TMPDIR"
fi

build_docker_image() {
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}Docker not found, so Docker-to-SIF fallback is unavailable.${NC}"
        return 1
    fi

    echo -e "\n${BLUE}Building Docker image ${DOCKER_TAG}...${NC}"
    docker build -t "${DOCKER_TAG}" -f containers/docker/Dockerfile .
    docker tag "${DOCKER_TAG}" gvclass:latest
    echo -e "${GREEN}✓ Docker image built successfully${NC}"
    docker images "${DOCKER_TAG}"
}

build_sif_from_docker() {
    if [ -z "${CONTAINER_RUNTIME}" ]; then
        echo -e "${RED}Neither apptainer nor singularity found; cannot create ${SIF_IMAGE}.${NC}"
        return 1
    fi

    echo -e "\n${BLUE}Converting ${DOCKER_TAG} to ${SIF_IMAGE}...${NC}"
    env "${RUNTIME_TMPDIR_VAR}=${BUILD_TMPDIR}" \
        "${CONTAINER_RUNTIME}" build --force "${SIF_IMAGE}" "docker-daemon://${DOCKER_TAG}"
    echo -e "${GREEN}✓ SIF image built successfully${NC}"
    ls -lah "${SIF_IMAGE}"
}

direct_build_ok=false
if [ -n "${CONTAINER_RUNTIME}" ]; then
    echo -e "\n${BLUE}Building ${SIF_IMAGE} directly from Apptainer definition (includes v2.0.0 database)...${NC}"
    if env "${RUNTIME_TMPDIR_VAR}=${BUILD_TMPDIR}" \
        "${CONTAINER_RUNTIME}" build --force "${SIF_IMAGE}" containers/apptainer/gvclass.def; then
        direct_build_ok=true
        echo -e "${GREEN}✓ SIF image built successfully${NC}"
        ls -lah "${SIF_IMAGE}"
    else
        echo -e "${RED}✗ Direct definition build failed; trying Docker-to-SIF fallback.${NC}"
    fi
else
    echo -e "\n${RED}Neither apptainer nor singularity found.${NC}"
    echo -e "Docker image can still be built, but SIF conversion requires Apptainer or Singularity."
fi

if [ "${direct_build_ok}" = false ]; then
    build_docker_image
    build_sif_from_docker
fi

echo -e "\n${GREEN}Container build complete!${NC}"
echo -e "\nUsage:"
echo -e "  Apptainer:  apptainer run -B /path/to/data:/data ${SIF_IMAGE} /data/input_dir -t 32"
echo -e "  Docker:     docker run -v /path/to/data:/data ${DOCKER_TAG} /data -t 32"
echo -e "\nNote: The container includes the v2.0.0 compact database."
