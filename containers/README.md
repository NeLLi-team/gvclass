# GVClass Container Build Files

This directory contains all files needed to build GVClass containers.

## Directory Structure

```
containers/
├── build.sh           # Main build script for all containers
├── docker/           
│   ├── Dockerfile     # Docker image definition
│   ├── docker-compose.yml  # Docker Compose configuration
│   └── build_containers.sh # Legacy build script
└── apptainer/
    └── gvclass.def    # Apptainer/Singularity definition (optional)
```

## Building Containers

### Quick Build (Recommended)
```bash
# From project root
bash containers/build.sh
```

This will:
1. Build Docker image `gvclass:1.1.0` (4.38GB)
2. Build Apptainer image `gvclass.sif` (1.1GB) if apptainer/singularity is available

### Manual Build

#### Docker
```bash
# From project root
docker build -t gvclass:1.1.0 -f containers/docker/Dockerfile .
docker tag gvclass:1.1.0 gvclass:latest
```

#### Apptainer/Singularity
```bash
# From project root (requires Docker image first)
apptainer build gvclass.sif docker-daemon://gvclass:1.1.0
# or
singularity build gvclass.sif docker-daemon://gvclass:1.1.0
```

## Container Sizes

- **Docker**: 4.38GB (uncompressed layers)
- **Apptainer**: 1.1GB (compressed SIF format)

Both contain:
- Complete Python environment (~2GB)
- GVClass database (~832MB)
- All dependencies and tools

## Notes

- The `.dockerignore` file must remain in the project root for Docker builds
- Apptainer images are built from the Docker image for consistency
- Both containers include the full database, no download needed at runtime