# GVClass Container Images

This directory contains container definitions for GVClass v1.4.0.

## Directory Structure

```
containers/
├── build.sh           # Main build script
├── apptainer/        
│   └── gvclass.def    # Apptainer/Singularity definition (primary)
└── docker/           
    ├── Dockerfile     # Docker image definition (alternative)
    └── docker-compose.yml  # Docker Compose configuration
```

## Apptainer/Singularity (HPC)

Common container format for HPC environments.

### Quick Build

```bash
# From project root
bash containers/build.sh

# Or build directly
apptainer build gvclass.sif containers/apptainer/gvclass.def
```

This creates a ~992MB self-contained image including:
- Complete Pixi environment
- All Python dependencies  
- Full reference database downloaded during the image build
- GVClass source code

### Usage

```bash
# Basic usage
singularity run -B /path/to/data:/data gvclass.sif /data/input_dir -t 32

# With custom output location
singularity run -B /data:/data -B /results:/results \
    gvclass.sif /data/input -o /results -t 32

# Override with external database (optional)
singularity run \
    -B /path/to/queries:/input \
    -B /path/to/database:/opt/gvclass/resources \
    gvclass.sif /input -t 16
```

### Publish to Apptainer Library (for `apptainer pull`)

To support `apptainer pull library://nelligroup-jgi/gvclass/gvclass:1.4.0`:

```bash
# Build the SIF
apptainer build gvclass.sif containers/apptainer/gvclass.def

# Authenticate to the Sylabs library (one-time)
apptainer remote login

# Push the image
apptainer push gvclass.sif library://nelligroup-jgi/gvclass/gvclass:1.4.0
```

## Alternative: Docker

Docker images are available but less suitable for HPC clusters.

### Build
```bash
# From project root
docker build -t gvclass:1.4.0 -f containers/docker/Dockerfile .
```

### Usage
```bash
docker run -v /path/to/data:/data gvclass:1.4.0 /data -t 32
```

## Container Features

Both formats include:
- **Database**: full reference database embedded in the final image; fetched during image creation rather than tracked in git
- **Dependencies**: All tools installed via Pixi
- **Python packages**: pyhmmer, pyrodigal, veryfasttree, etc.
- **Workflow**: direct GVClass pipeline execution with multi-threaded query parallelism

## Container Sizes

- **Apptainer**: ~992MB (compressed SIF, recommended)
- **Docker**: ~4.4GB (uncompressed layers)

## Notes

- No sudo/root required for Apptainer
- The final image embeds the runtime database, but the repo checkout keeps `resources/` local-only
- Supports bind mounting for input/output
- Compatible with HPC schedulers (SLURM, PBS, SGE)
- The `.dockerignore` file must remain in project root for Docker builds
