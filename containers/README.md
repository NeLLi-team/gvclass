# GVClass Container Images

This directory contains container definitions for the temporary GVClass `v2.0-dev` test image.

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

# If direct def-file builds are unavailable, build Docker then convert
docker build -t gvclass:v2.0-dev -f containers/docker/Dockerfile .
apptainer build --force gvclass_v2.0-dev.sif docker-daemon://gvclass:v2.0-dev
```

This creates a self-contained image including:
- Complete Pixi environment
- All Python dependencies  
- Compact v2.0.0 reference database downloaded during the image build
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
    -B /path/to/cache:/tmp/gvclass-resource-cache \
    --env GVCLASS_RESOURCE_CACHE=/tmp/gvclass-resource-cache \
    gvclass.sif /input -t 16
```

The v2.0.0 resource bundle stores labels and marker proteins as Parquet inside
the read-only SIF. GVClass materializes the TSV/FASTA views it needs under
`GVCLASS_RESOURCE_CACHE`. Bind a host cache directory to
`/tmp/gvclass-resource-cache` for a persistent warm cache; otherwise the cache
falls back to container `/tmp` and is rebuilt on later runs.

### Temporary Published SIF

The current `gvclass-dev` test image is published here:

```text
https://dl.newlineages.com/gvclass/gvclass_v2.0-dev.sif
```

It is a temporary v2.0-dev SIF with the v2.0.0 database embedded. The
`gvclass-a` wrapper downloads this URL into `~/.cache/gvclass/images/` and
reuses it on later runs.

### Optional Sylabs Mirror

To support `apptainer pull library://nelligroup-jgi/gvclass/gvclass:v2.0-dev`
after a Sylabs token is available:

```bash
# Build the SIF
apptainer build gvclass.sif containers/apptainer/gvclass.def

# Authenticate to the Sylabs library (one-time)
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
apptainer remote login SylabsCloud

# Push the image
apptainer push -U --library https://library.sylabs.io \
    gvclass.sif library://nelligroup-jgi/gvclass/gvclass:v2.0-dev
```

## Alternative: Docker

Docker images are available but less suitable for HPC clusters.

### Build
```bash
# From project root
docker build -t gvclass:v2.0-dev -f containers/docker/Dockerfile .
```

### Usage
```bash
docker run -v /path/to/data:/data gvclass:v2.0-dev /data -t 32
```

## Container Features

Both formats include:
- **Database**: compact v2.0.0 reference database embedded in the final image; fetched during image creation rather than tracked in git
- **Dependencies**: All tools installed via Pixi
- **Python packages**: pyhmmer, pyrodigal, veryfasttree, etc.
- **Workflow**: direct GVClass pipeline execution with multi-threaded query parallelism

## Container Sizes

- **Apptainer**: ~2.1GB compressed SIF, recommended for HPC
- **Docker**: ~2.9GB content size

## Notes

- No sudo/root required to run the published Apptainer image on a normal HPC Apptainer installation
- Building from `containers/apptainer/gvclass.def` may require setuid/fakeroot support; Docker-to-SIF conversion is the fallback
- The final image embeds the runtime database, but the repo checkout keeps `resources/` local-only
- Supports bind mounting for input/output
- Compatible with HPC schedulers (SLURM, PBS, SGE)
- The `.dockerignore` file must remain in project root for Docker builds
