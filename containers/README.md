# GVClass Container Images

This directory contains container definitions prepared for the GVClass `2.2.0` release image.

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
docker build -t gvclass:2.2.0 -f containers/docker/Dockerfile .
apptainer build --force gvclass_2.2.0.sif docker-daemon://gvclass:2.2.0
```

The recipes build a self-contained image that includes:

- The Pixi environment with all Python dependencies (pyhmmer, pyrodigal, veryfasttree and others)
- The compact v2.1.0 reference database, downloaded during the image build
- GVClass source code, including the contig-screen code (the screen add-on is not in the image)

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
    -B /path/to/cache:/resource-cache \
    --env GVCLASS_RESOURCE_CACHE=/resource-cache \
    gvclass.sif /input -t 16

# The screen add-on is not in the image. Install it and run classifications and
# screens through the launcher (./gvclass ... --image gvclass.sif), which mounts it:
#   ./gvclass screen --install --image gvclass.sif
# With a plain `singularity run`, mount it where the image looks for it:
#   -B /path/to/cache/contig_screen:/opt/gvclass/resources/contig_screen:ro
# The Docker container runs as root and writes the add-on into its own database
# folder, which is lost with the container.
```

The v2.1.0 resource bundle stores labels and marker proteins as Parquet inside
the read-only SIF. GVClass materializes the TSV/FASTA views it needs under
`GVCLASS_RESOURCE_CACHE`. Bind a host cache directory to
`/resource-cache` for a persistent warm cache; otherwise the cache
falls back to container `/tmp` and is rebuilt on later runs.

### Release SIF

After the v2.1.0 resource and model checks pass, the workflow publishes the
approved image to GHCR when a release tag is pushed, as the version and as `latest`:

```text
oras://ghcr.io/nelli-team/gvclass:2.2.0
```

Running `./gvclass --apptainer` (or a standalone `./gvclass` downloaded without a repo checkout) pulls the image from this address into
`~/.cache/gvclass/images/` and reuses it on later runs. The address is unavailable until publication completes.
Sylabs retired its library on 2026-10-22, so older `library://` addresses no longer work.

The workflow then pulls the version back without credentials and compares its
SHA-256 with the built file. A new GHCR package is private until its visibility is
set to public once in the package settings; until then that check fails, because
users cannot pull a private package.

A manually dispatched workflow builds and tests the image without publishing,
including when the selected ref is a release tag.

### GHCR Publish

To publish a release SIF manually:

```bash
# Build the SIF
apptainer build gvclass.sif containers/apptainer/gvclass.def

# Authenticate to GHCR (one-time): a token with the write:packages scope
printf '%s' "$GHCR_TOKEN" | apptainer registry login \
    --username "$GITHUB_USER" --password-stdin oras://ghcr.io

# Push the image
apptainer push gvclass.sif oras://ghcr.io/nelli-team/gvclass:2.2.0
```

## Alternative: Docker

Docker images are available but less suitable for HPC clusters.

### Build
```bash
# From project root
docker build -t gvclass:2.2.0 -f containers/docker/Dockerfile .
```

### Usage
```bash
docker run -v /path/to/data:/data gvclass:2.2.0 /data -t 32
```

## Container Sizes

- **Apptainer**: ~2.1GB compressed SIF, recommended for HPC
- **Docker**: ~2.9GB content size

## Notes

- Both recipes copy only the runtime source and the bundled example. They exclude tests, host wrappers, images, and build files
- Both recipes use the verified Pixi 0.72.0 base-image digest
- Docker Compose builds from the repository root
- No sudo/root required to run the published Apptainer image on a normal HPC Apptainer installation
- Building from `containers/apptainer/gvclass.def` may require setuid/fakeroot support; Docker-to-SIF conversion is the fallback
- The final image embeds the runtime database, but the repo checkout keeps `resources/` local-only
- Compatible with HPC schedulers (SLURM, PBS, SGE)
- The `.dockerignore` file must remain in project root for Docker builds
