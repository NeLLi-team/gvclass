# Run on an HPC cluster

GVClass ships as a single Apptainer (Singularity) image that bundles the database and every dependency, so a cluster run needs no local install and no `setup-db` step. The `gvclass-a` wrapper handles the image pull, the bind mounts, and the container call for you.

## Get the wrapper

Download the wrapper from the `gvclass-dev` branch and make it executable.

```bash
wget https://raw.githubusercontent.com/NeLLi-team/gvclass/gvclass-dev/gvclass-a
chmod +x gvclass-a
```

## Run it

Pass the query directory first and the results directory second, then set the thread count.

```bash
./gvclass-a QUERY_DIR RESULTS_DIR -t 32
```

The positional order is query then results. This differs from the pixi CLI, where the output is the `-o` flag (see the [CLI reference](../reference/cli.md)). `QUERY_DIR` holds your `.fna` or `.faa` bins; `RESULTS_DIR` is created if it does not exist.

On first use the wrapper pulls `library://nelligroup-jgi/gvclass/gvclass:2.0.0` from the public Sylabs library (no access token) and caches the SIF under `~/.cache/gvclass/images/`. Later runs reuse the cached image. The image carries the database (about 700 MB) and all tools, so you skip database setup entirely.

!!! note
    The wrapper calls `apptainer`, which must be on your `PATH`. Many clusters expose it through a module, for example `module load apptainer` or `module load singularity`.

## Wrapper options

The wrapper exposes the flags most runs need.

| Flag | Effect |
|------|--------|
| `-t, --threads N` | Total threads (wrapper default 16; the pixi CLI uses the config default of 4). Match your scheduler allocation. |
| `-j, --max-workers N` | Genomes classified in parallel (default auto). |
| `--mode-fast` | Enable fast mode, skipping order-level marker trees for a 2-3x speedup. Fast mode is already the default. |
| `--tree-method iqtree` | Use IQ-TREE instead of the default VeryFastTree (slower, more accurate). |
| `--sensitive` | Loosen the HMM search to `E=1e-5` and skip GA cutoffs. |
| `--contigs` | Treat each contig in an FNA as an independent genome. |

!!! tip
    Split your threads across genomes for a directory of many bins. `-t 32 -j 4` runs 4 genomes at once with 8 threads each. The [speed and accuracy guide](../how-to/tune-speed-and-accuracy.md) covers the tradeoffs.

## Run it from anywhere

Copy the wrapper into your personal `bin` and put that directory on your `PATH`.

```bash
mkdir -p "$HOME/bin"
cp gvclass-a "$HOME/bin/"
echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
source "$HOME/.bashrc"
```

Now `gvclass-a` runs from any working directory.

## Submit to the scheduler

Send the work to a compute node through a batch script. Request CPUs and pass the same number to `-t` so the allocation and the run agree.

!!! warning
    Never run a classification job on a login node. Tree building and HMM search saturate every core you give them. Submit the work and let the scheduler place it on a compute node.

```bash
#!/bin/bash
#SBATCH --job-name=gvclass
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=gvclass_%j.log

# Load the container runtime if your cluster uses modules
# module load apptainer

gvclass-a /scratch/$USER/bins /scratch/$USER/gvclass_results -t "$SLURM_CPUS_PER_TASK"
```

Save it as `gvclass.sbatch` and submit with `sbatch gvclass.sbatch`. Using `$SLURM_CPUS_PER_TASK` for `-t` keeps the thread count tied to `--cpus-per-task`, so editing one value updates both.

!!! note
    The `--cluster-type slurm`, `--cluster-queue`, `--cluster-project`, and `--cluster-walltime` flags belong to the pixi CLI (`./gvclass`), not the `gvclass-a` wrapper. Use them when GVClass itself submits work to the scheduler. See the [CLI reference](../reference/cli.md).

## Run the image yourself

For full control, pull the SIF and run it without the wrapper. Pull once from the public library, then bind your input and output directories into the container.

```bash
apptainer pull --library https://library.sylabs.io \
  gvclass_2.0.0.sif library://nelligroup-jgi/gvclass/gvclass:2.0.0

apptainer run -B /path/to/bins:/input -B /path/to/results:/output \
  gvclass_2.0.0.sif /input -o /output -t 32
```

Inside the container the query path is positional (`/input`) and the output uses `-o`, matching the pixi CLI. The wrapper does this for you and is the simpler path for most runs.

## See also

- [Tune speed and accuracy](../how-to/tune-speed-and-accuracy.md) for thread layout, fast mode, and tree method choices.
- [CLI reference](../reference/cli.md) for the full flag set, including the cluster submission options.
- [Classify bins](../how-to/classify-bins.md) for input preparation and the standard bin workflow.
- [Output reference](../reference/output.md) for the columns and files written to `RESULTS_DIR`.
