# Run on an HPC cluster

Submit GVClass as one batch job. It processes queries within that allocation; it does not submit separate scheduler jobs.

## Run with Pixi and Slurm

### 1. Prepare the repository and inputs

Complete [installation and database setup](../tutorials/getting-started.md). Run the following commands from the GVClass repository. For a test, use the bundled nucleotide assemblies:

```bash
mkdir -p query_genomes
cp example/*.fna query_genomes/
```

### 2. Create the batch script

```bash
cat > gvclass.sbatch <<'BATCH'
#!/usr/bin/env bash
#SBATCH --job-name=gvclass
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=gvclass_%j.log

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
pixi run gvclass query_genomes -o hpc_results -t "$SLURM_CPUS_PER_TASK"
BATCH
```

The memory and time requests are starting values; adjust them for your inputs and cluster. Make `pixi` available in the batch environment.

### 3. Submit and inspect the output

Add your cluster's account, partition and QOS options to the submission command where required.

```bash
sbatch --chdir="$PWD" gvclass.sbatch
```

After the job completes:

```bash
head -n 4 hpc_results/gvclass_summary.tsv
cat hpc_results/run.log
```

To include a shared species tree, add `--species-tree-combined` to the GVClass line in the script. See [Build a species tree](build-a-species-tree.md).

The CLI options `--cluster-type`, `--cluster-queue`, `--cluster-project` and `--cluster-walltime` do not submit jobs. Configure the scheduler through the batch script and `sbatch`.

## Use the Apptainer image

The image includes GVClass, its dependencies and database. The 2.1.0 image is prepared but not yet published: publication is held until the v2.1.0 candidate resources pass the release checks. The published image will carry the approved resources. The `gvclass` entrypoint requires Apptainer and Python 3.10 or later on the host when used in container mode. If your cluster uses environment modules, load Apptainer before running `gvclass`.

### 1. Download the launcher and image

If you do not have a repository checkout on the cluster, download the standalone `gvclass` launcher and pull the SIF image:

```bash
wget https://raw.githubusercontent.com/NeLLi-team/gvclass/main/gvclass
chmod +x gvclass
apptainer pull --library https://library.sylabs.io \
  gvclass_2.1.0.sif library://nelligroup-jgi/gvclass/gvclass:2.1.0
```

Pass the downloaded image to `gvclass` with `--image` (or pass `--apptainer` to pull and cache the default image automatically).

### 2. Run within a compute allocation

Put one `.fna` or `.faa` file per genome in `query_genomes/`, then run:

```bash
./gvclass query_genomes -o apptainer_results -t 8 --image "$PWD/gvclass_2.1.0.sif"
```

For a Slurm job, replace the `pixi run gvclass` line in the batch script with:

```bash
./gvclass query_genomes -o apptainer_results -t "$SLURM_CPUS_PER_TASK" \
  --image "$PWD/gvclass_2.1.0.sif"
```

The container mode accepts the output with `-o` / `--output-dir` (or as a second positional path).

### 3. Set the cache location if needed

During release review, the database cache stays at the public v2.0.0 path, `~/.cache/gvclass/resource-cache/v2.0.0/`. To put it on your allocated scratch space, set `--resource-cache-dir` to that directory.

| Option | Purpose |
| --- | --- |
| `-t N` | Total threads; match the CPU allocation |
| `-j N` | Number of queries processed concurrently |
| `--contigs` | Classify input contigs separately |
| `--apptainer` | Run inside the Apptainer SIF image |
| `--image PATH_OR_URI` | Select a local SIF or another image URI (implies `--apptainer`) |
| `--resource-cache-dir PATH` | Select a writable database cache directory |

All standard pipeline, quality-model and species-tree flags (`--quality-models`, `--species-tree`, `--species-tree-combined`, `--resume`, etc.) are forwarded into the container; `--config`, `--database` and the cluster options are native-only. See the [CLI reference](../reference/cli.md) for supported options.
