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

The image includes GVClass, its dependencies and database. The wrapper requires Apptainer and Python 3.10 or later on the host. If your cluster uses environment modules, load Apptainer before running the wrapper.

### 1. Download the wrapper and image

```bash
wget https://raw.githubusercontent.com/NeLLi-team/gvclass/main/gvclass-a
chmod +x gvclass-a
apptainer pull --library https://library.sylabs.io \
  gvclass_2.0.3.sif library://nelligroup-jgi/gvclass/gvclass:2.0.3
```

Pass the downloaded image to the wrapper with `--image` in each run.

### 2. Run within a compute allocation

Put one `.fna` or `.faa` file per genome in `query_genomes/`, then run:

```bash
./gvclass-a query_genomes apptainer_results -t 8 --image "$PWD/gvclass_2.0.3.sif"
```

For a Slurm job, replace the `pixi run gvclass` line in the batch script with:

```bash
./gvclass-a query_genomes apptainer_results -t "$SLURM_CPUS_PER_TASK" \
  --image "$PWD/gvclass_2.0.3.sif"
```

The wrapper accepts the output as a second positional path or with `-o`.

### 3. Set the cache location if needed

The database cache is stored in `~/.cache/gvclass/resource-cache/v2.0.0/`. To put it on your allocated scratch space, set `--resource-cache-dir` to that directory.

| Option | Purpose |
| --- | --- |
| `-t N` | Total threads; match the CPU allocation |
| `-j N` | Number of queries processed concurrently |
| `--contigs` | Classify input contigs separately |
| `--image PATH_OR_URI` | Select a local SIF or another image URI |
| `--resource-cache-dir PATH` | Select a writable database cache directory |

For species trees, use the Pixi steps above; the wrapper does not accept species-tree flags. See the [CLI reference](../reference/cli.md) for supported options.
