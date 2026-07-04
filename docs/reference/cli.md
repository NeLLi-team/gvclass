# Command-line interface

GVClass runs from the repository through pixi, or from the `gvclass-a` Apptainer wrapper. The entry point is `src.bin.gvclass_cli:main`, launched by the `./gvclass` script.

```bash
pixi run gvclass QUERY_DIR -o OUTPUT_DIR -t THREADS [options]
```

Run this from the repository directory so the `./gvclass` launcher can resolve `src/`. The `gvclass-a` Apptainer wrapper takes the input as the first positional argument and the output directory either as an optional second positional argument or with `-o`/`--output-dir` (default `<query>_results`). It exposes `-t`/`--threads`, `--tree-method`, `--mode-fast`, `-e`/`--extended`, `-j`/`--max-workers`, `--sensitive`, `--contigs`, `--image` (override the default published SIF), and `--resource-cache-dir` (host cache for compact Parquet materialization).

CLI arguments override config keys, which override built-in defaults.

## Input and output

Config keys are listed in [configuration](configuration.md); database download and location are covered in [configure the database](../how-to/configure-the-database.md).

| Option | Argument/Default | Description |
| --- | --- | --- |
| `query_dir` | path (positional) | Directory with `.fna`/`.faa` inputs (also `.fasta`/`.fas`), or a single FASTA file when `--contigs` is set. |
| `-o`, `--output-dir` | path; default `<query_dir>_results` | Output directory. |
| `-c`, `--config` | path; default `config/gvclass_config.yaml` | Configuration file. |
| `-d`, `--database` | path | Database path. Overrides `GVCLASS_DB` and config. |

## Parallelism

| Option | Argument/Default | Description |
| --- | --- | --- |
| `-t`, `--threads` | int; default `4` (config) | Total threads. Overrides config. |
| `-j`, `--max-workers` | int; auto if unset | Parallel workers. |
| `--threads-per-worker` | int; auto if unset | Threads per worker. |

## Pipeline and markers

For choosing between these settings, see [tune speed and accuracy](../how-to/tune-speed-and-accuracy.md).

| Option | Argument/Default | Description |
| --- | --- | --- |
| `--tree-method` | `{veryfasttree,iqtree,fasttree}`; default `veryfasttree` | Tree builder. `fasttree` is an alias for `veryfasttree`. |
| `--iqtree-mode` | `{fast,ufboot}`; default `fast` | Species-tree IQ-TREE search mode. Per-marker trees always use `--fast`. |
| `--mode-fast` | flag; on by default | Enable fast mode. |
| `-e`, `--extended` | flag | Build trees for all markers. Turns fast mode off (slower). |
| `--sensitive` | flag; on by default | Force `E=1e-5`, `domE=1e-5`; skip GA cutoffs. |
| `--completeness-mode` | `{legacy,novelty-aware}`; default `novelty-aware` | Completeness estimator surfaced as the primary estimate. |

## Species tree

See [build a species tree](../how-to/build-a-species-tree.md) and [the species tree](../explanation/species-tree.md).

| Option | Argument/Default | Description |
| --- | --- | --- |
| `--species-tree` | flag | Build one supermatrix species tree per query. Writes `species_tree/<query>/` and fills the four `species_tree_*` summary columns. |
| `--species-tree-combined` | flag | Implies `--species-tree`. Also builds one combined tree over all queries (`species_tree/combined.*`). |
| `--species-tree-trim` | `{witchi,pytrimal,none}`; default `witchi` | Supermatrix column trimming. |

## Inputs handling

| Option | Argument/Default | Description |
| --- | --- | --- |
| `-C`, `--contigs` | flag | Split FNA into per-contig queries. Accepts a file or a directory. |
| `--allow-short` | flag | Accept FNA shorter than the 20 kb minimum. |
| `--resume` | flag | Skip completed queries recorded in `run_status.json`. Older `.SUCCESS` sentinels and summary+archive pairs are still accepted for existing output directories. |
| `--plain-output` | flag | Disable emojis and ANSI colors. Also env `GVCLASS_PLAIN_OUTPUT=1`. |

## Cluster

See [run on HPC](../how-to/run-on-hpc.md).

!!! note
    These flags are accepted but currently have no effect: GVClass runs locally and parallelizes across the cores of the machine or allocation it is launched on. To use a scheduler, submit a GVClass run as a single batch job, as shown in [run on HPC](../how-to/run-on-hpc.md).

| Option | Argument/Default | Description |
| --- | --- | --- |
| `--cluster-type` | `{local,slurm,pbs,sge}`; default `local` | Intended scheduler; parsed but currently inert. |
| `--cluster-queue` | string | Queue/partition for cluster jobs. |
| `--cluster-project` | string | Project/account for cluster billing. |
| `--cluster-walltime` | string; default `04:00:00` | Walltime for cluster jobs. |

## Other

| Option | Argument/Default | Description |
| --- | --- | --- |
| `-v`, `--verbose` | flag | Verbose output. |
| `--version` | flag | Print version information and exit. |

## Notes

!!! note "Fast and sensitive modes default on"
    Fast mode (`mode_fast: true`) and sensitive mode (`sensitive_mode: true`) are both on by default. `-e`/`--extended` turns fast mode off and builds trees for all markers.

!!! note "Database path precedence"
    The database path resolves in this order: `--database`, then the `GVCLASS_DB` environment variable, then `database.path` in the config, then the default `<repo>/resources`.

!!! note "Compact Parquet resources"
    Compact bundles can store labels and reference proteins under `parquet/`. GVClass materializes the needed TSV and marker FASTA views into `<database.path>/.gvclass_cache/` by default. Set `database.cache_path` or `GVCLASS_RESOURCE_CACHE` to choose another cache directory; the environment variable has precedence. The Apptainer wrapper sets `GVCLASS_RESOURCE_CACHE=/tmp/gvclass-resource-cache` inside the container and bind-mounts that path from a host cache directory.

## Version output

```text
GVClass Pipeline
  Software version: v2.0.0
  Database version: v2.0.0
```

The database version is read from the installed bundle. The current public
setup-download archive installs `v2.0.0`.

## See also

- [Configuration](configuration.md)
- [Output files](output.md)
- [Marker panels](markers.md)
- [Getting started](../tutorials/getting-started.md)
- [Classify contigs](../how-to/classify-contigs.md)
