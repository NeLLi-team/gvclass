# Command-line interface

Run commands from the GVClass repository directory. Each file in `query_genomes/` is treated as one query:

```bash
pixi run gvclass query_genomes -o results -t 8
```

CLI flags override [configuration values](configuration.md), which override built-in defaults. Defaults below refer to the shipped configuration.

## Input and output

| Option | Argument/default | Description |
| --- | --- | --- |
| `query_dir` | path (positional) | Directory of `.fna`, `.faa`, `.fasta`, or `.fas` files. A single FASTA file is accepted with `--contigs`. |
| `-o`, `--output-dir` | path; `<query_dir>_results` | Output directory. |
| `-c`, `--config` | path; `config/gvclass_config.yaml` | YAML configuration file. |
| `-d`, `--database` | path | Database directory. Overrides `GVCLASS_DB` and `database.path`. |

## Threads and workers

| Option | Argument/default | Description |
| --- | --- | --- |
| `-t`, `--threads` | integer; `4` | Total thread budget. |
| `-j`, `--max-workers` | integer; automatic | Maximum queries processed in parallel. |
| `--threads-per-worker` | integer; automatic | Threads allocated to each query. |

## Classification

| Option | Argument/default | Description |
| --- | --- | --- |
| `--tree-method` | `veryfasttree`, `iqtree`, `fasttree`; `veryfasttree` | Tree builder. `fasttree` is an alias for `veryfasttree`. |
| `--iqtree-mode` | `fast`, `ufboot`; `fast` | IQ-TREE search mode for species trees. Per-marker trees use `--fast`. |
| `--mode-fast` | flag; enabled by default | Skip order-level marker trees; order-level HMM searches still run. |
| `-e`, `--extended` | flag | Include order-level marker trees. Overrides fast mode. |
| `--sensitive` | flag; enabled by default | Use `E=1e-5` and `domE=1e-5` instead of GA cutoffs. |
| `--completeness-mode` | `legacy`, `novelty-aware`; `novelty-aware` | Estimator used for `estimated_completeness`. |

See [speed and tree settings](../how-to/tune-speed-and-accuracy.md) for examples.

## Species trees

```bash
pixi run gvclass query_genomes -o species_tree_results -t 8 --species-tree-combined
```

| Option | Argument/default | Description |
| --- | --- | --- |
| `--species-tree` | flag; disabled | Build a concatenated-marker tree for each eligible query. Write `species_tree/<query>/` and four `species_tree_*` summary columns. |
| `--species-tree-combined` | flag; disabled | Enable `--species-tree` and add a combined tree of eligible queries and references for each viral panel. |
| `--species-tree-trim` | `witchi`, `pytrimal`, `none`; `witchi` | Alignment trimming method before marker concatenation. |

NCLDV, PPV, and MIRUS queries use separate panels. A combined tree contains only eligible queries processed in the current run; queries skipped by `--resume` are excluded. See [build a species tree](../how-to/build-a-species-tree.md) for marker requirements and output paths.

## Input handling and resume

| Option | Argument/default | Description |
| --- | --- | --- |
| `-C`, `--contigs` | flag | Split nucleotide FASTA files into one query per contig. Accepts a file or directory. |
| `--min-length` | integer; `20000` bp | Minimum total nucleotide length per input in bin/MAG mode. `0` disables the length threshold. |
| `--contigs-min-length` | integer; `10000` bp | Minimum contig length in `--contigs` mode. |
| `--allow-short` | flag | Bypass the nucleotide length threshold; retain other FASTA validation. |
| `--resume` | flag | Skip completed queries using `run_status.json`. Also accepts `.SUCCESS` files and valid summary/archive pairs from older runs. |
| `--plain-output` | flag | Disable emojis and ANSI colors. Also set by `GVCLASS_PLAIN_OUTPUT=1`. |

## Cluster options

GVClass runs on the machine where it is launched. The following options are accepted but do not submit scheduler jobs. See [run on HPC](../how-to/run-on-hpc.md) for a batch script.

| Option | Argument/default |
| --- | --- |
| `--cluster-type` | `local`, `slurm`, `pbs`, `sge`; `local` |
| `--cluster-queue` | queue or partition name |
| `--cluster-project` | account name |
| `--cluster-walltime` | time; `04:00:00` |

## Other options

| Option | Description |
| --- | --- |
| `-h`, `--help` | Show command-line help. |
| `-v`, `--verbose` | Enable verbose logging. |
| `--version` | Report the software version and installed database version. |

```bash
pixi run gvclass --help
pixi run gvclass --version
```

## Apptainer wrapper

`gvclass-a` accepts an input path and either a second positional output path or `-o`/`--output-dir`. Its default output is `<query>_results`; its default thread count is `16`.

Supported options are `-t`, `-j`, `--tree-method`, `--mode-fast`, `-e`/`--extended`, `--sensitive`, `-C`/`--contigs`, `--min-length`, and `--contigs-min-length`. `--image` selects an Apptainer image; `--resource-cache-dir` sets the host resource cache. The wrapper does not expose species-tree options. See [run on HPC](../how-to/run-on-hpc.md) for container commands.
