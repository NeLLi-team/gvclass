# Command-line interface

Run commands from the GVClass repository directory. Each file in `query_genomes/` is treated as one query:

```bash
pixi run gvclass query_genomes -o results -t 8
```

CLI flags override [configuration values](configuration.md), which override built-in defaults. Defaults below refer to the shipped configuration.

The command resolves the options and the configuration file into one run configuration, writes it to `<output>/run_config.json` and starts the runner on that file. The runner writes its log to `<output>/gvclass.log`. The exit status is 0 only if the runner exits 0 and the output validation finds no issue; the validation issues are printed before a non-zero exit. In the default mode the command keeps the last 200 lines of runner output (progress lines excluded) and prints them when the run fails; `--verbose` prints the runner output as it comes. Ctrl-C or SIGTERM stops the runner and its child processes and removes the temporary contig files.

The `gvclass` script starts itself again in the Python of the repository (`.pixi/envs/default/bin/python`), or under `pixi run` when that environment does not exist, so the Python on your `PATH` needs no GVClass packages. It sets `GVCLASS_RUNTIME_SELECTED=1` for the new process so that this happens once. The same applies to the copy that `pixi run install-cli` puts in `~/bin`.

## Input and output

| Option | Argument/default | Description |
| --- | --- | --- |
| `query_dir` | path (positional) | Directory of `.fna`, `.faa`, `.fa`, `.fasta`, or `.fas` files. A single FASTA file is accepted with `--contigs`. |
| `-o`, `--output-dir` | path; `<name>_results` in the current directory | Output directory. `<name>` is the name of the input directory, or the file name without its suffix for a single file. Native and container runs use the same default. |
| `-c`, `--config` | path; none | YAML configuration file. A path that does not exist is an error (exit status 2). Without `-c`, GVClass reads `gvclass_config.yaml` in the working directory, then `config/gvclass_config.yaml` in the repository, then uses the built-in defaults. See [configuration](configuration.md). |
| `-d`, `--database` | path | Database directory. Overrides `GVCLASS_DB` and `database.path`. |
| `--resource-cache-dir` | path; `GVCLASS_RESOURCE_CACHE`, then `database.cache_path`, else `<database>/.gvclass_cache` | Native Parquet materialization cache. Precedence: flag, environment variable, configuration, default. A relative configured path is resolved under the database directory. If the selected directory is not writable, use `<system temp>/gvclass_resource_cache_<uid>/<database-path hash>` (mode 0700, owned by the current user). |

## Threads and workers

| Option | Argument/default | Description |
| --- | --- | --- |
| `-t`, `--threads` | integer of at least 1; `4` | Total thread budget. |
| `-j`, `--max-workers` | integer of at least 1; automatic | Maximum queries processed in parallel. |
| `--threads-per-worker` | integer of at least 1; automatic | Threads allocated to each query. |

The worker count never exceeds the number of queries to process, and the threads are divided after that cap: `-j 8 -t 16` with two queries runs two workers with eight threads each. The startup banner shows this plan. The numeric libraries (OpenMP, OpenBLAS, MKL, NumExpr) get half of the threads of a worker, unless the environment already sets a smaller integer.

## Classification

| Option | Argument/default | Description |
| --- | --- | --- |
| `--tree-method` | `veryfasttree`, `iqtree`, `fasttree`; `veryfasttree` | Tree builder. `fasttree` is an alias for `veryfasttree`. |
| `--iqtree-mode` | `fast`, `ufboot`; `fast` | IQ-TREE search mode for species trees. Per-marker trees use `--fast`. |
| `--mode-fast` | flag; enabled by default | Skip order-level marker trees; order-level HMM searches still run. |
| `-e`, `--extended` | flag | Include order-level marker trees. Overrides `--mode-fast` regardless of argument order, natively and in container mode. |
| `--sensitive` | flag; enabled by default | Use `E=1e-5` and `domE=1e-5` instead of GA cutoffs. |
| `--no-sensitive` | flag | Disable sensitive HMM mode and use model GA cutoffs when available. Overrides `--sensitive` and configuration regardless of argument order. |
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
| `--resume` | flag | Skip a query only if `results.jsonl` holds a complete record whose input file name and hash, result settings, database version and GVClass version equal the current ones and whose archive is still in place. Other queries run again, and the log names the field that differed. See [classify bins](../how-to/classify-bins.md). |
| `--plain-output` | flag | Disable emojis and ANSI colors. Also set by `GVCLASS_PLAIN_OUTPUT=1`. Low-level search, tree, output-validation and post-processing diagnostics use plain text in all modes. |

## Other options

| Option | Description |
| --- | --- |
| `-h`, `--help` | Show command-line help. |
| `-v`, `--verbose` | Enable verbose logging and print the runner output as it comes. |
| `--version` | Report the software version and installed database version. |

```bash
pixi run gvclass --help
pixi run gvclass --version
```

`--version` reports software version 2.1.1 and the installed database version. Database setup installs v2.1.0.

## Contig screen

| Option | Argument/default | Description |
| --- | --- | --- |
| `--no-screen` | flag | Skip the contig screen. `screen_status` is `disabled`, and the contamination columns hold the marker-model values. Overrides `pipeline.contig_screen`. With `--resume`, saved screen results stay in place. |

The screen runs after the queries when the add-on is installed. See
[Screen bins for contaminating contigs](../how-to/screen-contigs.md).

## `gvclass screen`

```bash
pixi run gvclass screen <bin.fna | folder> [-o DIR] [-d DATABASE] [-t THREADS]
pixi run gvclass screen --install [-d DATABASE]
```

Runs the cellular screen and the giant-virus screen on every bin. Add `--image`
or `--apptainer` to use the container. See
[Screen bins for contaminating contigs](../how-to/screen-contigs.md).

| Option | Argument/default | Description |
| --- | --- | --- |
| `input` | FASTA file or folder | A bin, or a folder of `.fna`, `.fa`, `.fasta` or `.fas` files (also gzipped). |
| `-o`, `--output-dir` | `<input>_screen` | Directory for the tables and summaries. |
| `-d`, `--database` | `GVCLASS_DB`, then `resources` | Database directory that holds `contig_screen/`. |
| `-t`, `--threads` | 4 | Threads for gene calling, DIAMOND and the embedding model. |
| `--install` | flag | Download and install `contig_screen/` into the database, then exit. |

The exit code is 1 when no FASTA file is found, when the screen resources are not
installed, or when a screen fails for a bin; the other bins still run.

In container mode, `--install` writes the add-on to the host resource cache when
neither `-d` nor `GVCLASS_DB` is set. Screening mounts the `contig_screen/` folder of
that cache into the database of the image. `-d`, then `GVCLASS_DB`, overrides the
cache: its `contig_screen/` folder is mounted the same way, whether the directory holds
only the add-on or a complete database, and the database of the image stays in place.
Installation adds only `contig_screen/` to the selected directory.
The container reads inputs and resources through read-only mounts and writes to
the selected output directory. The thread count (`-t`) goes to the container as typed,
so the default of four threads and the `<input>_screen` output are those of a native run.

## `setup-db`

```bash
pixi run setup-db [-c CONFIG] [-d DATABASE] [--with-screen] [--update]
```

Downloads and installs the database when it is missing or incomplete. `-c` and `-d`
select the configuration file and the database directory, as for `gvclass`.

| Option | Description |
| --- | --- |
| `--update` | Replace a complete database with the pinned source, whatever its installed version. Without this option, setup asks before it replaces a complete database that is older than the pinned source (the default answer is No), and keeps it when no terminal is available. The contig-screen add-on is kept. |
| `--with-screen` | Also install the contig-screen add-on, as `gvclass screen --install` does. |

## Apptainer container mode

`gvclass` supports running inside the published Apptainer SIF container via `--apptainer` or `--image` (or automatically when `gvclass` is downloaded as a standalone script without a repository checkout). The launcher reads only the options that it needs for the mounts: the input path, `-o`/`--output-dir`, `-C`/`--contigs`, `--image`, `--resource-cache-dir`, and for `screen` also `-d` and `--install`. It passes every other option to GVClass in the container unchanged and sets no default of its own, so the thread count, the output directory (`<name>_results` in the current directory) and all other defaults are those of a native run. A missing input path is an error (exit code 2).

| Option | Argument/default | Description |
| --- | --- | --- |
| `--apptainer` | flag | Execute inside the Apptainer SIF image (`oras://ghcr.io/nelli-team/gvclass:2.1.1` by default). |
| `--image` | path or address | Local `.sif` path, or an `oras://`, `docker://`, `library://` or `https://` address (implies `--apptainer`). A `library://` address needs a library endpoint that Apptainer knows. |
| `--resource-cache-dir` | path; `~/.cache/gvclass/resource-cache/v2.1.0` | Classification: writable materialization cache mounted at `/resource-cache`. Screen: host add-on cache used when neither `-d` nor `GVCLASS_DB` is set; writable during installation and mounted read-only during screening. |

A plain `http://` image address is rejected; use `https://`. Remote images are cached under `~/.cache/gvclass/images/` as regular files, in a
file whose name starts with a digest of the whole address, so two addresses never
share a file; symlink cache entries are rejected. Before reuse, the launcher checks the SIF
object table, primary system partition and file bounds.
Downloads use a separate temporary path and enter the cache only after the pull
and those checks succeed. Invalid existing cache entries are left in place; move
the path named in the error before retrying. These checks detect malformed or
truncated files; they do not authenticate the publisher or verify every payload
byte. The cache filesystem must support hard links for atomic installation.
An explicit local image path is passed directly to Apptainer.

All other pipeline, species-tree and resume options (`-t`, `-j`, `--threads-per-worker`, `--tree-method`, `--iqtree-mode`, `--mode-fast`, `-e`/`--extended`, `--sensitive`, `--no-sensitive`, `--completeness-mode`, `--species-tree`, `--species-tree-combined`, `--species-tree-trim`, `--min-length`, `--contigs-min-length`, `--allow-short`, `--no-screen`, `--resume`, `--plain-output`, `-v`/`--verbose`) go to the container as typed. For classification, `--config` and `--database` are native-only and the launcher rejects them: the container uses the database of the image. When both `--sensitive` and `--no-sensitive` are given, `--no-sensitive` wins in both modes. A classification run also mounts the add-on read-only into the database of the image, from the
same place that `screen` uses (`GVCLASS_DB` if set, else the resource cache), so that the contig
screen runs in the container. Without an installed add-on, `screen_status` is `not_installed`.
Inside the container the database is always the one in the image, whatever `GVCLASS_DB` the
host has: that variable selects only where the add-on comes from.
See [run on HPC](../how-to/run-on-hpc.md) for container commands.
