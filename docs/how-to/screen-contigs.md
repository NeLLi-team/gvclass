# Screen bins for contaminating contigs

`gvclass screen` flags sequence in a giant-virus bin that does not belong to the
giant virus. Two models run on each bin:

- the **cellular screen** marks bacterial, archaeal and eukaryotic sequence;
- the **giant-virus screen** marks contigs of a second giant virus.

A flag is a review flag, not a contamination estimate, and not a rule to delete a
contig. Read [How the contig screen works](../explanation/contig-screen.md) for the
methods and the limits before you act on a call.

## Install the screen resources

The screen needs the 2.1.0 database and an extra resource archive: a reduced protein
reference (2.8 GB), marker references, the fitted models and the ESM-2 8M weights.
The archive is 3.2 GB to download (48 parts on [Zenodo](https://doi.org/10.5281/zenodo.23464780)) and 6 GiB on disk. Install both once:

```bash
pixi run setup-db
./gvclass screen --install
```

The second command downloads the archive in parts, checks the SHA-256 of every
part and of the joined archive, and installs it as `contig_screen/` inside the
database directory. A failed download leaves a previous installation unchanged.

## Run the screen

```bash
./gvclass screen bins/ -o bins_screen -t 8
```

`bins/` holds `.fna`, `.fa` or `.fasta` files, optionally gzipped; a single file also
works. The database comes from `-d`, then `GVCLASS_DB`, then `resources`.

DIAMOND searches against the reference dominate the run time. In a test with the
research prototype, a bin of 90 kb took 106 s on four threads with the 10% reference.

## Run with Apptainer

Select a local SIF image and install the optional resources in the host cache:

```bash
./gvclass screen --install --image "$PWD/gvclass_2.1.1.sif"
./gvclass screen bins/ -o bins_screen -t 4 --image "$PWD/gvclass_2.1.1.sif"
```

With `GVCLASS_DB` unset, the installer writes to
`~/.cache/gvclass/resource-cache/v2.1.0/contig_screen/`. Use
`--resource-cache-dir /path/to/cache` on both commands to select another add-on
cache. Screening mounts that add-on read-only alongside the database inside the
image. The input is mounted read-only; results go to the host output directory.

`-d /path/to/resources`, then `GVCLASS_DB`, overrides the cache and selects a
complete host database. When installing with either setting, select a directory
that already holds the base database if you will use it for screening. The
installer adds only `contig_screen/`; it does not install the base database.
For a directory that will hold only the add-on, use `--resource-cache-dir` and
unset `GVCLASS_DB`.

Installation writes to the selected host directory without changing the SIF. The
installer reports `/screen-resources`, which is the container mount of that host
directory. Keep `--image` and any cache or database flags on installation commands.

## Read the results

For each bin `<name>`, the output directory holds:

| File | Content |
| --- | --- |
| `<name>.cellular.tsv` | One row per 10 kb window (5 kb step) with the stage 1 score `p_cellular` and the bin-context score `p_context`. |
| `<name>.giant_virus.tsv` | One row per contig of 5 kb or more with `p_foreign` and `called_foreign`. |
| `<name>.summary.txt` | The printed summaries of both screens, including the union of the calls. |
| `contig_screen_summary.tsv` | One row per bin with the status of each screen. |

In `<name>.cellular.tsv`, `called_cellular` is true when `p_context` reaches the
stage 2 threshold. The summary counts a base pair as called when the mean score of
the windows over it reaches the threshold. To list the windows called:

```bash
awk -F'\t' 'NR==1 || $NF=="True"' bins_screen/bin1.cellular.tsv | cut -f1,4,8,11,12
```

Treat a bin with 80% or more of its base pairs called cellular as a review case.
The bin may not be a giant-virus bin, or its family may be absent from training.

## Requirements

- Contigs under 3 kb are not scored. The giant-virus screen needs two contigs of 5 kb
  or more with genes; it writes no table for a bin with fewer.
- The random draw of proteins for the embedding depends on the sequence content of the
  bin, not on the file name or the contig order. Ties between equal distances inside
  the bin-context features can still follow the row order, so a reordered copy of a
  bin can change a score slightly. Two contigs with an identical sequence can change
  the draw.
- `torch` and `transformers` come with the pixi environment. The screens read the ESM-2
  weights from the resource archive and use no network.

See [`gvclass screen`](../reference/cli.md#gvclass-screen) for all options.

## Validate a SIF before release

The container CI job for tags and manual runs installs the pinned add-on through
`./gvclass screen --install --image ... --resource-cache-dir ...` and screens
`example/GVMAG-S-1096109-37.fna` through the host launcher. The gate requires both
screens to succeed, 33 cellular windows, 10 giant-virus contigs, complete input
contig coverage, and finite probabilities between zero and one. The validation
artifact holds the classification and screen tables for comparison with native
results, plus logs, a disk report and the add-on manifest. Its SIF SHA-256
identifies the tested image from that run.

Classification runs with `--plain-output --verbose`. Its combined log must contain
marker debug messages and search diagnostics, without ANSI escapes or decorative
warning, checkmark and error symbols. Scientific text such as `×` is allowed.
The artifact also records host CPU information from `lscpu` and the SIF's PyTorch
version, CPU capability, thread counts, loaded thread pools and build configuration
from a fresh metadata process, before CLI per-task thread limits.
These records support runtime comparisons; they do not establish the cause of a
previous numerical difference.

The job allows 180 minutes for the build, public downloads and validation. The
installation and screening step has a 120-minute limit to leave time to save
evidence if a download stalls. Before installing the add-on, the job requires
12 GiB of free disk space: the joined archive and extracted files coexist during
installation. The add-on stays in the host cache and is not included in the SIF.
