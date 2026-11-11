# Screen bins for contaminating contigs

A plain `gvclass` run screens every bin for sequence that does not belong to the
giant virus when the screen add-on is installed. The shares of base pairs flagged
appear in four columns of `gvclass_summary.tsv`, and the tables of each bin are
written to `contig_screen/<bin>/` in the output directory. The cellular share is
also the contamination estimate of the bin. `gvclass screen` runs the
same screens outside a classification run.

Two models run on each bin:

- the **cellular screen** marks bacterial, archaeal and eukaryotic sequence;
- the **giant-virus screen** marks contigs of a second giant virus.

A bin that the cellular screen scored takes its cellular share as
`estimated_contamination`, with `estimated_contamination_strategy`
`contig_screen_cellular_v1`. `screen_giant_virus_pct` and `screen_flagged_pct` are
review flags, not contamination estimates, and no call is a rule to delete a
contig. Read [How the contig screen works](../explanation/contig-screen.md) for the
methods and the limits before you act on a call. `estimated_completeness` does not
depend on the screen.

## Install the add-on

The add-on holds a reduced protein reference (2.8 GB), marker references, the fitted
models and the ESM-2 8M weights. The archive is 3.2 GB to download (48 parts on
[Zenodo](https://doi.org/10.5281/zenodo.23464780)) and 6 GiB on disk. GVClass checks
the SHA-256 of every part and of the joined archive. A failed download leaves a
previous installation unchanged. The install needs about 9 GiB free while it runs
(the 3 GiB archive and the 6 GiB add-on) and about 15 GiB when it replaces an
installed add-on.

The add-on is read from `<database>/contig_screen/`. Install it once, with the
command for your install type.

### Repository checkout with Pixi

```bash
pixi run setup-db --with-screen
```

This installs the database, then the add-on. With a database that is already
installed, install the add-on alone:

```bash
./gvclass screen --install
```

The add-on goes into the database directory, which must be writable.

### Apptainer image or standalone launcher

The standalone launcher, downloaded without a repository, always runs the image. In a
checkout, add `--apptainer` or `--image PATH_OR_URI`. The add-on goes into a host
directory that the image mounts, not into the image;
[Run with Apptainer](#run-with-apptainer) gives the details.

```bash
wget https://raw.githubusercontent.com/NeLLi-team/gvclass/main/gvclass
chmod +x gvclass
./gvclass screen --install --resource-cache-dir /scratch/gvclass-cache
```

Use a cache directory on a filesystem with 9 GiB free, and pass the same
`--resource-cache-dir` to every later run, classification runs included.

### Docker image

The Docker container runs as root, so its database folder is writable and
`screen --install` writes `contig_screen/` to `/opt/gvclass/resources/` inside the
container. The add-on does not persist there: it is gone when the container is
removed. Use the Apptainer image or a repository checkout for a lasting install.

## Screen during a classification

```bash
./gvclass bins -o results -t 8
./gvclass bins -o results -t 8 --apptainer --resource-cache-dir /scratch/gvclass-cache
```

After the queries finish, GVClass screens each nucleotide bin and adds these
columns to the end of `gvclass_summary.tsv`:

| Column | Content |
| --- | --- |
| `screen_status` | `ok`, `partial`, `failed`, `disabled`, `not_installed` or `not_applicable`. |
| `screen_cellular_pct` | Percent of base pairs called cellular. |
| `screen_giant_virus_pct` | Percent of base pairs in contigs called as another giant virus. |
| `screen_flagged_pct` | Percent of base pairs called by either screen. |

A percentage that could not be computed reads `NaN`. Missing evidence is not a clean
bin: `NaN` means the screen did not score it. [Output files](../reference/output.md#screen-columns)
lists what each status means.

When the cellular screen scored the bin (`screen_status` is `ok`, or `partial` with a
number in `screen_cellular_pct`), these columns change as well:

| Column | Value |
| --- | --- |
| `estimated_contamination` | `screen_cellular_pct`. |
| `estimated_contamination_strategy` | `contig_screen_cellular_v1`. |
| `contamination_query_support` | `supported`. |
| `contamination_model_reliability` | `simulated_benchmark`. |
| `contamination_type` | `cellular` if the cellular share is at least 10; `mixed_viral` if it is below 10 and `screen_giant_virus_pct` is at least 10; otherwise `clean`. |

For any other status the five columns keep the values of the marker model. See
[Completeness and contamination](../explanation/quality-metrics.md#contamination)
for the error of each estimate.

A bin of 176 kb takes about 2 minutes on 8 threads. DIAMOND searches against the
reference take nearly all of that time; loading the models takes seconds. The screen
runs one process per bin, as many at once as the worker plan allows.

In a run of `example/` at 8 threads, two bins (176 kb and 241 kb) were screened at
once with 4 threads each. In three runs they added 5 to 8 minutes to the run, and
peak memory was 16 GB, against 4 GB with `--no-screen`. Request memory for the number
of bins that run at once: about 16 GB for two.

A screen that runs for more than 4 hours is stopped together with the processes it
started, and the bin reads `failed`. The output of each screen is in
`contig_screen/<bin>/screen.log`.

DIAMOND 2.1.23 sometimes stops without error and never ends. A DIAMOND search that
runs for more than 30 minutes is killed and started again, up to 3 runs in all, and
`screen.log` records a warning with the run number. A search that exits with an error
is not repeated. If all 3 runs hit the limit, the screen fails and its message names
the 30-minute limit.

Without the add-on, `screen_status` is `not_installed`, the three percentages are
`NaN`, and the classification, completeness and contamination columns are the same
as before.

### Turn the screen off

```bash
./gvclass bins -o results -t 8 --no-screen
```

`screen_status` is then `disabled` and the contamination columns hold the marker-model
values. To turn it off for every run, set
`pipeline.contig_screen: false` in the [configuration file](../reference/configuration.md).

### Resume

`--resume` screens every bin of the run, also the queries that it skips for the
classification, and reuses the saved result of a bin when all of these hold:

- the folder `contig_screen/<bin>/` holds a `contig_screen_summary.tsv` with the
  totals;
- neither screen failed;
- the row names this bin, and its `input_sha256` is the SHA-256 of the bin file as it
  is now;
- its `manifest_sha256` is the digest of the add-on that is installed now.

Any other bin is screened again. A bin that reported `not_installed` is screened
once the add-on is installed, every bin is screened again after the add-on changes,
and a bin whose file changed is screened again. With `--no-screen`, `--resume` reports
`disabled` and leaves the saved results in place for a later run.

## Screen without a classification

```bash
./gvclass screen bins/ -o bins_screen -t 8
```

`bins/` holds `.fna`, `.fa` or `.fasta` files, optionally gzipped; a single file also
works. The database comes from `-d`, then `GVCLASS_DB`, then `resources`. With
`--apptainer`, `--image` or the standalone launcher the command runs in the image;
see [Run with Apptainer](#run-with-apptainer).

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

`-d /path/to/resources`, then `GVCLASS_DB`, overrides the cache. Screening mounts the
`contig_screen/` folder of that directory into the database of the image, so the
directory may hold only the add-on or a complete database; the base database of the
image is always used. The installer adds only `contig_screen/` to the directory.

Installation writes to the selected host directory without changing the SIF. The
installer reports `/screen-resources`, which is the container mount of that host
directory. Keep `--image` and any cache or database flags on installation commands.

A classification run with `--image` or `--apptainer` mounts the same add-on read-only
and screens each nucleotide bin. It takes the add-on from `GVCLASS_DB` if that is set,
else from the resource cache, as `screen` does; without one, `screen_status` is
`not_installed`. Add `--no-screen` to skip the screen.

## Read the results

For each bin `<name>`, the screen output directory holds:

| File | Content |
| --- | --- |
| `<name>.cellular.tsv` | One row per 10 kb window (5 kb step) with the stage 1 score `p_cellular` and the bin-context score `p_context`. |
| `<name>.giant_virus.tsv` | One row per contig of 5 kb or more with `p_foreign` and `called_foreign`. |
| `<name>.summary.txt` | The printed summaries of both screens, including the union of the calls. |
| `contig_screen_summary.tsv` | One row per bin with the status of each screen and the base pairs called. |

In a classification run, this directory is `contig_screen/<bin>/` in the output
directory.

In `<name>.cellular.tsv`, `called_cellular` is true when `p_context` reaches the
stage 2 threshold. The summary counts a base pair as called when the mean score of
the windows over it reaches the threshold. To list the windows called:

```bash
awk -F'\t' 'NR==1 || $NF=="True"' bins_screen/bin1.cellular.tsv | cut -f1,4,8,11,12
```

Treat a bin with 80% or more of its base pairs called cellular as a review case.
The bin may not be a giant-virus bin, or its family may be absent from training.

## Requirements

- Contigs under 3 kb are not scored. A bin with no contig of 3 kb or more has no
  cellular result: the status of that screen is `not_scored`.
- Two records with the same FASTA identifier make the bin `failed` with no numbers: the
  screen cannot total such a bin, and it does not rename or merge contigs. The reason
  is in `<name>.summary.txt`. A bin that cannot be read fails the same way, and the
  other bins of the batch still run.
- The giant-virus screen needs two contigs of 5 kb or more with genes; for a bin with
  fewer it writes no table and its status is `not_scored`.
- The cellular screen needs at least one called gene. For a bin with contigs of 3 kb
  or more but no called gene it writes no table and prints that the bin is not scored.
- The scores do not depend on the file name or the contig order. The random draw of
  proteins for the embedding depends on the sequence content of the bin. Windows or
  contigs that tie (equal distances, equal marker counts) are ordered by their
  content, and two contigs with an identical sequence share one draw.
- The module command lines (`python -m src.contig_screen.cellular` and
  `python -m src.contig_screen.giant_virus`) reject a FASTA file with two records of
  one identifier and read `.gz` in any case. A run that scores nothing removes the
  table of an earlier run at `--out`.
- `torch` and `transformers` come with the pixi environment and the image. The
  screens read the ESM-2 weights from the add-on and use no network.

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
The automatic-screen, resume and disabled-screen runs also check their combined
logs for ANSI escapes and decorative symbols, including any validation warnings.
The artifact also records host CPU information from `lscpu` and the SIF's PyTorch
version, CPU capability, thread counts, loaded thread pools and build configuration
from a fresh metadata process, before CLI per-task thread limits.
These records support runtime comparisons; they do not establish the cause of a
previous numerical difference.

The SIF classification results are also resumed through three launcher routes:
`gvclass --apptainer --image IMAGE`, `gvclass --image=IMAGE`, and the standalone
`gvclass` script outside a checkout. Each must skip all three queries without
starting scientific processing. The checks compare the TSV and CSV bytes and
archive member contents after each route, ignoring gzip and tar header metadata. The original
query execution timestamps must stay unchanged. Each route must have a fresh run
ID, complete status, matching database/software versions and stable settings
(except resume and worker allocation). Observed status, fingerprints and exit
codes are saved before assertions, alongside per-route logs, in the validation
artifact.

The job allows 180 minutes for the build, public downloads and validation. The
installation and screening step has a 120-minute limit to leave time to save
evidence if a download stalls. Before installing the add-on, the job requires
12 GiB of free disk space: the joined archive and extracted files coexist during
installation. The add-on stays in the host cache and is not included in the SIF.
