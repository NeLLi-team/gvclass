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
