# Classify a directory of bins

Metagenomic binning gives you one FASTA per putative genome. GVClass takes a directory of those bins and returns a taxonomy call and a quality table for each one. This is the recommended way to run the tool.

## Prepare the input directory

Put one FASTA file per putative genome in a single directory. A file may hold several contigs that belong to the same genome; GVClass treats the whole file as one query. Use nucleotide FASTA (`.fna`) or protein FASTA (`.faa`).

For reliable giant-virus calls, mind the assembled length and gene content:

- Default minimum: 20 kb of total assembled sequence per nucleotide file. GVClass rejects shorter `.fna` inputs unless you lower the floor with `--min-length`, set `quality.min_length` in the config, or pass `--allow-short`.
- Better reliability: at or above 30 kb.
- Preferred: at or above 50 kb.

Length is really a proxy for gene content. GVClass infers taxonomy by placing marker genes in reference trees, so a query needs to carry several genes for the method to work. A short fragment with only one or two predicted proteins rarely hits enough markers, and GVClass cannot assign a taxonomy when no markers are found.

Keep filenames clean. The filename becomes the query name, so avoid `.`, `;`, and `:`. Use `_` or `-` instead. For protein input, write headers as `filename|proteinid`.

!!! tip
    Filter short contigs (below a few kb) out of each bin before you run. For giant viruses, prefer bins assembled to at least 50 kb; short fragments add noise and weaken the marker signal.

To adjust the nucleotide length gate for a run, pass `--min-length`:

```bash
pixi run gvclass my_bins -o my_results -t 32 --min-length 30000
```

For a persistent default, set the same value in the config:

```yaml
quality:
  min_length: 30000
```

Use `--min-length 0` only when you want to keep normal FASTA validation but remove the nucleotide length floor. Use `--allow-short` when you want to bypass the length gate for exploratory runs while still seeing a warning.

## Run the classification

Run from the cloned repository so the launcher can find `src/`:

```bash
pixi run gvclass my_bins -o my_results -t 32
```

Here `my_bins` is your input directory, `-o my_results` is the output directory, and `-t 32` sets the total thread budget. Omit `-o` and results go to `<query_dir>_results` (for this command, `my_bins_results`).

!!! note
    To run from any directory, install the CLI wrapper or use the Apptainer wrapper on a cluster. See [Configure the database](configure-the-database.md) for shared database setups and [Run on an HPC cluster](run-on-hpc.md) for batch submission.

## Choose parallelism

Two flags control throughput:

- `-t` sets the total number of threads.
- `-j` sets how many bins run at once (workers). GVClass picks a worker count automatically when you leave `-j` unset.

For a directory of many small bins, more workers help; for a few large genomes, give each worker more threads. See [Tune speed and accuracy](tune-speed-and-accuracy.md) for the trade-offs.

## Read the results

GVClass writes a combined table for the whole run plus per-query files:

- `gvclass_summary.tsv` and `gvclass_summary.csv`: one row per bin with the taxonomy call and quality metrics.
- `<query>.tar.gz`: per-query artifacts, including the per-query summary rows.
- `gvclass_summary.extended.tar.gz`: archived extended diagnostics.

See [Output files and columns](../reference/output.md) for every file and the full column layout. To turn the table into a curation decision, read [Assess genome quality](assess-genome-quality.md). For what happens between FASTA and taxonomy, see [How GVClass works](../explanation/how-it-works.md).
