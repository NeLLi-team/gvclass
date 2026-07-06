# Classify individual contigs

Most GVClass runs treat one FASTA file as one putative genome. The `--contigs` mode flips that: every sequence in a multi-contig FNA becomes its own query. Use it to screen raw assembly output before any binning has happened.

## Run it

Point GVClass at a single multi-contig FNA and add `--contigs` (short form `-C`):

```bash
pixi run gvclass --contigs metagenome_contigs.fna -o results -t 32
```

A directory of FNA files also works, but a single multi-contig file is the primary case. The combined result lands in `results/gvclass_summary.tsv`, with one row per retained contig.

## What happens to your input

In `--contigs` mode GVClass:

1. Splits each sequence in the input into its own file.
2. Sanitises contig IDs for filenames (spaces and `/\:*?"<>|` become `_`).
3. Drops contigs shorter than the contig-mode minimum length (default `10000` bp).
4. Classifies each retained contig independently through the full pipeline (see [How GVClass works](../explanation/how-it-works.md)).
5. Combines every per-contig result into `gvclass_summary.tsv`, then removes the temporary split files.

Each contig is classified in isolation, so the taxonomy, completeness, and contamination on a row describe that one sequence, not the assembly as a whole.

## When to use it

Reach for `--contigs` when:

- You want to screen metagenome contigs before binning.
- You are evaluating candidate viral contigs as standalone genomes.
- Bins are not ready and you need rapid triage.

For finished bins, use the default directory workflow in [Classify a directory of bins](classify-bins.md) instead. That workflow keeps all contigs of one genome together as a single query, which gives a genome-level completeness and contamination call rather than a per-contig one.

## Tune the length filter for giant viruses

The default 10 kb cutoff is generous. Giant-virus genomes are large, and short contigs rarely carry enough markers for a confident placement. For one run, pass `--contigs-min-length`:

```bash
pixi run gvclass --contigs metagenome_contigs.fna -o results -t 32 --contigs-min-length 50000
```

For a persistent default, set `pipeline.contigs_min_length` to at least `30000`, preferably `50000`:

```yaml
pipeline:
  contigs_min_length: 50000
```

See [Configure the database](configure-the-database.md) and the [configuration reference](../reference/configuration.md) for where this key lives and how CLI flags, config, and defaults take precedence.

The same contig-mode minimum is used when validating the split contig files. If you set `--contigs-min-length 5000`, retained 5 kb contigs are accepted by the length gate. Use `--contigs-min-length 0` to remove the contig split filter for exploratory runs.

The full flag list, including `-C`, `--contigs-min-length`, `--allow-short`, and the thread options, is in the [CLI reference](../reference/cli.md).
