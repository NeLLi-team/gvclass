# Classify a directory of bins

Use one FASTA file per genome or bin. All contigs in a file are analysed together
as one query.

## 1. Prepare the input

Create a directory containing nucleotide FASTA (`.fna`) or protein FASTA (`.faa`)
files. `.fasta` and `.fas` are also accepted; GVClass infers their sequence type
from the content. Give each file a distinct name, using letters, numbers, `_`, or
`-` before the extension.

Copy the two bundled nucleotide bins to try this workflow:

```bash
mkdir -p query_genomes
cp example/*.fna query_genomes/
```

For your own analysis, put your FASTA files in `query_genomes/` instead. Files
must be directly inside the directory, not in subdirectories. The default
minimum total nucleotide length is 20,000 bp per file. Protein inputs have no
nucleotide length requirement.

## 2. Run GVClass

Run from the GVClass repository after [installation and database setup](../tutorials/getting-started.md):

```bash
pixi run gvclass query_genomes -o bin_results -t 8
```

`-t 8` sets the total thread budget. GVClass chooses how many queries to process
at once. See [thread settings](tune-speed-and-accuracy.md#threads-and-workers)
for manual control.

To change the nucleotide length requirement, add `--min-length 30000` for a
30,000 bp minimum, or `--min-length 0` to disable it. `--allow-short` accepts
shorter inputs with a warning. A sequence that passes the length check still
needs marker genes for classification.

To add eligible queries to a combined species tree, add
`--species-tree-combined`. See the [species-tree guide](build-a-species-tree.md)
for the full command and outputs.

## 3. Read the results

```bash
cut -f1,2,4 bin_results/gvclass_summary.tsv
```

The copied example contains two bins, so the summary should have two data rows.
Each row reports the query, assigned lineage, and classification confidence.
The full TSV and CSV summaries include quality estimates. Each query's files
are stored in a `<query>.tar.gz` archive inside `bin_results/`.

If a query fails, inspect `bin_results/gvclass_failed_queries.tsv` and
`bin_results/run.log`. After correcting the cause, repeat the command with
`--resume` to skip completed queries:

```bash
pixi run gvclass query_genomes -o bin_results -t 8 --resume
```

Use a new output directory when changing analysis settings. For quality
interpretation, follow [Assess genome quality](assess-genome-quality.md); for
all output paths, see the [output reference](../reference/output.md).
