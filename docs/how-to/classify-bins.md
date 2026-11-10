# Classify a directory of bins

Use one FASTA file per genome or bin. All contigs in a file are analysed together
as one query.

## 1. Prepare the input

Create a directory containing nucleotide FASTA (`.fna`) or protein FASTA (`.faa`)
files. `.fasta`, `.fas` and `.fa` are also accepted; GVClass infers their
sequence type from the content. Extensions are matched without regard to case.
Give each file a distinct name, using letters, numbers, `_`, or `-` before the
extension. The names `gvclass_summary`, `gvclass_summary.extended`,
`gvclass_failed_queries`, `run_status`, `run`, `results`, `species_tree` and
`_combined` are reserved for GVClass outputs; GVClass stops and asks you to
rename a file that uses one.

Sequence IDs inside a file must be unique. GVClass renames each record to
`<file name>|<original ID>`, with every `|` and `:` in the original ID replaced
by `__` and empty fields at either end dropped. For example, `binA|contig_1` in
`q.fna` becomes `q|binA__contig_1`, and `gi|9626243|ref|NC_001416.1|` becomes
`q|gi__9626243__ref__NC_001416.1`. GVClass output used as input keeps its IDs.
If two IDs in a file become the same name, or a file has a repeated ID, that file
is rejected with a message that names the IDs.

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

A file that fails validation (for example an empty file, a repeated sequence ID,
or a total length below the minimum) does not stop the other files. It is
reported as a failed query with status `invalid_input`, the remaining files are
processed, and the run ends with a non-zero exit code. If no file is valid, the
run stops before it processes anything.

If a query fails, inspect `bin_results/gvclass_failed_queries.tsv` and
`bin_results/run.log`. After correcting the cause, repeat the command with
`--resume` to skip completed queries:

```bash
pixi run gvclass query_genomes -o bin_results -t 8 --resume
```

`--resume` reads `bin_results/results.jsonl`, the record of finished queries. It
skips a query only if its input file, the settings that change results, the database
version and the GVClass version equal those of the saved result and its archive is
still there; otherwise it runs the query again and the log says which field differed.
Changing only the number of threads does not stop a skip. For quality
interpretation, follow [Assess genome quality](assess-genome-quality.md); for
all output paths, see the [output reference](../reference/output.md).
