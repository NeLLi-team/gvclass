# Getting started

Classify the three bundled example genomes and inspect their summary. You need
Linux x86-64 (`linux-64`) and Git.

## 1. Install Pixi

Install Pixi using its [installation instructions](https://pixi.sh/latest/installation/).
Open a new shell, then check that it is available:

```bash
pixi --version
```

## 2. Install GVClass

```bash
git clone https://github.com/NeLLi-team/gvclass.git
cd gvclass
pixi install --frozen
```

Run the remaining commands from this directory. Pixi installs Python and the
other dependencies locally.

## 3. Download the reference database

```bash
pixi run setup-db
```

The default download is about 1.5 GB. Setup checks the archive checksum and
extracts it into `resources/`. Later runs reuse that directory. For another
location, follow [Configure the database](../how-to/configure-the-database.md).

## 4. Run the example

```bash
pixi run gvclass example -o example_results -t 8
```

The `example/` directory contains two nucleotide bins (`.fna`) and one protein
set (`.faa`). Each file is one query. GVClass uses eight threads as its total
budget and writes results to `example_results/`.

## 5. Inspect the summary

Print the query names, assigned lineages, and classification confidence:

```bash
cut -f1,2,4 example_results/gvclass_summary.tsv
```

The table should contain these three queries:

```text
AC3300027503___Ga0255182_1000024
GVMAG-S-1096109-37
PkV-RF01
```

Check the number of result rows:

```bash
awk 'END {print NR - 1, "queries"}' example_results/gvclass_summary.tsv
```

Expected output:

```text
3 queries
```

`taxonomy_majority` contains the lineage supported by the marker gene trees.
`taxonomy_confidence` describes its marker support. For completeness and
contamination estimates, see [Assess genome quality](../how-to/assess-genome-quality.md).

Open `example_results/gvclass_summary.csv` in a spreadsheet for the remaining
columns. `PkV-RF01` is a protein input, so its `ttable` value is `no_fna`.

To classify your own genomes, follow [Classify a directory of bins](../how-to/classify-bins.md).
To put queries together in a tree, follow [Build a species tree](../how-to/build-a-species-tree.md).
All files and columns are listed in the [output reference](../reference/output.md).
