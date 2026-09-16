# Classify individual contigs

Use `--contigs` to classify each sequence in a nucleotide FASTA file separately.
For a genome bin whose contigs should stay together, use
[Classify a directory of bins](classify-bins.md).

## 1. Select a nucleotide FASTA

Run from the GVClass repository after [installation and database setup](../tutorials/getting-started.md).
Copy a bundled bin to try per-contig classification:

```bash
cp example/GVMAG-S-1096109-37.fna assembly.fna
```

For your own run, replace `assembly.fna` with your assembly file. Each FASTA
record must have a distinct identifier. GVClass converts filename-unsafe
characters in identifiers to underscores and stops if two retained contigs
would receive the same name.

## 2. Classify the contigs

```bash
pixi run gvclass assembly.fna --contigs -o contig_results -t 8
```

The default minimum is 10,000 bp per contig. Shorter sequences are excluded.
To change it, add `--contigs-min-length 30000` for a 30,000 bp minimum, or
`--contigs-min-length 0` to keep all contigs. This setting also controls length
validation of the split files.

A directory of nucleotide FASTA files can replace `assembly.fna`. GVClass
prefixes each contig's name with its source filename in that case.

## 3. Read the per-contig results

```bash
cut -f1,2,4 contig_results/gvclass_summary.tsv
```

Each successfully processed contig has its own row. Completeness and
contamination estimates describe that contig; they do not describe the full
assembly. Contigs without sufficient marker evidence can remain unclassified.

GVClass removes temporary split files after the run. Results remain in
`contig_results/`. See [Assess genome quality](assess-genome-quality.md) for
interpretation and the [output reference](../reference/output.md) for all files.
