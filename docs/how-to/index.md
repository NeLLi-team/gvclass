# How-to guides

Each guide here solves a single goal and assumes you already have the basics from [Getting started](../tutorials/getting-started.md).

- [Classify a directory of bins](classify-bins.md): the recommended workflow, one bin per FASTA file.
- [Run on an HPC cluster](run-on-hpc.md): run with the Apptainer wrapper, no local install needed.
- [Classify contigs separately](classify-contigs.md): treat each contig of a metagenome as its own genome.
- [Build a species tree](build-a-species-tree.md): a genome-level supermatrix tree and nearest-reference placement.
- [Tune speed and accuracy](tune-speed-and-accuracy.md): fast mode, tree method, and thread layout.
- [Configure the database](configure-the-database.md): set the database location, share it, and keep it current.
- [Assess genome quality](assess-genome-quality.md): decide whether a GVMAG is high quality.

!!! tip

    For exact flag definitions see the [CLI reference](../reference/cli.md) and the [output columns](../reference/output.md). For the reasoning behind each step, follow into [Explanation](../explanation/index.md).
