# GVClass

GVClass classifies giant-virus genomes and metagenome-assembled genomes using marker-gene trees. It supports Nucleocytoviricota (NCLDV), Mirusviricota and Preplasmiviricota (PPV), and reports taxonomy, marker counts, completeness and contamination estimates.

## Install and run an example

On Linux x86-64 with [Pixi](https://pixi.sh/latest/installation/) and Git installed:

```bash
git clone https://github.com/NeLLi-team/gvclass.git
cd gvclass
pixi install --frozen
pixi run setup-db
pixi run gvclass example -o example_results -t 8
head -n 4 example_results/gvclass_summary.tsv
```

The example contains two nucleotide assemblies and one protein set. Follow [Getting started](tutorials/getting-started.md) for input preparation and output interpretation. For clusters, use the [Apptainer and Slurm instructions](how-to/run-on-hpc.md).

## Run your genomes

Create an input directory in the repository:

```bash
mkdir -p query_genomes
```

Place one genome per `.fna` file or one proteome per `.faa` file in this directory, then run:

```bash
pixi run gvclass query_genomes -o results -t 8
```

To also build a shared species tree:

```bash
pixi run gvclass query_genomes -o species_tree_results -t 8 --species-tree-combined
```

See [Build a species tree](how-to/build-a-species-tree.md) for input requirements and tree files.

## Documentation

| Task | Page |
| --- | --- |
| Install and run the bundled examples | [Getting started](tutorials/getting-started.md) |
| Classify genome bins | [Classify a directory of bins](how-to/classify-bins.md) |
| Classify contigs separately | [Classify individual contigs](how-to/classify-contigs.md) |
| Run a batch job | [HPC instructions](how-to/run-on-hpc.md) |
| Read the results | [Output files and columns](reference/output.md) |
| Look up an option | [Command-line reference](reference/cli.md) |
| Understand the method | [How GVClass works](explanation/how-it-works.md) |

## Citation and license

Cite Pitot et al. (2024), [*Conservative taxonomy and quality assessment of giant virus genomes with GVClass*](https://www.nature.com/articles/s44298-024-00069-7), npj Viruses.

GVClass is licensed for non-commercial use only. See [LICENCE](https://github.com/NeLLi-team/gvclass/blob/main/LICENCE). Report problems through [GitHub Issues](https://github.com/NeLLi-team/gvclass/issues).
