<p align="center">
  <img src="images/GVClass_logo.png" alt="GVClass logo" width="50%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-2.0.3-blue.svg" alt="Version">
  <img src="https://img.shields.io/badge/license-non--commercial-orange.svg" alt="License: non-commercial use only">
</p>

# GVClass

GVClass classifies giant-virus genomes and metagenome-assembled genomes using marker-gene trees. It supports Nucleocytoviricota (NCLDV), Mirusviricota and Preplasmiviricota (PPV), and reports taxonomy, marker counts, completeness and contamination estimates.

## Install and run

On Linux x86-64 with [Pixi](https://pixi.sh/latest/installation/) and Git installed:

```bash
git clone https://github.com/NeLLi-team/gvclass.git
cd gvclass
pixi install --frozen
pixi run setup-db
pixi run gvclass example -o example_results -t 8
head -n 4 example_results/gvclass_summary.tsv
```

Create a directory for your genomes:

```bash
mkdir -p query_genomes
```

Put one genome per `.fna` file or one proteome per `.faa` file in this directory, then run:

```bash
pixi run gvclass query_genomes -o results -t 8
```

To build a shared species tree as well:

```bash
pixi run gvclass query_genomes -o species_tree_results -t 8 --species-tree-combined
```

## Documentation

- [Getting started](https://NeLLi-team.github.io/gvclass/tutorials/getting-started/): installation, bundled examples and output.
- [Classify genome bins](https://NeLLi-team.github.io/gvclass/how-to/classify-bins/): input preparation and resume.
- [Build a species tree](https://NeLLi-team.github.io/gvclass/how-to/build-a-species-tree/): per-query and combined trees.
- [Run on a cluster](https://NeLLi-team.github.io/gvclass/how-to/run-on-hpc/): Slurm and Apptainer.
- [Command-line options](https://NeLLi-team.github.io/gvclass/reference/cli/) and [output columns](https://NeLLi-team.github.io/gvclass/reference/output/).
- [Methods](https://NeLLi-team.github.io/gvclass/explanation/how-it-works/): taxonomy and quality estimates.

## Citation

> Pitot et al. (2024): Conservative taxonomy and quality assessment of giant virus genomes with GVClass. npj Viruses. https://www.nature.com/articles/s44298-024-00069-7

## Database sources

The trained contamination model is shipped in the runtime resource bundle at
`resources/contamination/model.joblib`.

The v2.0.0 runtime resource bundle is archived on Zenodo:
https://doi.org/10.5281/zenodo.21225457

Reference data sources:

> Medvedeva S, Guyet U, Pelletier E, et al. (2026): Widespread and intron-rich mirusviruses are predicted to reproduce in nuclei of unicellular eukaryotes. Nature Microbiology 11:228-239. https://doi.org/10.1038/s41564-025-02190-6

> Roux S, Fischer MG, Hackl T, Katz LA, Schulz F, Yutin N (2023): Updated Virophage Taxonomy and Distinction from Polinton-like Viruses. Biomolecules 13(2):204. https://doi.org/10.3390/biom13020204

> Fiamenghi MB, Camargo AP, Chasapi IN, et al. (2025): Meta-virus resource (MetaVR): expanding the frontiers of viral diversity with 24 million uncultivated virus genomes. Nucleic Acids Research gkaf1283. https://doi.org/10.1093/nar/gkaf1283

> Vasquez YM, Nardi T, Terasaki GM, et al. (2025): Genome-resolved expansion of Nucleocytoviricota and Mirusviricota reveals new diversity, functional potential, and biotechnological applications. bioRxiv 2025.09.26.678796. https://doi.org/10.1101/2025.09.26.678796

> Bellas CM, Sommaruga R (2026): A framework for Polinton-like virus diversity across aquatic microbiomes reveals links to multiple viral classes and Nucleocytoviricota. bioRxiv 2026.06.19.733378. https://doi.org/10.64898/2026.06.19.733378

## License and contact

Non-commercial use only (see `LICENCE`). Report issues at https://github.com/NeLLi-team/gvclass/issues or contact fschulz@lbl.gov.
