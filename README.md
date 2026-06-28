<p align="center">
  <img src="images/GVClass_logo.png" alt="GVClass logo" width="50%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-2.0.0-blue.svg" alt="Version">
  <img src="https://img.shields.io/badge/license-BSD--3--Clause-green.svg" alt="License">
  <img src="https://img.shields.io/badge/pixi-enabled-orange.svg" alt="Pixi">
</p>

# GVClass — Giant Virus Classification

GVClass assigns taxonomy to giant virus genomes and metagenome-assembled genomes by phylogenetic placement against a comprehensive reference genome database. It covers Nucleocytoviricota (NCLDV), Mirusviricota, and Preplasmiviricota (PPV). It predicts taxonomy from domain down to genus and species level, assigning a nearest-reference label. Each genome returns a majority-vote taxonomy with a confidence flag, plus completeness and contamination metrics tuned for giant viruses.

## Documentation

**Full documentation: https://NeLLi-team.github.io/gvclass/**

- [Tutorials](https://NeLLi-team.github.io/gvclass/tutorials/) — learn GVClass on the bundled example.
- [How-to guides](https://NeLLi-team.github.io/gvclass/how-to/) — bins, contigs, HPC, species trees, speed/accuracy, quality.
- [Reference](https://NeLLi-team.github.io/gvclass/reference/) — every CLI flag, config key, output column, and marker panel.
- [Explanation](https://NeLLi-team.github.io/gvclass/explanation/) — how placement, taxonomy, and the quality models work.

## Quick start

Pixi (local):

```bash
git clone -b gvclass-dev https://github.com/NeLLi-team/gvclass.git
cd gvclass
pixi install
pixi run setup-db
pixi run example
```

Apptainer (HPC):

```bash
wget https://raw.githubusercontent.com/NeLLi-team/gvclass/gvclass-dev/gvclass-a
chmod +x gvclass-a
./gvclass-a QUERY_DIR RESULTS_DIR -t 32
```

The Apptainer image bundles the database and dependencies. See [Getting started](https://NeLLi-team.github.io/gvclass/tutorials/getting-started/) for the full walkthrough.

## Input

GVClass works best on bins after metagenomic binning: a directory of one or more FASTA files (`.fna` or `.faa`), one file per putative genome. For giant virus discovery, filter contigs to >=30 kb (>=50 kb preferred). Use `--contigs` to classify each contig in a multi-contig `.fna` independently. Details are in [the how-to guides](https://NeLLi-team.github.io/gvclass/how-to/).

## Citation

> Pitot et al. (2024): Conservative taxonomy and quality assessment of giant virus genomes with GVClass. npj Viruses. https://www.nature.com/articles/s44298-024-00069-7

## Database sources

The GVClass runtime resources include genomes and models derived from:

> Medvedeva S, Guyet U, Pelletier E, et al. (2026): Widespread and intron-rich mirusviruses are predicted to reproduce in nuclei of unicellular eukaryotes. Nature Microbiology 11:228-239. https://doi.org/10.1038/s41564-025-01906-2

> Roux S, Fischer MG, Hackl T, Katz LA, Schulz F, Yutin N (2023): Updated Virophage Taxonomy and Distinction from Polinton-like Viruses. Biomolecules 13(2):204. https://doi.org/10.3390/biom13020204

> Fiamenghi MB, Camargo AP, Chasapi IN, et al. (2025): Meta-virus resource (MetaVR): expanding the frontiers of viral diversity with 24 million uncultivated virus genomes. Nucleic Acids Research gkaf1283. https://doi.org/10.1093/nar/gkaf1283

> Vasquez YM, Nardi T, Terasaki GM, et al. (2025): Genome-resolved expansion of Nucleocytoviricota and Mirusviricota reveals new diversity, functional potential, and biotechnological applications. bioRxiv 2025.09.26.678796. https://doi.org/10.1101/2025.09.26.678796

> Bellas CM, Sommaruga R (2026): A framework for Polinton-like virus diversity across aquatic microbiomes reveals links to multiple viral classes and Nucleocytoviricota. bioRxiv 2026.06.19.733378. https://doi.org/10.64898/2026.06.19.733378

## License and contact

BSD 3-Clause (see `LICENCE`). Report issues at https://github.com/NeLLi-team/gvclass/issues or contact fschulz@lbl.gov.
