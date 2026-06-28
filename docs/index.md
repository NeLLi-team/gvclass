# GVClass

GVClass assigns taxonomy to giant virus genomes and metagenome-assembled genomes by phylogenetic placement against a comprehensive reference genome database. It covers Nucleocytoviricota (NCLDV), Mirusviricota, and Preplasmiviricota (PPV). It predicts taxonomy from domain down to genus and species level, assigning a [nearest-reference label](explanation/taxonomy.md).

Each genome returns a majority-vote taxonomy with a [confidence flag](reference/output.md), plus [completeness and contamination metrics](explanation/quality-metrics.md) tuned for giant viruses. The taxonomy comes from a per-marker tree nearest-neighbor vote across the GVOG panels.

## Where to go

- [Tutorials](tutorials/index.md): learn GVClass by running it on the bundled example.
- [How-to guides](how-to/index.md): task recipes for bins, contigs, HPC, and tuning.
- [Reference](reference/index.md): every flag, the 44 output columns, and config keys.
- [Explanation](explanation/index.md): how placement, taxonomy, and the quality models work.

## Quick install

=== "Pixi (local)"

    ```bash
    git clone -b gvclass-dev https://github.com/NeLLi-team/gvclass.git
    cd gvclass
    pixi install
    pixi run setup-db
    pixi run example
    ```

=== "Apptainer (HPC)"

    ```bash
    wget https://raw.githubusercontent.com/NeLLi-team/gvclass/gvclass-dev/gvclass-a
    chmod +x gvclass-a
    ./gvclass-a QUERY_DIR RESULTS_DIR -t 32
    ```

For the full walkthrough, see [Getting started](tutorials/getting-started.md).

!!! note "Versions on this branch"

    Software 2.0.0, database bundle v1.7.1. Downloads and changelogs are on the [GitHub Releases](https://github.com/NeLLi-team/gvclass/releases) page.

## Citation, license, contact

Cite Pitot et al. (2024), [*Conservative taxonomy and quality assessment of giant virus genomes with GVClass*](https://www.nature.com/articles/s44298-024-00069-7), npj Viruses.

Licensed under BSD-3-Clause. Questions and bugs go to [GitHub Issues](https://github.com/NeLLi-team/gvclass/issues) or fschulz@lbl.gov.
