![screenshot](GVClass_logo.png)


# gvclass
_version 1.0 8 July 2024_

Giant viruses are abundant and diverse and frequently found in environmental microbiomes. gvclass assigns taxonomy to putative giant virus contigs or metagenome assembled genomes ([GVMAGs](https://doi.org/10.1038/s41586-020-1957-x)). It uses a conservative approach based on the consensus of single protein trees built from giant virus orthologous groups ([GVOGs](https://doi.org/10.1371/journal.pbio.3001430)), additional Mirusvirus, Mryavirus and Poxvirus hallmark genes and cellular single copy panorthologs. Genome completeness and contamination is then estimated based on copy numbers of a larger set of genes typically conserved in single copy at order-level.

## Running gvclass

### Requirements
* Conda environment wih snakemake, check here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
* Input is a directory that contains single contigs or MAGs as nucleic acid (fna) or proteins (faa)
* File extensions .fna or .faa
* Recommended length for assembly size is 50kb, but at least 20kb
* No special characters (".", ";", ":") in filebase name, "\_" or "-" are okay
* Recommended sequence header format if faa provided: <filenamebase>|<proteinid>
* Input will be checked and reformatted if necessary

### Settings
* Config file allows to specify options for MAFFT (default is mafft-linsi), iqtree (default) or fasttree
* fast_mode (default) can be set to False in config file, in that case single protein trees are also built for all conserved order-level marker genes

### IMG/VR
* Upload you metagenome assembled genome or single contig to [IMG/VR](https://img.jgi.doe.gov/vr/) using the gvclass feature

### Snakemake workflow

* The preferred way to run it is with apptainer, for snakemake our modified version of pyrodigal has to be set up first

### Docker / Apptainer

* Run it like this with apptainer using the providing gvclass_apptainer.sh

```
bash gvclass_apptainer.sh <querydir> <n processes>
```

* Or directly

```
PROCESSES=<number of processes, e.g. 8>
QUERYDIR=<dir with query genomes, e.g. example>

apptainer run docker://docker.io/doejgi/gvclass:latest \
  snakemake --snakefile /gvclass/workflow/Snakefile \
           -j $PROCESSES \
           --use-conda \
           --conda-frontend mamba \
           --conda-prefix /gvclass/.snakemake/conda \
           --config querydir="$QUERYDIR"
           database_path="/gvclass/resources"
```

* Alternatively pull the image and use docker or shifter

```
docker pull doejgi/gvclass:latest
```
  
* Clone the repository
```
git clone --recurse-submodules https://github.com/NeLLi-team/gvclass
```

* Activate snakemake (8.14.0) conda environment, install cython and pyrodigal
```
pip install cython
cd workflow/scripts/
pip install --user ./pyrodigal
```
* Test gvclass using the provided giant virus assemblies
```
snakemake -j 24 --use-conda --config querydir="example"
```
* If this completes successfully, run it using your own directory of query genomes
```
snakemake -j <number of processes> --use-conda --config querydir="<path to query dir>"
```

## Interpretation of the results
* The classification result is summarized in a tab separated file in a subdir "results" in the the query dir

### Gene calling
* Different genetic codes are tested and evaluated based on hmmsearch using the general models
* Genetic code that yields the largest number of matches to general models with the highest average bitscore and the highest coding density is selected

### Taxonomy assignments
* Taxonomy assignments are provided on different taxonomic levels
* To yield an assignments all nearest neighbors in GVOG phylogenetic trees have to be in agreement

### Contamination
* Giant virus genomes typically have less than 10 out of a set of 56 universal cellular housekeeping genes (UNI56). Higher UNI56 counts indicate cellular contamination, or giant virus sequences that are located on host contigs.
  * UNI56u (unique counts), UNI56t(total counts), UNI56df (duplication factor) are provided and can be used for further quality filtering
* Giant virus genomes typically have a duplication factor of GVOG7 and  GVOG9 of below 3. Higher GVOG7 duplication factors indicate the presence mixed viral populations.
  * GVOG8u, GVOG4u (unique counts), GVOG8t, GVOG4t (total counts), GVOG8df (duplication factor) are provided and can be used for further quality filtering
     * GVOG8df < 2 and order_dup < 1.5: low chance of representing mixed bin [high quality]
     * GVOG8df 2-3 and order_dup 1.5-2: medium chance of representing mixed bin [medium quality]
     * GVOG8df >3 and order_dup >3: high chance of representing mixed bin [low quality]
### Completeness
* Genome completeness estimate based on count of genes conserved in 50% of genomes of the respective Nucleocytoviricota order. 
  * \< 30%: low completeness  [low quality]
  * 30-70%: medium completeness [medium quality]
  * \> 70% high completeness [high quality]

## Benchmarking

* Will be provided soon

## Citation

* Will be provided soon

## References
1. [Schulz F, Roux S, Paez-Espino D, Jungbluth S, Walsh DA, Denef VJ, McMahon KD, Konstantinidis KT, Eloe-Fadrosh EA, Kyrpides NC, Woyke T. Giant virus diversity and host interactions through global metagenomics. Nature. 2020 Feb;578(7795):432-6.](https://doi.org/10.1038/s41586-020-1957-x)
2. [Aylward FO, Moniruzzaman M, Ha AD, Koonin EV. A phylogenomic framework for charting the diversity and evolution of giant viruses. PLoS biology. 2021 Oct 27;19(10):e3001430.](https://doi.org/10.1371/journal.pbio.3001430)

## Acknowledgements
GVClass was developed by the [New Lineages of Life Group](https://jgi.doe.gov/our-science/scientists-jgi/new-lineages-of-life/) at the DOE Joint Genome Institute supported by the Office of Science of the U.S. Department of Energy under contract no. DE-AC02-05CH11231.


