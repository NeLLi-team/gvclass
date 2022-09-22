n environmental microbiomes. GVClass assigns taxonomy to putative giant virus contigs or metagenome assembled genomes ([GVMAGs](https://doi.org/10.1038/s41586-020-1957-x)). It uses a conservative approach based on the consensus of single protein trees built from up to 9 giant virus orthologous groups ([GVOGs](https://doi.org/10.1371/journal.pbio.3001430)). These 9 GVOGs are often conserved across different lineages in the viral phylum Nucleocytoviricota.
## Running GVClass

### Requirements
* Input is a directory that contains single contigs or MAGs as nucleic acid (fna) or proteins (faa)
* File extensions .fna or .faa
* Recommended length for single contigs is 50kb, but at least 20kb
* No special characters (".", ";", ":") in filebase name, "\_" or "-" are okay
* No whitespace in sequence headers
* Recommended sequence header format if faa provided: <filenamebase>|<proteinid>
* Input will be checked and reformatted

### IMG/VR
* Upload you metagenome assembled genome or single contig to [IMG/VR](https://img.jgi.doe.gov/vr/) using the GVClass feature

### Docker / Shifter
* The recommended way to run GVClass (especially on Windows or Mac) is to use the Docker image [fschulzjgi/gvclass](https://hub.docker.com/repository/docker/fschulzjgi/gvclass)

```
docker pull fschulzjgi/gvclass:0.9.2
```
* Query genomes to test GVClass are in this example in a directory called "test"
```
docker run -t -i -v $(pwd)/test:/gvclass/querydir --user $(id -u):$(id -g) fschulzjgi/gvclass:0.9.2 snakemake --use-conda -j 4 --config querydir="querydir"
```
* If this completes successfully, run it using your own directory of query genomes
```
docker run -t -i -v <path to query dir>:/gvclass/querydir --user $(id -u):$(id -g) fschulzjgi/gvclass:0.9.2 snakemake --use-conda -j 64 --config querydir="querydir"
```
* An alternative way to run GVClass is to use Shifter instead of Docker
```
shifterimg pull fschulzjgi/gvclass:0.9.2
shifter \
  --volume=<path to query dir>:/gvclass/example \
  --image=fschulzjgi/gvclass:0.9.2 \
  bash -c \
  "snakemake \
  --snakefile /gvclass/workflow/Snakefile \
  -j 64 \
  --use-conda \
  --config querydir="/gvclass/example" \
  --conda-prefix /gvclass/.snakemake/conda "
```

### Snakemake workflow
* GVClass can also be run using the snakemake workflow directly
```
git clone https://github.com/NeLLi-team/gvclass
```
* Test GVClass using the provided giant virus assemblies
```
snakemake -j 64 --use-conda --config querydir="example"
```
* If this completes successfully, run it using your own directory of query genomes
```
snakemake -j <number of processes> --use-conda --config querydir="<path to query dir>"
```

## Interpretation of the results
* The classification result is summarized in a tab separated file \<query name\>.summary.tab

### Taxonomy assignments
*  Taxonomy assignments are provided on different taxonomic levels. 
*  To yield an assignments all nearest neighbors in GVOG phylogenetic tree have to be in agreement. 
*  Depending on the number of identified unique GVOGs an assignment is provided if at least 1 GVOG (stringency "gte1"), 2 GVOGs (stringency "gte2") or 3 GVOGs (stringency "gte3") were found.
*  Less stringent is the "majority" approach, here more than 50% of identified markers have to yield the same taxonomy to enable and assignment. 
*  If taxonomy assignments are not in agreement at a low taxonomy level (e.g. species, genus, or family) then the next higher taxonomy level will be evaluated, up to the domain level.

### Contamination
* Giant virus genomes typically have less than 10 out of a set of 56 universal cellular housekeeping genes (UNI56). Higher UNI56 counts indicate cellular contamination, or giant virus sequences that are located on host contigs.
  * UNI56u (unique counts), UNI56t(total counts), UNI56df (duplication factor) are provided and can be used for further quality filtering
* Giant virus genomes typically have a duplication factor of GVOG7 (as subset of GVOG9) of below 3. Higher GVOG7 duplication factors indicate the presence mixed viral populations.
  * GVOG7u (unique counts), GVOG7t(total counts), GVOG7df (duplication factor) are provided and can be used for further quality filtering 
* Giant viruses may break any of these rules and novel giant virus genomes are often full of surprises. Thus, GVClass does not perform automatic quality filtering based on marker gene counts.

## Benchmarking

## References
1. [Schulz F, Roux S, Paez-Espino D, Jungbluth S, Walsh DA, Denef VJ, McMahon KD, Konstantinidis KT, Eloe-Fadrosh EA, Kyrpides NC, Woyke T. Giant virus diversity and host interactions through global metagenomics. Nature. 2020 Feb;578(7795):432-6.](https://doi.org/10.1038/s41586-020-1957-x)
2. [Aylward FO, Moniruzzaman M, Ha AD, Koonin EV. A phylogenomic framework for charting the diversity and evolution of giant viruses. PLoS biology. 2021 Oct 27;19(10):e3001430.](https://doi.org/10.1371/journal.pbio.3001430)

## Acknowledgements
GVClass was developed by the [New Lineages of Life Group](https://jgi.doe.gov/our-science/scientists-jgi/new-lineages-of-life/) at the DOE Joint Genome Institute supported by the Office of Science of the U.S. Department of Energy under contract no. DE-AC02-05CH11231.


