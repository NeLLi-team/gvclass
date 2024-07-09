```
                                                       ``,:'-                     
                                                    `-'~^*r?^^_'`                 
                                                -',~";*>>*r|r?rr:,-              
                                            `.',,"_:>\>\?>|r\r|r?^>;:'`          
                                         ``',~~_:;;>LL?/||\r|*?*?>r*rr^r;,-       
                                      -':;;^>*\cziy9Uhqqh66PXeauujzf7vLL\\|?;:-`  
                                   `;>\LLL//c77z1f6N6UH9KRPR6URPRZHw966kPP6PPeet' 
                                   `*\Li///cvvzJ76Ne&6ZdRw8UqWw&PdU9K6OPOhqhqUKz` 
                                   `>\\\\L77vz1vHR6#ODWwXgUK&PH6OK9d6H6H6OPOPKei  
                                ______ _    _ _______        _______ _______ _______ 
                               |  ____  \  /  |       |      |_____| |______ |______
                               |_____|   \/   |_____  |_____ |     | ______| ______| 
                                    ^zzztzza9RdH8KgdHDqNOdN99R9DqHHUDUH6H6OPqOu`  
                                    ;tttJzjq89&DON9RRdNU8W6NdwOKK&hR689H9d6KOq/   
                                    ;Jttzj68UNdOg9NHOg6BOKQ6H#mKgXghDKqO699RKU^   
                                     :zyUR9R6gDHN9BRDBU#qd8HRDROWqhqUW9DqkgOgq?    
                                       ,S9D6H9gRDNKQKRN6gKRHDODDqRhd89OOORhDRar-     
                                        `:7{N9B9NNKQqNgqNOHHRK&OORkgOON96RUqS_        
                                          `:7U@@@@@@@@@@@@Q@QQQBQgQNW#P\:-         `
                                             -;LOB@@@@@@@@@@@@@@@QDJ:`             
                                               `:vwQ@@@@@Q@@Bk^'                 
                                                  `;1hQQKv:`              `     
```

# gvclass
_version 1.0 8 July 2024_

Giant viruses are abundant and diverse and frequently found in environmental microbiomes. gvclass assigns taxonomy to putative giant virus contigs or metagenome assembled genomes ([GVMAGs](https://doi.org/10.1038/s41586-020-1957-x)). It uses a conservative approach based on the consensus of single protein trees built from giant virus orthologous groups ([GVOGs](https://doi.org/10.1371/journal.pbio.3001430)), additional Mirusvirus, Mryavirus and Poxvirus hallmark genes and cellular single copy panorthologs. Further, a gene content based classifier predicts giant virus order-level affiliation and genome completeness and contamination is then estimated based on copy numbers of a larger set of genes typically conserved at order-level.

## Running gvclass

### Requirements
* Conda environment wih snakemake, check here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
* Input is a directory that contains single contigs or MAGs as nucleic acid (fna) or proteins (faa)
* File extensions .fna or .faa
* Recommended length for single contigs is 50kb, but at least 20kb
* No special characters (".", ";", ":") in filebase name, "\_" or "-" are okay
* No whitespace in sequence headers
* Recommended sequence header format if faa provided: <filenamebase>|<proteinid>
* Input will be checked and reformatted if necessary
* Domain- and order-level classifier will only be used on genomes that have at least 3 features


### Settings
* Config file allows to specify options for MAFFT (default is mafft-linsi), iqtree (default) or fasttree
* fast_mode (default) can be set to False in config file, in that case single protein trees are also built for all conserved order-level marker genes

### IMG/VR
* Upload you metagenome assembled genome or single contig to [IMG/VR](https://img.jgi.doe.gov/vr/) using the gvclass feature

### Snakemake workflow
* Clone the repository
```
git clone --recurse-submodules https://github.com/NeLLi-team/gvclass
```
* Test gvclass using the provided giant virus assemblies
```
snakemake -j 24 --use-conda --config querydir="example"
```
* If this completes successfully, run it using your own directory of query genomes
```
snakemake -j <number of processes> --use-conda --config querydir="<path to query dir>"
```

### Docker / Apptainer

* Get the sif from
```
wget https://portal.nersc.gov/cfs/nelli/gvclass.sif
```

* Run it like this

```
bash gvclass_apptainer.sh <querydir>
```

* Or 

## Interpretation of the results
* The classification result is summarized in a tab separated file \<query name\>.summary.tab

### Gene calling
* Different genetic codes are tested and evaluated based on hmmsearch using the general models
* Genetic code that yields the largest number of matches to general models with the highest average bitscore and the highest coding density is selected

### Taxonomy assignments
* Taxonomy assignments are provided on different taxonomic levels
* To yield an assignments all nearest neighbors in GVOG phylogenetic tree have to be in agreement

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


