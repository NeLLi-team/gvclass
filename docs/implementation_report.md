# GVClass Pipeline Implementation Report

## Executive Summary

This report analyzes the old Snakemake pipeline implementation and compares it with the current Prefect implementation, documenting each step, thresholds, parameters, and identifying areas for performance improvement.

## 1. Pipeline Overview

The GVClass pipeline classifies giant virus contigs/MAGs through phylogenetic analysis based on giant virus orthologous groups (GVOGs). Both implementations follow the same core workflow:

1. **Input Validation & Reformatting**
2. **Gene Calling with Genetic Code Selection**
3. **Database Search (HMM & BLAST)**
4. **Alignment & Phylogeny**
5. **Taxonomy Assignment**
6. **Quality Assessment**

## 2. Detailed Step Analysis

### 2.1 Input Validation & Reformatting

#### Old Implementation (reformat.py)
- **Input**: `.fna` (nucleic acid) or `.faa` (protein) files
- **Processing**:
  - Reformats sequence headers to `<filename>|<proteinid>` format
  - Replaces colons with double underscores in IDs
  - Calculates basic statistics (GC%, length, contig count)
- **Output**: Reformatted sequences and stats file

#### Current Implementation
- Functionality preserved in Prefect implementation
- Same header formatting rules applied

### 2.2 Gene Calling with Genetic Code Selection (opgecall.py)

#### Old Implementation
- **Tool**: Pyrodigal (Python version of Prodigal)
- **Genetic Codes Tested**: [0, 1, 4, 6, 15, 11, 29]
  - Code 0: Metagenomic mode
  - Code 1: Standard
  - Code 4: Mycoplasma/Spiroplasma
  - Code 6: Ciliate
  - Code 11: Bacterial/Archaeal/Plant Plastid
  - Code 15: Blepharisma
  - Code 29: Mold/Protozoan/Coelenterate Mitochondrial
- **Selection Criteria**:
  1. Run HMM search with each genetic code
  2. Calculate coding density percentage
  3. Select code with highest complete best hits count
  4. If tied, select code with highest average best hit score
  5. If still tied, select code with highest coding density (>2% improvement)
- **Thresholds**:
  - Completeness cutoff: 0.66 (66% HMM coverage)
  - Coding density improvement: 1.02x (2% better than metagenomic)

#### Performance Notes
- Runs pyhmmer search 7 times (once per genetic code)
- Could be parallelized across genetic codes

### 2.3 HMM Search (hmmsearch.py)

#### Old Implementation
- **Tool**: pyhmmer (Python HMMER implementation)
- **Parameters**:
  - E-value threshold: 10.0
  - Domain E-value threshold: 10.0
  - Threads: 4 (default)
- **Cutoffs**: Loaded from `models_APRIL24--databaseApril24.cutoffs`
  - Model-specific score and length cutoffs
- **Processing**:
  - Filters hits by score and length cutoffs
  - For general models: keeps best hit per protein
  - For order-level models (OG*): keeps best hit per protein-model pair
- **Output**: Count and score matrices for all models

#### Current Implementation
- Same pyhmmer parameters and cutoff logic
- Parallelized across models in Prefect

### 2.4 BLAST Search & Reference Selection (blastp_reduce_merge.py)

#### Old Implementation
- **Tool**: Originally Diamond BLASTP, now pyswrd
- **Parameters**:
  - E-value threshold: 1e-5
  - Top hits: Up to 100 references per query
  - Scoring matrix: BLOSUM62
- **Processing**:
  1. BLAST query hits against reference database
  2. Select up to 100/n_queries references per query
  3. Merge query and reference sequences
- **Note**: Skip hits with 100% identity (removed in current code)

#### Performance Issues
- pyswrd implementation lacks detailed alignment info
- Falls back to simplified output format
- Could benefit from using Diamond or MMseqs2

### 2.5 Alignment & Trimming (align_trim.py)

#### Old Implementation
- **Alignment Tool**: pyfamsa (Python FAMSA implementation)
- **Alignment Options**:
  - `linsi`: More accurate (UPGMA guide tree)
  - `ginsi`: Global alignment (UPGMA guide tree)
  - `auto`: Default parameters
- **Trimming Tool**: pytrimal
- **Trimming Parameters**:
  - Gap threshold: 0.1 (keep columns with ≤10% gaps)
- **Requirements**: Minimum 3 sequences for alignment

#### Current Implementation
- Same tools and parameters
- Parallelized across markers

### 2.6 Tree Building (build_tree.py)

#### Old Implementation
- **Methods**:
  1. **FastTree** (default):
     - Tool: VeryFastTree (Python implementation)
     - Parameters: `-lg -spr 4 -mlacc 2 -slownni`
     - Falls back to original FastTree if VeryFastTree fails
  2. **IQ-TREE**:
     - Model: LG4X
     - Options: `-fast` (fast mode)
     - Threads: 4 (default)
- **Cleanup**: Removes all IQ-TREE temporary files except .treefile

#### Performance Notes
- Tree building is computationally intensive
- Could benefit from GPU acceleration (IQ-TREE supports CUDA)

### 2.7 Nearest Neighbor Analysis (get_nn_tree.py)

#### Old Implementation
- **Tool**: ETE3 for tree parsing
- **Algorithm**:
  1. For each query paralog in tree:
     - Find parent node
     - Collect all leaves under parent (excluding query)
     - Find closest relative by branch distance
  2. Extract taxonomic level (order) from closest relatives
  3. Generate consensus taxonomy based on frequency
- **Consensus Rules**:
  - Strict: 100% agreement required
  - Majority: >50% agreement required
  - Domain level: Hyphenated list if no consensus

### 2.8 Result Summarization (summarize.py)

#### Old Implementation
- **Marker Set Analysis**:
  - GVOG4M: 4 core markers (completeness indicator)
  - GVOG8M: 8 extended markers
  - MIRUS: 5 Mirusviricota-specific markers
  - MCP: 5 Major Capsid Protein variants
  - MRYA: 6 Metagenomic Russian Yokohama markers
  - UNI56: 56 universal single-copy markers
  - BUSCO: 255 eukaryotic markers
  - PHAGE: 20 phage markers
- **Quality Metrics**:
  - Completeness: Percentage of markers found
  - Duplication: Average copies per marker
  - Contamination: Based on cellular marker presence
- **Order-specific Completeness**:
  - Uses `order_completeness.tab` lookup table
  - Calculates completeness based on order-specific marker sets

## 3. Key Parameters and Thresholds

### 3.1 Configuration Parameters (config.yml)
- `database_path`: "resources" (default)
- `keep_temp`: True (retain intermediate files)
- `mafftoption`: "linsi" (alignment method)
- `treeoption`: "fasttree" (tree method)
- `mode_fast`: True (skip order-level marker trees)

### 3.2 Critical Thresholds
- **Minimum contig length**: 50kb recommended (20kb absolute minimum)
- **HMM coverage for complete hit**: 66%
- **Alignment minimum sequences**: 3
- **Trimming gap threshold**: 10%
- **BLAST E-value**: 1e-5
- **HMM E-value**: 10.0
- **Coding density improvement**: 2% over metagenomic mode

## 4. Comparison with Current Prefect Implementation

### 4.1 Preserved Functionality
- All core analysis steps maintained
- Same tools and parameters
- Identical output format
- Same quality metrics

### 4.2 Improvements in Prefect
- **Parallelization**: Tasks run concurrently with Dask
- **Error Handling**: Better error recovery and logging
- **Monitoring**: Real-time progress tracking
- **Modularity**: Cleaner separation of concerns

### 4.3 Potential Issues
- pyswrd implementation less feature-rich than Diamond
- Missing some optimization flags from original

## 5. Performance Improvement Recommendations

### 5.1 High-Impact Optimizations

1. **Parallel Genetic Code Testing**
   - Run all 7 genetic codes in parallel
   - Potential speedup: 3-5x for gene calling step

2. **Replace pyswrd with Diamond/MMseqs2**
   - Diamond is 20,000x faster than BLAST
   - MMseqs2 offers similar speed with better sensitivity
   - Provides detailed alignment statistics

3. **GPU Acceleration for Tree Building**
   - IQ-TREE supports CUDA acceleration
   - Potential speedup: 5-10x for large trees

4. **Caching Strategy**
   - Cache HMM search results by sequence hash
   - Cache tree files for identical alignments
   - Implement database versioning

### 5.2 Algorithmic Optimizations

1. **Early Termination in Gene Calling**
   - Stop testing genetic codes once clear winner emerges
   - Implement progressive testing (most likely codes first)

2. **Adaptive Reference Selection**
   - Use clustering to reduce redundant references
   - Implement taxonomy-aware reference selection

3. **Batch Processing**
   - Group similar-sized sequences for better load balancing
   - Implement dynamic batching based on sequence length

### 5.3 Infrastructure Optimizations

1. **Memory Management**
   - Stream large files instead of loading fully
   - Implement sequence chunking for very large contigs

2. **I/O Optimization**
   - Use binary formats where possible (e.g., HDF5 for matrices)
   - Implement parallel I/O for large datasets

3. **Distributed Computing**
   - Full Dask cluster support for HPC environments
   - Implement checkpointing for long-running analyses

## 6. Conclusion

The GVClass pipeline implementation is scientifically sound with well-chosen thresholds and parameters. The migration to Prefect has improved parallelization and error handling while maintaining scientific accuracy. The main opportunities for performance improvement lie in:

1. Better parallelization of genetic code selection
2. Replacing pyswrd with faster alternatives
3. Implementing intelligent caching
4. GPU acceleration for computationally intensive steps

These optimizations could potentially reduce runtime by 5-10x for typical datasets while maintaining or improving accuracy.