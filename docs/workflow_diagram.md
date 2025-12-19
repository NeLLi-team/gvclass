# GVClass Workflow Diagram

## Pipeline Overview

```mermaid
graph TD
    A[Input Files<br/>.fna/.faa] --> B{File Type?}
    
    B -->|Nucleotide .fna| C[Gene Calling<br/>pyrodigal]
    B -->|Protein .faa| D[Reformat Headers]
    
    C --> E[Test Genetic Codes<br/>0,1,4,6,11,15,29]
    E --> F[Select Best Code<br/>+2% coding density]
    F --> D
    
    D --> G[HMM Search<br/>pyhmmer]
    G --> H[Apply Cutoffs<br/>E-value & coverage]
    
    H --> I[Extract Markers<br/>99 GVOGs]
    
    I --> J[BLAST Search<br/>pyswrd/Diamond]
    J --> K[Get References<br/>E-value 1e-5]
    
    K --> L[Merge Sequences<br/>Query + References]
    L --> M[Align<br/>pyfamsa]
    
    M --> N[Trim<br/>pytrimal 10% gaps]
    N --> O{Tree Method?}
    
    O -->|FastTree| P[VeryFastTree<br/>Default]
    O -->|IQ-TREE| Q[IQ-TREE<br/>LG4X model]
    
    P --> R[Get Nearest Neighbors]
    Q --> R
    
    R --> S[Assign Taxonomy<br/>Consensus voting]
    S --> T[Quality Metrics<br/>Completeness/Contamination]
    
    T --> U[Summary Report<br/>gvclass_summary.tsv & csv]
```

## Parallel Processing

```mermaid
graph LR
    subgraph "Per Sample Processing"
        A1[Sample 1] --> B1[Gene Call] --> C1[HMM Search]
        A2[Sample 2] --> B2[Gene Call] --> C2[HMM Search]
        A3[Sample N] --> B3[Gene Call] --> C3[HMM Search]
    end
    
    subgraph "Per Marker Processing"
        C1 --> D1[Marker 1] --> E1[BLAST] --> F1[Tree]
        C1 --> D2[Marker 2] --> E2[BLAST] --> F2[Tree]
        C1 --> D3[Marker N] --> E3[BLAST] --> F3[Tree]
    end
    
    F1 --> G[Combine Results]
    F2 --> G
    F3 --> G
```

## Key Decision Points

1. **File Type Detection**
   - `.fna` → Gene calling required
   - `.faa` → Direct to HMM search

2. **Genetic Code Selection**
   - Tests codes: 0, 1, 4, 6, 11, 15, 29
   - Selects code with >2% better coding density than code 0

3. **Marker Filtering**
   - Must pass E-value threshold
   - Must have ≥66% HMM coverage
   - Must have ≥3 sequences for alignment

4. **Tree Method Selection**
   - FastTree: Default, faster
   - IQ-TREE: More accurate, slower

5. **Taxonomy Assignment**
   - Requires consensus among all nearest neighbors
   - Falls back to higher taxonomic levels if no consensus

## Performance Bottlenecks

```mermaid
graph TD
    A[Gene Calling<br/>🐌 Serial] -->|Bottleneck| B[Can parallelize<br/>genetic codes]
    
    C[BLAST Search<br/>🐌 pyswrd slow] -->|Bottleneck| D[Replace with<br/>Diamond/MMseqs2]
    
    E[Tree Building<br/>🐌 CPU bound] -->|Bottleneck| F[GPU acceleration<br/>for IQ-TREE]
    
    G[File I/O<br/>🐌 Many small files] -->|Bottleneck| H[Batch operations<br/>Archive format]
```
