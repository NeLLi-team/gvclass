<p align="center">
  <img src="images/GVClass_logo.png" alt="GVClass Logo" width="50%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-v1.2.2-blue.svg" alt="Version">
  <img src="https://img.shields.io/badge/license-BSD--3--Clause-green.svg" alt="License">
  <img src="https://img.shields.io/badge/python-3.11-blue.svg" alt="Python">
  <img src="https://img.shields.io/badge/pixi-enabled-orange.svg" alt="Pixi">
</p>

# GVClass - Giant Virus Classification Tool

GVClass assigns taxonomy to giant virus contigs and metagenome-assembled genomes (GVMAGs). It uses phylogenetic analysis based on giant virus orthologous groups (GVOGs) to provide accurate classification from domain to species level.

## Quick Start

### Option 1: Pixi (Local/Development)

**Best for:** Contributing code, modifying the pipeline, or running locally

```bash
# 1. Install Pixi (one-time)
curl -fsSL https://pixi.sh/install.sh | bash

# 2. Clone repository
git clone https://github.com/NeLLi-team/gvclass.git
cd gvclass

# 3. Install dependencies (pixi handles everything)
pixi install

# 4. Run GVClass (must run from repo directory)
pixi run gvclass <input_directory> -t 16

# Optional: install CLI wrappers into ~/bin
pixi run install-cli

# Test installation with example data
pixi run run-example
```

### Option 2: Apptainer/Singularity (HPC)

**Best for:** Running on HPC clusters, no installation needed

```bash
# Download the wrapper script
wget https://raw.githubusercontent.com/NeLLi-team/gvclass/main/gvclass-a
chmod +x gvclass-a

# Run from anywhere (image auto-downloads on first use)
./gvclass-a /path/to/query_genomes /path/to/results -t 32

# With options
./gvclass-a my_data my_results -t 32 --tree-method iqtree --mode-fast
```

`gvclass-a` auto-pulls `library://` images through the public Sylabs endpoint (`https://library.sylabs.io`), so users do **not** need an access token for normal pulls/runs.
The Apptainer image includes the database (~700MB) and all dependencies. No setup needed!

## Input Requirements

GVClass works best **after metagenomic binning**.

- **Recommended input**: A directory containing one or more bin files (`.fna` or `.faa`)
- **Bin definition**: One FASTA file with one or more contigs representing the same putative genome
- **Alternative mode**: Use `--contigs` to treat contigs in a multi-contig `.fna` as independent viral genomes
- **Length guidance for giant virus discovery**:
  - Minimum supported: ~20 kb
  - Better reliability: filter contigs to `>=30 kb`
  - Preferred when possible: `>=50 kb`
- **Clean filenames**: Avoid special characters (`.` `;` `:`), use `_` or `-` instead
- **Protein headers**: Format as `filename|proteinid` for best results

## Example Usage

### Using Pixi (from repo directory)

```bash
# Recommended: run on bins after metagenomic binning
pixi run gvclass my_genomes -o my_results -t 32

# With options
pixi run gvclass my_genomes -t 32 --mode-fast --tree-method iqtree -j 4

# Sensitive HMM search mode (uses E-value 1e-5 instead of GA cutoffs)
pixi run gvclass my_genomes -t 32 --sensitive

# Alternative: classify each contig in a multi-contig FNA separately
pixi run gvclass --contigs my_genome.fna -o results -t 32
```

### Using Apptainer (gvclass-a)

```bash
# Recommended: run on bins after metagenomic binning
./gvclass-a my_genomes my_results -t 32

# Fast mode (skip order-level markers for 2-3x speedup)
./gvclass-a my_genomes my_results -t 32 --mode-fast

# Use IQ-TREE for more accurate phylogeny (slower)
./gvclass-a my_genomes my_results -t 32 --tree-method iqtree

# Sensitive HMM search mode
./gvclass-a my_genomes my_results -t 32 --sensitive

# Control parallelization (4 workers × 8 threads = 32 total)
./gvclass-a my_genomes my_results -t 32 -j 4

# Alternative: classify each contig in a multi-contig FNA separately
./gvclass-a --contigs my_genome.fna -o results -t 32
```
## Output

Results are saved to `<input_name>_results/` containing:
- `gvclass_summary.csv`, `gvclass_summary.tsv`
- Individual query subdirectories with detailed analysis

### Output Columns Explained

| Column | Description |
|--------|-------------|
| query | Input filename |
| taxonomy_majority | Full taxonomy based on majority rule |
| taxonomy_strict | Conservative taxonomy (100% agreement) |
| species → domain | Individual taxonomic levels with taxon counts |
| avgdist | Average tree distance to references |
| order_dup | Duplication factor indicating contamination level |
| order_completeness | Order-specific completeness (% unique markers found) |
| gvog4_unique | Count of unique GVOG4 markers found |
| gvog8_unique/total/dup | GVOG8 marker counts and duplication |
| ncldv_mcp_total | NCLDV-specific MCP marker count |
| mcp_total | All MCP marker count (NCLDV + Mirus) |
| vp_completeness | Virophage completeness (n/4 core markers: MCP, Penton, ATPase, Protease) |
| vp_mcp | Count of proteins with VP MCP marker hits |
| plv | Count of proteins with PLV marker hits (single PLV marker; values can be 0..N) |
| vp_df | Virophage duplication factor (total VP hits / 4) |
| mirus_completeness | Mirusviricota completeness (n/4 core markers: MCP, ATPase, Portal, Triplex) |
| mirus_df | Mirusviricota duplication factor |
| mrya_unique/total | Mryavirus-specific marker counts |
| phage_unique/total | Phage marker counts |
| cellular_unique/total/dup | Cellular contamination markers |
| contigs | Number of contigs |
| LENbp | Total length in base pairs |
| GCperc | GC content percentage |
| genecount | Number of predicted genes |
| CODINGperc | Coding density percentage |
| ttable | Genetic code used |
| weighted_order_completeness | Weighted completeness score considering marker importance |

## Configuration (Optional)

Create `gvclass_config.yaml` to set defaults:

```yaml
database:
  path: resources                    # Relative path: <gvclass_repo>/resources
  # path: /media/shared-expansion/dbs/gvclass_resources  # Absolute path on shared storage

pipeline:
  tree_method: fasttree             # or 'iqtree' for more accuracy
  mode_fast: false                  # Skip order-level marker trees when true (speeds up analysis)
  sensitive_mode: false             # Use E-value 1e-5 for pyhmmer instead of GA cutoffs
  contigs_min_length: 10000         # In --contigs mode, skip contigs shorter than this (bp)
  threads: 16                       # Default thread count
```

Database path precedence:
- `--database` CLI flag
- `GVCLASS_DB` environment variable
- `database.path` in config
- default `resources` in the GVClass repo

Examples:
```bash
# Use shared database location for all runs
export GVCLASS_DB=/media/shared-expansion/dbs/gvclass_resources
pixi run setup-db
pixi run gvclass example -o example_results
```

## What's New in v1.2.2

- **HPC Runtime Robustness**: Removed runtime dependence on Prefect task orchestration to avoid ephemeral Prefect API startup failures on restricted systems.
- **Container Release**: Updated public Apptainer image tag to `1.2.2`.
- **Database Version**: Unchanged (`v1.2.1`, tarball `resources_v1_2_1.tar.gz`).

## Advanced Usage

### Advanced Container Usage

The `gvclass-a` wrapper handles container execution automatically. For manual control:

```bash
# Pull the image manually (works without auth token for public images)
apptainer pull --library https://library.sylabs.io \
  gvclass_1.2.2.sif library://nelligroup-jgi/gvclass/gvclass:1.2.2

# Run with manual bind mounts
apptainer run -B /path/to/data:/input -B /path/to/results:/output \
  gvclass_1.2.2.sif /input -o /output -t 32
```

The wrapper is simpler and handles bind mounts automatically.

#### Publishing the Apptainer Image (library://)

To make `apptainer pull --library https://library.sylabs.io ...` work for users, you must build and push the SIF to the Sylabs library:

```bash
# Build the SIF from the definition file
apptainer build gvclass.sif containers/apptainer/gvclass.def

# Authenticate to the Sylabs library (one-time)
apptainer remote login

# Push the image to the library
apptainer push gvclass.sif library://nelligroup-jgi/gvclass/gvclass:1.2.2
```

### Full CLI Reference (gvclass)

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `query_dir` | | Input directory or file | Required |
| `--output-dir` | `-o` | Output directory | `<query>_results` |
| `--threads` | `-t` | Total threads | 16 |
| `--max-workers` | `-j` | Parallel workers | Auto |
| `--threads-per-worker` | | Threads per worker | Auto |
| `--database` | `-d` | Override database path | `GVCLASS_DB` → config → `<repo>/resources` |
| `--tree-method` | | `fasttree` or `iqtree` | fasttree |
| `--mode-fast` | `-f` | Fast mode: core markers only | True |
| `--extended` | `-e` | Extended mode: all marker trees | False |
| `--sensitive` | | Sensitive HMM mode (`E=1e-5`, `domE=1e-5`, skip GA cutoffs) | False |
| `--contigs` | `-C` | Treat each contig as an independent query genome | False |
| `--resume` | | Resume interrupted run | False |
| `--verbose` | `-v` | Enable debug output | False |
| `--version` | | Show version info | |
| `--cluster-type` | | `local`, `slurm`, `pbs`, `sge` | local |
| `--cluster-queue` | | HPC queue/partition name | |
| `--cluster-project` | | HPC project/account | |
| `--cluster-walltime` | | HPC time limit | 04:00:00 |

### Contig Splitting Mode (`--contigs` / `-C`)

By default, GVClass expects bins (directory input). Use `--contigs` when you want to process contigs as separate query genomes.

Primary use case:

```bash
# Split contigs and classify each independently
gvclass --contigs metagenome_contigs.fna -o results -t 32
```

**How it works:**
1. Accepts a **single multi-contig FNA file** (primary use case; directory input is also supported)
2. Splits each sequence into individual contig files in a temporary directory
3. Sanitizes contig IDs for filenames (replaces `/\:*?"<>|` and spaces with `_`)
4. Applies minimum length filter from config (`pipeline.contigs_min_length`, default: `10000` bp)
5. Processes each retained contig independently through the full pipeline
6. Combines results into `gvclass_summary.tsv`
7. Cleans up temporary files automatically

**Use cases:**
- Screening metagenome contigs before binning
- Evaluating candidate viral contigs as independent genomes
- Rapid triage when bins are not yet available

For giant virus-focused analyses, set `pipeline.contigs_min_length` to at least `30000` (preferably `50000`) for more reliable classifications.

### Genetic Code Selection

GVClass tests **9 genetic codes** to find optimal gene predictions:

- **Code 0**: Meta mode using pretrained models (pyrodigal metagenomic mode)
- **Codes 1, 4, 11**: Standard translation tables
- **Codes 6, 15, 29, 106, 129**: Additional translation tables

**Selection logic:**
1. Start with meta mode (code 0) as baseline
2. Override if another code has:
   - More complete marker hits (>66% HMM coverage), OR
   - Same hits but >5% better average hit score, OR
   - Same hits but >5% better coding density

The selected code is reported in the `ttable` output column.

### Sensitive HMM Mode (`--sensitive`)

Use sensitive mode when you want a more permissive marker search:

- Forces pyhmmer to use `E=1e-5` and `domE=1e-5`
- Disables GA/TC/NC model cutoff filtering (`--cut_ga`-style behavior)
- Can recover weaker hits, but may increase false positives

## Interpreting Quality Metrics

### Genome completeness
- `order_completeness` shows the percentage of order-specific markers detected in single copy. Values above ~70% generally indicate near-complete genomes, 30–70% partial, and below 30% highly fragmented assemblies.
- `weighted_order_completeness` applies conservation-based weights; large gaps here usually point to missing hallmark genes even if raw counts look acceptable.

### Contamination and mixed populations
- `order_dup` and `gvog8_dup` summarize marker duplication. Values above ~2 suggest multiple populations or assembly chimeras; below ~1.5 is typically clean.
- `gvog8_total` and `gvog8_unique` help distinguish true gene expansions (high total, moderate duplication) from assembly artefacts (high duplication, low uniqueness).
- `ncldv_mcp_total`, `mirus_df`, `mrya_total` provide additional lineage-specific duplication hints.
- `vp_completeness` and `mirus_completeness` show core marker coverage (n/4) for virophages and Mirusviricota respectively.
- `plv` count helps distinguish PLV from virophages (PLVs share VP markers but have additional PLV-specific marker; count is not binary).

### Cellular carry-over
- `uni56_total` (UNI56) counts universal cellular markers; more than ~10 unique hits point to host contamination or bins that include cellular contigs.
- The wrapper also reports `busco_*` fields when available; elevated counts complement UNI56 for spotting non-viral material.

Use these fields together: a high completeness score with low duplication and low UNI56 is characteristic of a high-quality GVMAG; any combination of low completeness plus high duplication or high UNI56 warrants manual curation.

## Performance Optimization

### Speed Up Analysis

1. **Enable Fast Mode** - Skip order-level marker trees (OG markers):
   ```bash
   # Command line option
   pixi run gvclass <input_directory> --mode-fast
   
   # Or in config file
   pipeline:
     mode_fast: true  # Skips ~100 order-specific markers, 2-3x faster
   ```

2. **Use FastTree Instead of IQ-TREE**:
   ```bash
   # Default (faster)
   pixi run gvclass <input_directory> --tree-method fasttree
   
   # IQ-TREE (more accurate but slower)
   pixi run gvclass <input_directory> --tree-method iqtree
   ```

3. **Optimize Thread Usage**:
   ```bash
   # Use all available cores
   pixi run gvclass <input_directory> -t 32

   # Control parallelization: 4 parallel workers, 8 threads each (= 32 total)
   pixi run gvclass <input_directory> -t 32 -j 4
   ```

### IQ-TREE Specific Options

When using IQ-TREE (`--tree-method iqtree`), the pipeline automatically uses:
- Model: `LG+F+G` (fast protein model)
- `-fast` flag for faster tree search
- Single thread per marker (parallelization happens at marker level)

To modify IQ-TREE behavior, edit `src/core/marker_processing.py`.

### Understanding Markers

- **Core markers**: Always processed (GVOG4, GVOG8, MCP, etc.)
- **Order-level markers**: 576 OG markers conserved in different viral orders
  - Processed when `mode_fast: false` (default)
  - Skipped when `mode_fast: true` (faster but less precise order assignment)

## Making GVClass Available Globally

Both wrappers can be made globally accessible:

### Option 1: Apptainer (gvclass-a) - Recommended for HPC

```bash
# Copy gvclass-a to your personal bin
mkdir -p "$HOME/bin"
cp gvclass-a "$HOME/bin/"
chmod +x "$HOME/bin/gvclass-a"

# Add to PATH (if not already)
echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
source "$HOME/.bashrc"

# Now use from anywhere!
gvclass-a /data/genomes /data/results -t 32
```

### Option 2: Symlink gvclass (for Pixi users)

```bash
# From the gvclass repo directory
cd /path/to/gvclass  # Navigate to your cloned repo first

# Create symlink (REQUIRED - copying won't work!)
mkdir -p "$HOME/bin"
ln -s "$(pwd)/gvclass" "$HOME/bin/gvclass"

# Add to PATH (if not already)
echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
source "$HOME/.bashrc"

# Now run from anywhere - symlink allows script to find its repo
cd /anywhere
gvclass my_data -o results -t 32
```

## How It Works

```mermaid
flowchart TD
    subgraph Input
        FNA[".fna<br/>nucleic acid"]
        FAA[".faa<br/>amino acid"]
    end
    
    subgraph Database
        DB[(Reference<br/>Database)]
        MODELS[GVOG HMMs<br/>+ marker sets]
        REF[Reference<br/>sequences]
    end
    
    FNA --> OPGC{Optimized<br/>Gene Calling}
    FAA --> ID1[Identify Markers]
    
    subgraph "Gene Calling (pyrodigal)"
        OPGC --> META[Meta mode<br/>code 0<br/>pretrained models]
        OPGC --> STD[Standard codes<br/>1, 4, 11]
        OPGC --> ADD[Additional codes<br/>6, 15, 29, 106, 129]
        META --> RANK[Select best by<br/>marker hits &<br/>coding density]
        STD --> RANK
        ADD --> RANK
    end
    
    RANK --> ID2[Identify Markers]
    
    subgraph "Marker Detection (pyhmmer)"
        ID1 --> HMM1[HMM search<br/>against GVOGs]
        ID2 --> HMM2[HMM search<br/>against GVOGs]
        HMM1 --> HITS1[Extract hits<br/>E-value cutoffs]
        HMM2 --> HITS2[Extract hits<br/>E-value cutoffs]
    end
    
    HITS1 --> BLAST[BLAST/pyswrd<br/>top 100 hits]
    HITS2 --> BLAST
    
    subgraph "Alignment & Trees"
        BLAST --> ALIGN[MAFFT/pyfamsa<br/>alignment]
        ALIGN --> TRIM[TrimAl/pytrimal<br/>trimming]
        TRIM --> TREE{Tree Building}
        TREE -->|FastTree| FT[veryfasttree<br/>LG4X model]
        TREE -->|IQ-TREE| IQ[iqtree<br/>LG+F+G -fast]
    end
    
    FT --> NN[Get nearest<br/>neighbors]
    IQ --> NN
    
    subgraph "Classification"
        NN --> TAX[Majority/strict<br/>taxonomy<br/>assignment]
        NN --> QC[Quality metrics:<br/>completeness,<br/>contamination]
    end
    
    TAX --> OUT[Results:<br/>taxonomy,<br/>QC metrics]
    QC --> OUT
    
    MODELS -.-> HMM1
    MODELS -.-> HMM2
    REF -.-> BLAST
    DB -.-> NN
    
    style FNA fill:#e8f4f8
    style FAA fill:#e8f4f8
    style OUT fill:#d4edda
    style DB fill:#f8d7da
    style MODELS fill:#f8d7da
    style REF fill:#f8d7da
```

## Citation

If you use GVClass, please cite:

> Pitot et al. (2024): Conservative taxonomy and quality assessment of giant virus genomes with GVClass. npj Viruses. https://www.nature.com/articles/s44298-024-00069-7

## Database References

The GVClass v1.2.1 reference database includes genomes from the following sources:

> Medvedeva S, Guyet U, Pelletier E, et al. (2026): Widespread and intron-rich mirusviruses are predicted to reproduce in nuclei of unicellular eukaryotes. Nature Microbiology 11:228-239. https://doi.org/10.1038/s41564-025-01906-2

> Roux S, Fischer MG, Hackl T, Katz LA, Schulz F, Yutin N (2023): Updated Virophage Taxonomy and Distinction from Polinton-like Viruses. Biomolecules 13(2):204. https://doi.org/10.3390/biom13020204

> Fiamenghi MB, Camargo AP, Chasapi IN, et al. (2025): Meta-virus resource (MetaVR): expanding the frontiers of viral diversity with 24 million uncultivated virus genomes. Nucleic Acids Research gkaf1283. https://doi.org/10.1093/nar/gkaf1283

> Vasquez YM, Nardi T, Terasaki GM, et al. (2025): Genome-resolved expansion of Nucleocytoviricota and Mirusviricota reveals new diversity, functional potential, and biotechnological applications. bioRxiv 2025.09.26.678796. https://doi.org/10.1101/2025.09.26.678796

## Support

- **Issues**: [GitHub Issues](https://github.com/NeLLi-team/gvclass/issues)
- **Contact**: fschulz@lbl.gov

## License

BSD 3-Clause License - see LICENSE file for details

---
<sub>Version 1.2.2 - February 2026</sub>
