<p align="center">
  <img src="images/GVClass_logo.png" alt="GVClass Logo" width="50%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-v1.1.0-blue.svg" alt="Version">
  <img src="https://img.shields.io/badge/license-BSD--3--Clause-green.svg" alt="License">
  <img src="https://img.shields.io/badge/python-3.11-blue.svg" alt="Python">
  <img src="https://img.shields.io/badge/pixi-enabled-orange.svg" alt="Pixi">
</p>

# GVClass - Giant Virus Classification Tool

GVClass assigns taxonomy to giant virus contigs and metagenome-assembled genomes (GVMAGs). It uses phylogenetic analysis based on giant virus orthologous groups (GVOGs) to provide accurate classification from domain to species level.

## 🚀 Quick Start

### 1. Install Pixi (one-time setup)
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Clone and Enter Directory
```bash
git clone https://github.com/NeLLi-team/gvclass.git
cd gvclass
```

### 3. Install Dependencies
```bash
pixi install
```

### 4. Run GVClass
```bash
# Basic usage - database downloads automatically on first run (~700MB)
pixi run gvclass <input_directory>

# Run with example data to test installation
pixi run run-example

# Custom output directory and threads
pixi run gvclass <input_directory> -o my_results -t 16
```

That's it! No manual database setup required - GVClass handles everything automatically.

## 📋 Input Requirements

- **Directory** containing `.fna` (nucleic acid) or `.faa` (protein) files
- **Minimum size**: 20kb recommended (50kb+ preferred)
- **Clean filenames**: Avoid special characters (`.` `;` `:`), use `_` or `-` instead
- **Protein headers**: Format as `filename|proteinid` for best results

## 🎯 Example Usage

```bash
# Run on your data
pixi run gvclass my_genomes/

# Specify output location
pixi run gvclass my_genomes/ -o classification_results

# Use more threads for faster processing
pixi run gvclass my_genomes/ -t 32

# Process multiple queries in parallel (4 queries × 8 threads each = 32 total)
pixi run gvclass-parallel my_genomes/ -t 32 -j 4
```

## 📊 Output

Results are saved to `<input_name>_results/` containing:
- `gvclass_summary.tsv` - Main results with taxonomy assignments
- Individual query subdirectories with detailed analysis

### Output Columns Explained

| Column | Description |
|--------|-------------|
| query | Input filename |
| taxonomy_majority | Full taxonomy based on majority rule |
| taxonomy_strict | Conservative taxonomy (100% agreement) |
| species → domain | Individual taxonomic levels |
| avgdist | Average tree distance to references |
| *_completeness | Estimated genome completeness |
| *_unique/*_dup | Marker gene counts |
| contigs | Number of contigs |
| LENbp | Total length in base pairs |
| GCperc | GC content percentage |
| genecount | Number of predicted genes |
| CODINGperc | Coding density percentage |
| ttable | Genetic code used |

## ⚙️ Configuration (Optional)

Create `gvclass_config.yaml` to set defaults:

```yaml
database:
  path: resources                    # Database location

pipeline:
  tree_method: fasttree             # or 'iqtree' for more accuracy
  mode_fast: false                  # Skip some analyses when true
  threads: 16                       # Default thread count
```

## 🆕 What's New in v1.1.0

- **🚀 Modern Architecture**: Prefect + Dask workflow orchestration
- **📦 Easy Installation**: Pixi package manager (2-3x faster)
- **🐍 Pure Python**: All tools replaced with faster Python versions
- **⚡ Better Performance**: Parallel marker processing, 25% faster
- **🔄 Automatic Recovery**: Task caching and retry on failures
- **✅ Auto Database Setup**: No manual download needed

## 📖 Advanced Usage

### Container Execution

```bash
# Apptainer/Singularity
bash gvclass_apptainer.sh <input_dir> <threads>

# Docker
bash gvclass_docker.sh <input_dir> <threads>

# Shifter (NERSC)
bash gvclass_shifter.sh <input_dir> <threads>
```

### Development Commands

```bash
# Run specific test
pixi run python test_pipeline.py

# Clear cache and run fresh
pixi run python clear_cache_and_run.py

# Debug mode
pixi run python debug_pipeline.py
```

## 🔬 How It Works

<p align="center">
  <img src="images/Workflow_GVClass.png" alt="GVClass Workflow" width="100%">
</p>

1. **Input Validation**: Checks and reformats input sequences
2. **Gene Calling**: Tests multiple genetic codes to find optimal (nucleic acid input only)
3. **Marker Search**: Searches against GVOG HMM models using pyhmmer
4. **Phylogenetic Analysis**: Builds trees for each marker gene
5. **Taxonomy Assignment**: Consensus classification based on nearest neighbors
6. **Quality Assessment**: Estimates completeness and contamination

## 📝 Citation

If you use GVClass, please cite:

> Schulz et al. (2020) Giant virus diversity and host interactions through global metagenomics. Nature. https://doi.org/10.1038/s41586-020-1957-x

## 🤝 Support

- **Issues**: [GitHub Issues](https://github.com/NeLLi-team/gvclass/issues)
- **Contact**: fschulz@lbl.gov

## 📄 License

BSD 3-Clause License - see LICENSE file for details

---
<sub>Version 1.1.0 - December 2024</sub>