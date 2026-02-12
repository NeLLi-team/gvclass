# Changelog

All notable changes to GVClass will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.2.1] - 2026-02-12

### Highlights
- Updated bundled reference database to `resources_v1_2_1.tar.gz`
- Added models for Polinton-like virus MCP markers
- Extended and updated the eukaryotic reference database
- Updated labels file to `gvclassFeb26_labels.tsv`
- Added sensitive HMM mode (`--sensitive`) to run pyhmmer with `E=1e-5` / `domE=1e-5` instead of GA cutoffs
- Streamlined CLI and pipeline code paths for clearer flow and more consistent runtime error handling

### Database Updates
- **Database Version**: Added `resources_v1_2_1.tar.gz` as the primary release database
- **New Models**: Added Polinton-like virus MCP models
- **Reference Expansion**: Extended and refreshed eukaryotic references
- **Labels Update**: Replaced `gvclassJan26_labels.tsv` with `gvclassFeb26_labels.tsv`

### Versioning and Distribution
- Updated software version strings and wrappers to `v1.2.1`
- Updated Apptainer image defaults to `library://nelligroup-jgi/gvclass/gvclass:1.2.1`

## [v1.2.0] - 2026-01-08

### Highlights
- Updated reference database v1.2.0 with refined GA thresholds and model annotations
- Improved HMM model accuracy through updated gathering thresholds
- Added functional annotations to HMM models
- Bug fixes and documentation improvements

### Database Updates
- **Expanded Reference Database**: Added new reference genomes for improved classification coverage
  - Nucleocytoviricota and Mirusviricota genomes from Vasquez et al. (2025) bioRxiv
  - Mirusviricota models and genomes from Medvedeva et al. (2026) Nature Microbiology
  - Models for Polinton-like viruses (PLV) and virophages (PV) from Roux et al. (2023) Biomolecules
  - Extended VP, PLV and phage reference set from MetaVR database (Fiamenghi et al. 2025)
- **Updated GA Thresholds**: Refined gathering (GA) thresholds in HMM models for more accurate marker detection
- **Model Annotations**: Added functional annotations to HMM models for better interpretability
- **Database Version**: New `resources_v1_2_0.tar.gz` with all updates

### New Output Columns
- **VP (Virophage) Metrics**:
  - `vp_completeness`: Core marker completeness (n/4) based on MCP, Penton, ATPase, Protease
  - `vp_mcp`: Count of proteins with VP MCP marker hits
  - `vp_df`: Virophage duplication factor (total VP hits / 4)
- **PLV (Polinton-like Virus) Metrics**:
  - `plv`: Count of proteins with PLV marker hits (PLVs share VP markers but have additional PLV marker)
- **Mirus (Mirusviricota) Metrics**:
  - `mirus_completeness`: Core marker completeness (n/4) based on MCP, ATPase, Portal, Triplex
  - `mirus_df`: Mirusviricota duplication factor
- **NCLDV MCP**:
  - `ncldv_mcp_total`: NCLDV-specific MCP marker count (includes OG1352, OG484)
- **Removed**: Old columns `mirus_unique`, `mirus_total`, `mirus_dup` replaced by category-based completeness

### CLI Improvements
- **Renamed**: `--no-mode-fast` flag renamed to `--extended` / `-e` for clearer semantics
  - Use `-e` or `--extended` to build trees for all markers (more comprehensive analysis)
  - Default behavior remains fast mode (core markers only)

### Bug Fixes
- **Fixed**: Removed `time.sleep()` usage in `prefect_flow.py` that violated coding standards
- **Fixed**: Inconsistent bytes encoding in `genetic_code_optimizer.py` - now uses consistent string encoding for pyrodigal

### Documentation
- **Added**: Full CLI reference table documenting all command-line options
- **Added**: Detailed documentation of `--contigs` mode behavior (file splitting, ID sanitization)
- **Added**: Genetic code selection logic explanation (9 codes tested, selection criteria)
- **Updated**: Mermaid diagram now correctly shows code 0 (meta mode with pretrained models)
- **Removed**: Emoticons from README for cleaner presentation

### Code Quality
- Simplified directory creation logic (removed unnecessary retry loop)
- Consistent sequence encoding in gene calling module
- Cleaner, more maintainable codebase

---

## [v1.1.1] - 2025-10-01

### Highlights
- Switched the GVClass pipeline to Prefect + Dask orchestration (formerly Snakemake) for resilient scheduling and parallelism.
- Refreshed the bundled reference database to `resources_v1_1_1.tar.gz`, delivering corrected eukaryotic taxonomy strings and updated giant virus placements ([Schulz et al., 2025](https://doi.org/10.1101/2025.09.26.678796)).
- Updated README, pixi tasks, and developer tooling to align with the Prefect-based workflow.

## [v1.1.0] - 2024-12-XX

### 🚀 Major Changes

#### Complete Pipeline Modernization with Prefect + Dask
- **NEW**: Migrated from Snakemake to Prefect 2.0 for workflow orchestration
- **NEW**: Integrated Dask for distributed parallel processing
- **NEW**: Task caching and checkpointing for robust pipeline execution
- **NEW**: Rich monitoring and logging throughout the pipeline
- **NEW**: Automatic retry logic for failed tasks

#### Dependency Management Migration
- **BREAKING**: Migrated from conda to pixi for faster, more reliable dependency management
- **NEW**: Added `pixi.toml` configuration file with locked dependencies
- **NEW**: Simplified installation process with single `pixi install` command
- **IMPROVED**: 2-3x faster dependency installation compared to conda

#### Bioinformatics Tools Modernization
- **REPLACED**: HMMER with pyhmmer v0.10.15 for HMM searches (native Python, 2-3x faster)
- **REPLACED**: trimAl with pytrimal v0.8.0 for alignment trimming
- **REPLACED**: MAFFT with pyfamsa v0.5.3 for sequence alignment
- **REPLACED**: Diamond with pyswrd v0.1.0 for sequence similarity searches
- **REPLACED**: FastTree with VeryFastTree v4.0.3 for phylogenetic tree construction
- **IMPROVED**: All tools now use Python APIs for better integration and performance
- **UPDATED**: GVMAGs v2 reference database with refreshed orthogroups and taxonomy annotations
- **ADDED**: Support for both IQ-TREE and VeryFastTree methods
- **ADDED**: Custom genetic code support for codes 106 and 129 (giant virus-specific)

### ✨ New Features

- **Prefect-based Workflow**: Modern workflow orchestration with task dependencies and caching
- **Parallel Processing**: Dask-powered parallel execution of marker-specific analyses
- **Enhanced Error Handling**: Comprehensive error handling framework with structured logging
- **Input Validation**: Robust input validation for all parameters and file formats
- **Full Taxonomy Output**: Complete taxonomic classification from domain to species level
- **Improved Stats**: Better sequence statistics for both nucleotide and protein inputs
- **Marker-specific Processing**: Each HMM marker processed independently for better parallelization

### 🛠️ Technical Improvements

#### Security Enhancements
- **FIXED**: Command injection vulnerabilities in subprocess calls
- **ADDED**: Input validation for shell scripts
- **IMPLEMENTED**: Path traversal protection
- **ENHANCED**: Resource management with proper context managers

#### Code Quality
- **ADDED**: Comprehensive error handling with custom exception classes
- **IMPROVED**: Type hints and documentation throughout codebase
- **ENHANCED**: Logging with structured context and debugging information
- **STANDARDIZED**: Code style and security best practices

#### Dependencies
- **UPDATED**: Complete Python-based tool ecosystem
- **ADDED**: Prefect 2.0 for workflow orchestration
- **ADDED**: Dask and prefect-dask for distributed computing
- **ADDED**: pyhmmer 0.10.15 for HMM searches
- **ADDED**: pytrimal 0.8.0 for alignment trimming  
- **ADDED**: pyfamsa 0.5.3 for sequence alignment
- **ADDED**: pyswrd 0.1.0 for sequence similarity searches
- **ADDED**: VeryFastTree 4.0.3 for fast phylogenetic tree construction
- **ADDED**: pyrodigal 3.5.2 for gene calling with custom support for genetic codes 106 and 129
- **ADDED**: rich for enhanced CLI output
- **REMOVED**: External binary dependencies (HMMER, MAFFT, trimAl, Diamond)
- **MAINTAINED**: IQ-TREE as optional tree method

### 📖 Documentation

- **UPDATED**: Installation instructions for pixi
- **ADDED**: Migration guide for existing users
- **ENHANCED**: Development documentation in CLAUDE.md
- **CREATED**: Comprehensive security improvements documentation
- **ADDED**: Migration summary with detailed change explanations
- **ADDED**: Troubleshooting guide for common installation issues
- **CREATED**: Installation fix script for problematic systems

### 🔄 Migration Guide

#### For New Users
```bash
# Install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and setup
git clone https://github.com/NeLLi-team/gvclass.git
cd gvclass
pixi install

# Download database (only needed once)
pixi run setup-db

# Run GVClass
pixi run gvclass example/ output/ --threads 24
```

#### For Existing Users
```bash
# Update repository
git pull

# Install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Install dependencies
pixi install

# New command syntax (Prefect-based)
# OLD: snakemake -j 24 --use-conda --config querydir="example"
# NEW: pixi run gvclass example/ output/ --threads 24
```

### 🐛 Bug Fixes

- **FIXED**: Resource leaks in file operations
- **FIXED**: Potential security vulnerabilities in input handling
- **FIXED**: Error handling edge cases
- **IMPROVED**: Memory management throughout pipeline

### ⚠️ Breaking Changes

1. **Dependency Management**: 
   - Conda environment is no longer used
   - Must install pixi and run `pixi install`
   - Commands must use `pixi run` prefix

2. **Workflow System**:
   - Snakemake replaced with Prefect + Dask
   - New command-line interface: `pixi run gvclass`
   - Parallel marker processing with automatic load balancing

3. **Tool Replacements**:
   - All external binaries replaced with Python implementations
   - HMMER → pyhmmer (2-3x faster, native Python)
   - trimAl → pytrimal (same algorithm, Python API)
   - MAFFT → pyfamsa (comparable performance)
   - Diamond → pyswrd (better integration)
   - FastTree → VeryFastTree (faster, Python-based)

3. **Container Changes**:
   - Container scripts updated to use pixi internally
   - No change required for container users

### 📊 Performance Improvements

- **Faster Installation**: Pixi installs dependencies 2-3x faster than conda
- **Better Integration**: Python-based tools eliminate subprocess overhead
- **Reduced Memory Usage**: More efficient memory management
- **Improved Parallel Processing**: Better thread utilization

### 🧪 Testing

- **Validated**: All functionality against example datasets
- **Tested**: Backward compatibility for output formats
- **Verified**: Container execution workflows
- **Confirmed**: Performance improvements in real-world scenarios

### 📋 File Changes

#### New Structure
- `src/` - Complete rewrite with modular Python architecture
  - `src/pipeline/` - Prefect workflow definitions
  - `src/core/` - Core analysis modules
  - `src/utils/` - Utilities and error handling
  - `src/bin/` - Command-line interfaces
- `pixi.toml` - Pixi dependency configuration
- `pyproject.toml` - Python package configuration

#### Modified Files
- `gvclass_*.sh` - Updated to use new pipeline
- `README.md` - Complete overhaul with new usage
- Container scripts updated for pixi

### 📞 Support

For questions about migration or new features:
- Check the migration guide in `MIGRATION_SUMMARY.md`
- Review updated documentation in `README.md`
- Report issues on GitHub

---

## [v1.0.0] - 2024-07-08

### Initial Release

- Giant virus classification using phylogenetic analysis
- Support for Docker, Apptainer, and Shifter containers
- Conda-based dependency management
- HMMER, MAFFT, and trimAl integration
- Comprehensive taxonomy assignment pipeline
- Quality assessment and contamination detection
