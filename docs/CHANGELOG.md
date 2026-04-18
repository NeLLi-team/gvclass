# Changelog

All notable changes to GVClass will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.5.0] - 2026-04-18

### Highlights
- **Contamination model retrained on simulated MAGs.** Replaces the bundled
  ExtraTrees contamination regressor with one trained on a dataset of
  giant-virus MAG simulations (isolate genomes fragmented into MAG-like
  contigs with added bacterial, eukaryotic, NCLDV-host eukaryotic,
  and mixed-kingdom contamination) plus clean intact isolate genomes.
  Eliminates the novel-virus false-positive class seen on prior releases:
  mean predicted contamination on clean shredded NCLDV drops from ~29%
  to ~4%, while calibration on curated intact inputs is preserved.
  Contaminated-bin MAE improves 3–5× across scenarios; Pearson r on
  contaminated bins reaches 0.97.
- Closed three load-bearing correctness bugs surfaced by the v1.4.2 code
  review (multi-HMM dedup bypass, contamination model calibration under
  sensitive mode, resume accepting corrupt runs).
- **Per-contig taxonomic-purity classifier** added to the contamination
  pipeline so novel giant-virus bins are no longer falsely flagged as
  contaminated when their markers tree-resolve to scattered eukaryotic
  HGT-recipients. A contig is only called `cellular_coherent` when it
  carries ≥3 cellular-leaning proteins, zero viral markers, ≥60%
  order-level agreement, AND median cellular-hit BLAST identity ≥70%.
  Bins with the novel-virus signature (no cellular_coherent contigs +
  ≥3 viral_bearing contigs + viral_mixture rule-based source)
  downgrade from `mixed_viral` to `uncertain` so curators can triage.
- Added a SHA-256 gate on the bundled contamination model plus a model card
  so any swap of the joblib fails fast and goes through code review.
- Added an opt-in `--allow-short` flag and hardened input validation to
  fail-closed on short / malformed FASTA.
- Added `src/__version__.py` as the single source of truth, `pyproject.toml`
  for `pip install -e .`, committed `pixi.lock`, and introduced
  `.github/workflows/ci.yml` with lint/test + PR golden-file regression.

### Fixed (correctness)
- **Multi-HMM dedup**: `src/pipeline/query_processing_engine.py::_run_hmm_search`
  previously copied the raw per-domain HMM output verbatim into
  `models.out.filtered`. A new `dedup_domtbl_file` helper
  (`src/core/hmm_search_multi.py`) collapses the file to one row per
  `(protein, model)` pair, matching the single-HMM path. Domtbl rows with
  fewer than the emitted 22 fields are now rejected to prevent malformed
  input from slipping through the dedup pipeline.
- **Contamination model retrained on simulated MAGs**: final v1.5.0
  bundle (`src/bundled_models/contamination_model.joblib`,
  `training_profile: simulated_mag_plus_intact_v1_5_0`) is fitted on a
  250-row training set — 200 simulated giant-virus MAG bins
  (isolate genomes fragmented to MAG-like contigs with added
  bacterial, eukaryotic, NCLDV-host eukaryotic, and mixed-kingdom
  contamination) plus 50 clean intact isolate genomes. Selected
  estimator is `ExtraTreesRegressor` (`model_name: extra_trees`).
  Per-scenario training MAE: clean_intact 0.00%, clean_mag_canonical
  3.80%, clean_mag_hgt_rich 3.45%, cellular_bact 1.80%, cellular_euk
  2.08%, adversarial_euk_host 3.28%, mixed_cellular 2.23%, viral_mix
  3.32%. Contaminated Pearson r = 0.97. The bundled model is safe to
  apply under both sensitive and non-sensitive runs, and
  `FullSummarizer` no longer needs the transitional NaN /
  `uncertain_sensitive_mode` gate.
- **Resume reliability**: per-query tarballs are written atomically
  (`.tar.gz.part` + `is_tarfile` verification + `os.replace` + parent
  `fsync`). A new JSON `*.SUCCESS` sentinel (with summary/tar SHA-256,
  software version, completion timestamp) is the primary resume marker;
  the legacy `(summary.tab, tar.gz)` pair is accepted only if the tar
  passes `tarfile.is_tarfile`. On reruns, `_clear_prior_outputs` wipes
  every resumable artefact so a crash mid-regeneration cannot leave a
  silently-skippable legacy pair on disk. Database install also extracts
  into `<db_path>.new` and atomically swaps, with rollback on failure.
- **Taxonomy majority per marker**: `summarize_full.py::_collect_taxonomy_counts`
  now tracks a per-marker breakdown; the consensus rule uses one-vote-
  per-marker with a minimum number of supporting markers per level
  (order=3 normal, 2 fast-mode; others=2) and emits a new
  `taxonomy_confidence` column. Deterministic alphabetic tiebreak ensures
  two runs on the same data cannot disagree.
- **Marker extraction**: `marker_extraction.extract_marker_sequences` now
  indexes by both `record.id` and `record.description` so headers with
  trailing annotation (pyrodigal-style `"prot_1 # start # end"`) resolve.
- **Genetic-code deterministic tiebreak**: `select_best_code` uses an
  explicit `_CODE_PREFERENCE_RANK` map (11, 1, 4, 0, 6, 15, 29, 106, 129).
  A non-meta code must clear an absolute margin over meta (>=2 extra
  complete hits or 10% higher avg best-hit score) to override. Prevents
  noise-driven flips to exotic codes.
- **Thread coordination**: `build_runtime_env` now seeds
  `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS` at the per-marker budget, and
  `veryfasttree.run` actually receives the `threads` argument (previously
  single-threaded regardless of allocation).
- **Contig splitter**: now transactional. `split_contigs` and the
  directory-mode wrapper detect sanitized filename collisions in a plan
  phase before any files are written; a collision leaves the output dir
  empty.
- **`.fasta` / `.fas` content validation**: extension-agnostic inputs now
  route through content-alphabet inference so DNA/protein character
  checks run. Directory validation (`validate_query_directory`) re-raises
  per-file errors instead of the previous log-and-continue.
- **Novelty R² gating**: predicted completeness is now gated on
  `r2_holdout`. Below `0.5` the strategy-2 tier score is the primary
  estimate and the ML prediction surfaces only as
  `estimated_completeness_advisory`. New `estimated_completeness_quality`
  column (high / moderate / advisory_only) and `estimated_completeness_r2_holdout`.
  NaN r2 values are treated as missing.

### Fixed (security)
- SHA-256 gate on `src/bundled_models/contamination_model.joblib`. The
  constant (`CONTAMINATION_MODEL_SHA256`) is in Python source rather than
  a sidecar so rotation requires a `.py` change. Verification hashes and
  deserializes the same in-memory bytes to close the TOCTOU window
  between the hash and `joblib.load`. New rotation helper at
  `scripts/rotate_contamination_model.py`.

### Added
- `src/__version__.py` as the single source of truth.
- `pyproject.toml` with console script `gvclass = src.bin.gvclass_cli:main`.
- `.github/workflows/ci.yml`: lint + type-check + pytest on push/PR;
  golden-file regression (`tests/fixtures/expected_example_summary.tsv`)
  on PR with DB caching; container build + SBOM on tag push.
- `docs/quality_metrics.md` (now tracked) and model-card YAML at
  `src/bundled_models/contamination_model.yaml`.
- `--allow-short` CLI flag for below-minimum inputs.
- New output columns: `taxonomy_confidence`,
  `estimated_completeness_quality`, `estimated_completeness_advisory`,
  `estimated_completeness_r2_holdout`.
- New per-contig taxonomic-purity output columns:
  `cellular_coherent_contig_count`,
  `cellular_coherent_protein_fraction`,
  `cellular_coherent_bp_fraction`,
  `cellular_lineage_purity_median`,
  `cellular_hit_identity_median`,
  `viral_bearing_contig_count`,
  `contig_attribution_mode`.

### Changed
- `.gitignore` no longer masks `pixi.lock` or `docs/` (plans directory
  still ignored). `pixi.lock` is now committed.
- `resolve_completeness_mode` fallback matches `default_config`
  (`novelty-aware`), closing the documented/actual drift.
- README resolution claim walked back from "domain to species level" to
  "domain to family level; genus / species only when a near-identical
  reference exists".
- Query enumeration now accepts `.fasta` / `.fas` in addition to `.fna` / `.faa`.

### Removed
- Dead Prefect deployment helpers (`create_local_deployment`,
  `create_slurm_deployment`) from the orchestration module.
- Redundant per-query BLAST pass. The engine no longer emits
  `*.blastpout` files; the canonical `<marker>.m8` BLAST output from
  `MarkerProcessor` is the single source of truth for both
  contamination scoring and tree building. Halves per-query BLAST cost.

### Renamed
- `src/pipeline/prefect_flow.py` → `src/pipeline/parallel_runner.py`.
- `src/bin/gvclass_prefect.py` → `src/bin/gvclass_runner.py`. The
  module never used Prefect at runtime, so the rename removes the
  misleading dependency implication. The CLI subprocess spawn path,
  all internal imports, and logger names were updated accordingly.

### Deferred to v1.5.0
- Tree-NN quality filter (min tips, bootstrap / patristic gate).
- MRYA generic `ATPase` HMM retightening.

### Notes
- Software release is `v1.5.0`; runtime resource bundle remains `v1.4.0`.
- Test coverage grew from 44 passing tests at `v1.4.2` to 107 passing
  (+1 opt-in golden-file) on the v1.5.0 branch — 63 new
  pytest cases across ~64 new test functions. The golden-file test runs
  only when `GVCLASS_RUN_GOLDEN=1` (matches the PR-only CI job).
- Verified end-to-end: `pixi run gvclass example -t 8 --mode-fast`
  completes in 43s; three queries report `taxonomy_confidence = high`.

## [v1.4.1] - 2026-03-12

### Changed
- Enabled sensitive pyhmmer filtering (`E=1e-5` / `domE=1e-5`) by default in the
  shipped configuration and deployment templates.

### Fixed
- Collapsed repeated hits from the same protein to the same marker model before
  writing filtered HMM outputs and marker count tables.
- Corrected downstream marker summaries so GVOG and related totals inherit the
  deduplicated per-protein-per-model counts.

### Notes
- The GVClass software release is `v1.4.1`.
- The bundled runtime resource bundle remains `v1.4.0`.

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
