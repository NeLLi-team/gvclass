# Changelog

All notable changes to GVClass are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/), and this
project follows [Semantic Versioning](https://semver.org/). The entries below are
condensed; full per-release notes live on the
[GitHub Releases page](https://github.com/NeLLi-team/gvclass/releases), which also
covers releases prior to v1.2.2.

The runtime database/resource bundle is versioned separately from the software; the
compatible bundle is noted per release.

## [Unreleased]

## [1.6.1] - 2026-05-31

Compatible resource bundle: v1.6.0 (Zenodo DOI 10.5281/zenodo.20479524).

### Fixed
- `taxonomy_majority` now builds consensus top-down from full neighbor lineages:
  after a domain wins, each lower rank is voted only among hits under the
  already selected parent lineage. This prevents minority EUK or off-parent
  lower-rank labels from leaking into MIRUS/NCLDV calls.
- `taxonomy_strict` remains no more specific than the accepted majority path
  and blanks lower ranks once an upstream parent rank is unresolved.
- Taxonomy support-count columns now render placeholder and grouped labels
  consistently as `DOMAIN__unclassified` and `DOMAIN__other` for non-domain
  ranks, while domain counts remain plain `DOMAIN`. This removes output strings
  such as `NCLDV__NCLDV_unclassified`, `NCLDV_other`, and `PLV_other`.
- The colorized terminal banner frame now aligns correctly by computing padding
  from visible text after stripping ANSI escape sequences.

### Changed
- Default database download configuration now points to the harmonized
  `resources_v1_6_0.tar.gz` Zenodo record and SHA256.
- Container and wrapper references now target the `1.6.1` software image.

### Added
- Reusable label harmonization script and resource-integrity tests that compare
  FAA headers against active labels, with alias support available only when
  explicitly enabled.

### CI
- Hardened the Apptainer SIF smoke test: it now runs the bare `gvclass --version` via
  `PATH` and asserts `command -v gvclass` resolves to `/usr/local/bin/gvclass`, so CI
  exercises the exact entrypoint Snakemake/Nextflow use (the `/usr/local/bin/gvclass`
  symlink + non-root exec permissions) instead of the absolute `/opt/gvclass/gvclass`
  path. A regression of the container command-publishing fix now fails CI. The
  container-build job is tag-gated (`refs/tags/v*`), so this guards each release tag.
- Added a tag-gated, secret-guarded step that auto-publishes the SIF to the Sylabs
  library (`library://nelligroup-jgi/gvclass/gvclass:<version>`) after the artifact
  upload (a publish failure never loses the SIF) and self-skips when the `SYLABS_TOKEN`
  repo secret is absent (e.g. forks).

## [1.6.0] - 2026-05-29

Correctness, scientific-output, and security release. Compatible resource bundle:
v1.5.0. Some fixes change reported values (completeness, completeness quality,
per-contig contamination) — re-run analyses where those columns matter.

### Fixed
- Novelty-completeness R² gate now fires (reads the shipped `test_r2` column) instead
  of a never-present metadata key; `estimated_completeness_r2_holdout` is populated at
  full precision.
- Completeness no longer collapses to 0% on fallback (uses the strategy-1 estimate).
- `taxonomy_strict` is now always at least as conservative as `taxonomy_majority`.
- Per-contig contamination counts pipeless phage references (their genome ids were
  silently dropped), changing `viral_bearing_contig_count` and downstream signals.
- Gene calling: DNA `.fasta`/`.fas` inputs are gene-called; minus-strand `.fna` records
  are reverse-complemented; genetic-code selection gains a 5% coding-density override gate.
- Weighted completeness reads a real per-order marker-conservation table instead of an
  HMM size table (activates once a bundle ships `markers/marker_conservation.tsv`).
- BLAST honours the documented top-100 hits per query (pyswrd defaulted to 10);
  `--resume` reports skipped queries with missing summaries instead of dropping them;
  partial/empty config files no longer crash (defaults deep-merged); added
  `ErrorHandler.log_error`; Docker compose/example commands invoke the entrypoint correctly.

### Security
- Database downloads fail closed: an env-provided `GVCLASS_DB_URL` requires
  `GVCLASS_DB_SHA256`, sources must be `https`, and unpinned legacy NERSC mirrors were
  removed.
- `completeness/model.joblib` is now SHA-256-gated like the contamination model.

## [1.5.2] - 2026-05-27

### Fixed
- Apptainer/Singularity container entrypoint for workflow engines: publishes a stable
  `gvclass` command at `/usr/local/bin/gvclass`, prepends `/usr/local/bin` to `PATH`
  ahead of the Pixi environment, and normalizes `/opt/gvclass` permissions so non-root
  Apptainer/Singularity users can exec `gvclass`/`gvclass-a`. Compatible bundle: v1.5.0.

## [1.5.1] - 2026-05-18

### Fixed
- `--resume` rebuilds `gvclass_summary.tsv`/`.csv` from both newly processed and
  previously completed (skipped) queries (#13); completed queries now also write
  `<query>.final_summary.tsv` so future resumes preserve the full summary schema.
  Covered for `.faa`, direct `.fna`, all-skipped, and `--contigs` modes. Compatible
  bundle: v1.5.0 (no database update required).

## [1.5.0] - 2026-04-21

Retrained contamination model (simulated giant-virus MAGs) that eliminates the
novel-virus false-positive class, plus reorganised runtime resources.

### Added
- Per-contig taxonomic-purity classifier with new columns
  (`cellular_coherent_contig_count`, `cellular_lineage_purity_median`,
  `viral_bearing_contig_count`, `contig_attribution_mode`, …); novel viruses with
  scattered HGT markers no longer downgrade to `mixed_viral`.
- `taxonomy_confidence` column; opt-in `--allow-short`.

### Changed
- Contamination model retrained (ExtraTrees) on 200 simulated MAG bins + 50 clean
  isolates: mean predicted contamination on clean shredded NCLDV 28.8% → 4.3%,
  contaminated-bin MAE 9.4% → 2.55%, Pearson r ~0 → 0.97; external holdout max 0.00%.
- Reorganised `resources/` (`hmm/`, `markers/`, `completeness/`, `contamination/`,
  `database/`); tarball 3.2 GB → 1.7 GB. Removed a redundant per-query BLAST pass.
- Module renames: `prefect_flow.py` → `parallel_runner.py`,
  `gvclass_prefect.py` → `gvclass_runner.py`.

### Removed
- `database/dmnd/` (the pipeline uses `pyswrd` directly on `.faa` references).

### Security
- Bundled contamination model replaced and its SHA-256 constant rotated (loads only
  when the on-disk digest matches).

### Breaking
- Resources layout reorganised, so v1.4.x tarballs no longer pass validation — use
  `pixi run setup-db` to install the v1.5.0 bundle (Zenodo DOI 10.5281/zenodo.19674504).
  Internal module renames affect direct importers.

## [1.4.2] - 2026-03-26

### Changed
- Recalibrated the contamination estimate with a new ExtraTrees model (threshold 3.0)
  trained on the real-contig benchmark subset; moved the bundled model into
  `src/bundled_models/` so contamination-model updates no longer require a new
  `resources/` archive. Completeness bundle stays v1.4.0.

## [1.4.1] - 2026-03-13

### Fixed
- Collapse duplicate `(protein_id, model)` HMM hits before writing filtered outputs,
  correcting GVOG and category-level marker summary totals.

### Changed
- Sensitive pyhmmer filtering (`E=1e-5` / `domE=1e-5`) is now the default runtime mode
  (override with `sensitive_mode: false`). Compatible bundle: v1.4.0.

## [1.4.0] - 2026-03-09

### Added
- Novelty-aware completeness estimator (`--completeness-mode legacy|novelty-aware`)
  with new `order_completeness_v2*` and `estimated_completeness*` columns plus OOD flags.
- Trained contamination model bundle (`hist_gbm_v1`) replacing the rule-based score.

### Changed
- Runtime resources are distributed as `resources_v1_4_0.tar.gz` (Zenodo) and installed
  into a local `resources/` instead of the git tree; the final summary exposes single
  primary completeness/contamination fields while retaining duplication factors.

## [1.2.2] - 2026-02-17

### Fixed
- HPC/runtime failures caused by ephemeral Prefect API startup in restricted
  environments; simplified the execution path to run pipeline tasks directly while
  preserving parallel query processing.

### Changed
- Software/CLI and container image references updated to 1.2.2 (database remains
  v1.2.1). See GitHub Releases for notes on v1.2.1 and earlier.
