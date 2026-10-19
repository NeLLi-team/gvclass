# Changelog

All notable changes to GVClass are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/), and this
project follows [Semantic Versioning](https://semver.org/). The entries below are
condensed; full per-release notes live on the
[GitHub Releases page](https://github.com/NeLLi-team/gvclass/releases), which also
covers releases prior to v1.2.2.

The runtime database/resource bundle is versioned separately from the software; the
compatible bundle is noted per release.

## [2.1.0] - 2026-10-14

Compatible resource bundle: v2.1.0 (Zenodo DOI 10.5281/zenodo.23381313).

### Added

- `gvclass screen`: a cellular screen (bacterial, archaeal and eukaryotic sequence in a
  giant-virus bin) and a giant-virus screen (contigs of a second giant virus). The
  resources (a 10% protein reference, fitted models and ESM-2 8M weights) install
  with `gvclass screen --install`. Torch and transformers join the pixi environment.
  Guide: `docs/how-to/screen-contigs.md`.
- The database installer downloads archives split into parts, checks each part and
  the joined archive, and installs optional archives next to the database.

### Changed

- Unified `gvclass` and `gvclass-a` into a single `gvclass` entrypoint supporting
  both local execution and Apptainer containers (`--apptainer`, `--image`, and
  standalone `wget` fallback). Container mode forwards the analysis options;
  `--config`, `--database` and the cluster options work only in native mode.
  `gvclass-a` is removed from `main`: a job script that downloads it must
  download `gvclass` and keep its `--image` or `--apptainer` argument. In
  container mode `-t` defaults to 16 and the output directory defaults to
  `<input>_results` beside the input.
- Consolidated `src/test/` into `tests/`, moved benchmark helpers to `scripts/`,
  and removed legacy `src/bin/gvclass_*.sh` wrappers.
- Resource v2.1.0 removes 2,246 prokaryotic proteins (3,083 marker rows) from 319
  eukaryotic reference genomes. A screen of 49,260 single-cell assemblies named
  them by near-identical placement of confident prokaryotic genomes; no reference
  genome was removed. Prokaryotic genomes had matched these references and
  received eukaryotic or unassigned domain calls. On a panel of 44 genomes, of 9
  gut prokaryote genomes the old resource gives 3 eukaryotic calls, 5 unassigned
  and 1 bacterial; the new one gives 8 bacterial and 1 eukaryotic. Three bacterial
  genomes, one archaeal genome and 31 eukaryotic genomes (19 protist single cells,
  12 eukaryotic MAGs) keep their calls, and completeness and contamination do not
  change. A second screen of the bundle (2,190 more proteins with a near-identical
  bacterial or archaeal hit and no close eukaryotic homolog) changed no call on
  that panel and is not applied; its list is not part of the archive. Marker
  conservation tables and completeness baselines were not recomputed.
- A eukaryotic domain is now assigned only when its markers are at least 80% of the
  markers that vote for a cellular domain. Before, a plurality of two markers was
  enough, and 59 of 1,470 prokaryotic gut genomes received a eukaryotic call
  because distant homologs and contaminated eukaryotic references outvoted the
  correct domain (5 remain; 12 of 12 eukaryotic single cells keep their call). A
  domain that fails the test is left unassigned and `taxonomy_confidence`
  reports `mixed_cellular`. Bacterial, archaeal, and viral calls are unchanged.
- Bundled models report support explicitly. A completeness estimate for an input
  with low reference or query support is no longer floored at the higher of the
  two baseline estimates and can fall to the lower one. An input without an
  informative marker reports `NaN` and `completeness_query_support`
  `insufficient_evidence`, not 0. The bundled contamination model reports only for
  targets accepted as NCLDV with supported evidence; MIRUS, PPV and other targets,
  and inputs without evidence, report `NaN` and name the reason in
  `contamination_model_applicability`. The labels `moderate` and `high` are no
  longer emitted: the bundled validation cannot support them, so model-backed
  values read `unvalidated_legacy`.
- Weighted marker completeness (`weighted_order_completeness` in the per-query
  `.summary.tab`) averages repeated marker rows when it falls back to all
  orders, so a duplicated row no longer weighs twice. These values can differ
  from 2.0.3; the headline completeness columns are not computed from them.
- Database setup uses one source: the environment override, configured archive,
  or built-in default. It does not replace that source with a moving Zenodo
  release or fall back to an older resource schema. An automatically discovered
  Zenodo archive must have a valid SHA-256 sidecar.
- Compact database validation accepts marker proteins from
  `parquet/faa.parquet` or nonempty legacy `database/faa/*.faa` files.
- The Apptainer and Docker definitions copy the runtime source and bundled
  example without development tests, host wrappers, images, or build files.
  Both definitions use the verified Pixi 0.72.0 base-image digest. Docker
  Compose builds from the repository root.

### Fixed

- The `taxonomy` column of the per-query `.tree_nn` file took the order, family
  and genus of 7-rank database labels from the ranks above them (class, order
  and family). It now reads the right ranks.
- Protein identifiers such as `q|YP_009703.1` keep their whole accession in
  contig attribution. Before, every such protein of a `.faa` input fell into
  one contig called `YP`.
- GC content counts lowercase bases. Query files with the extension `.fa` are
  accepted. The per-query archive no longer holds `<query>_reformatted.fna`,
  which duplicated `query_fna/<query>.fna`.
- Required search, tree and summary failures propagate to query failure.
  Retried queries retain prior work in a separate directory and use only
  current-attempt evidence. Sequence and GFF files are available
  before quality scoring. GFF phases match the gene caller's output.
- Markers with two unequal-length sequences skip tree building before
  alignment validation. Automatic worker allocation respects
  requested thread budgets below four.
- Nearest-neighbor search in marker trees and species-tree representative
  selection no longer treat a reference leaf whose name starts with the query ID
  (for example `G10|c_1` for query `G1`) as a query leaf.
- A contig is counted once in the contamination features when its reformatted
  nucleotide record (`<parent>__<contig>|<contig>`) and its protein records name
  it differently. A single 24 kb contig with five cellular-marker proteins used to
  count as two contigs and reported 5.9% suspicious bases instead of 100%.
- Marker FASTA files list their hit sequences in sorted order. Temporary files of
  marker trimming and of species-tree reference subsets are removed when a step
  fails, and subset files carry a process and thread suffix.
- Materialized reference files are cached per resource identity, so a replaced
  resource never reads a cache made from an older one, and two resources can share
  one cache directory.
- A truncated or malformed row in a multi-HMM search output fails the query. The
  deduplication step used to drop such rows and let the run finish with less
  evidence.
- The database check no longer asks Zenodo for the latest release when the
  installed database is complete.
- A query whose working directory would contain its own input file stops with an
  error instead of moving or deleting the input.
- Cached container images from different HTTPS addresses no longer share a file
  name. `--sensitive` with `--no-sensitive` resolves to non-sensitive in container
  mode, as in native mode. `--apptainer .` names its output from the resolved
  directory. An unrelated `src/` directory next to the launcher no longer counts
  as a GVClass checkout. `--threads-per-worker` is capped at `--threads`.

### Performance

- Each worker thread reuses one summarizer, and contamination scoring parses the
  BLAST output once per query instead of twice. No benchmark figure is claimed.

### Release checks

- Container publication waits for lint, tests, repository hygiene, and the
  resource-backed golden regression. The workflow rejects a tag that differs
  from the software version, verifies the embedded resource and model
  identities, and retains the SIF and its SBOM before the Sylabs push. It
  requires resource v2.1.0.

## [2.0.3] - 2026-07-27

Compatible resource bundle: v2.0.0 (Zenodo DOI 10.5281/zenodo.21225457),
unchanged from 2.0.2.

### Fixed

- Fixed `No genetic codes produced valid results` during genetic-code selection
  when PyHMMer could not classify the first predicted protein or read it as DNA
  (#18). GVClass now opens predicted protein files with an explicit amino
  alphabet.

## [2.0.2] - 2026-07-24

Compatible resource bundle: v2.0.0 (Zenodo DOI 10.5281/zenodo.21225457), unchanged
from 2.0.1.

### Fixed

- Runs are reproducible. The same input now produces the same output every time.
  BLAST results were written in the order pyswrd's parallel search happened to
  return them, and equal-scoring hits were left in arrival order, so two runs of
  one input could select different reference sequences, build different gene
  trees, and report different per-rank taxonomy votes. Two runs of the bundled
  example differed in the `genus` and `species` columns. Hits are now written
  ordered by query, then score, e-value and target.

  This changes some reported values, because the previous ones depended on which
  of several equally-scoring references happened to be chosen. On the bundled
  example the assigned taxonomy, completeness and contamination are unchanged;
  what moves is per-rank vote detail, `avgdist`, and `capsid_group`. Most of the
  movement is toward stronger agreement: PkV-RF01's family and class consensus go
  from 95.24 percent to 100 percent, and the average nearest-reference distance
  for GVMAG-S-1096109-37 drops from 0.20 to 0.09.
- `.fasta` and `.fas` files are no longer misread as DNA when their sequences
  are protein. File type decided membership in the IUPAC nucleotide alphabet,
  which contains the ambiguity codes `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`,
  `H`, `V`, so a peptide spelled only from those letters was classified as DNA,
  validated against the DNA alphabet, and sent to nucleotide gene calling. Type
  is now decided by alphabet, base diversity, and base composition together,
  which also covers low-complexity protein repeats such as the silk-fibroin
  `GAGAGS` motif. Nucleotide contigs carrying ambiguity codes or assembly-gap
  `N` runs are unaffected.
- Input files between 50 MB and 500 MB are accepted. The limit was 50 MB, set on
  a constant named for sequence length in bp but compared against the file size
  in bytes, which rejected the assembly FASTA files that `--contigs` mode exists
  to process. The new ceiling is bounded by memory: validation and gene calling
  each hold a whole file in memory.

### Changed

- The built-in default for `pipeline.threads` is `4`, matching the shipped
  configuration file and the documentation. It was `4` in the file and `16` in
  the code, so a run without a configuration file used `16`.
- `config/gvclass_config.yaml` is the single source of default settings. The
  fallback in the CLI, used by installs without `config/`, is pinned to the
  shipped file by a test.
- Genetic-code selection is no longer a configuration key. GVClass always
  evaluates the full panel and keeps the best-scoring code per query, which is
  what it already did; the `genetic_codes` keys were never read.
- Removed configuration keys that no code read: `quality.min_markers`,
  `quality.recommended_length`, `database.expected_size`,
  `resources.memory_limit`, `resources.task_timeout`, `logging.level`,
  `logging.keep_temp`. Log verbosity follows `-v` / `--verbose`.

### Internal

- `pyrodigal` is pinned to a commit, matching how `witchi` is pinned, so the
  gene caller cannot move when the lock file is refreshed.
- Deleted unreachable validation and error-handling helpers, and moved benchmark
  and prototype scripts out of `src/` so they are no longer packaged.
- `FullSummarizer` is split from 1778 lines into 420, with the extracted logic
  under `src/core/summarize/`. Contamination scoring deliberately stays on the
  class: it is the only step that both reads and writes the result row.
- One worker-planning heuristic replaces the two that had diverged between the
  CLI and the runner.

## [2.0.1] - 2026-07-10

Compatible resource bundle: v2.0.0 (Zenodo DOI 10.5281/zenodo.21225457).

### Fixed

- Runs with one or more failed queries now exit with an error after preserving
  the partial summary and `run_status.json`. They no longer print a successful
  completion message or return status 0.
- Query files with the same stem, such as `sample.fna` and `sample.faa`, are
  rejected before output or database setup because they would overwrite the
  same per-query files.
- Concurrent Parquet label and marker materialization now uses a distinct
  temporary file for each worker.
- If both a database swap and its rollback fail, the previous database remains
  available as `<database>.old` for manual recovery.

### Release checks

- The compact Parquet resource audit now checks all label and protein records,
  and the audit runs alongside the example-genome golden regression on every CI
  event.
- The example-genome fixture now records the v2.0.0 database's current capsid
  group result.
- Container definitions and wrappers now identify the 2.0.1 software image.
  The embedded database and persistent resource-cache namespace remain v2.0.0.

## [2.0.0] - 2026-07-06

Compatible resource bundle: v2.0.0 (Zenodo DOI 10.5281/zenodo.21225457).

Capscan / Preplasmiviricota augmentation with the v2.0.0 resource bundle.

### Added — configurable minimum-length gates
- `--min-length` and `quality.min_length` set the total nucleotide length floor
  for bin/MAG `.fna` inputs. `--min-length 0` removes that length floor while
  keeping normal FASTA validation.
- `--contigs-min-length` overrides `pipeline.contigs_min_length` for `--contigs`
  runs. The same contig-mode floor is used when validating the split contig
  files, so retained contigs are governed by one threshold.

### Changed — cleaner output directory layout
- The CLI now renders one event-driven progress bar across all queries instead
  of polling loose summary files.
- New runs write one `run_status.json` manifest plus `run.log` for run metadata
  and resume state, replacing new per-query `.SUCCESS` files.
- Per-query `summary.tab`, `final_summary.tsv`, and contamination-candidate
  diagnostics are kept inside each `<query>.tar.gz` instead of being duplicated
  at the output-directory root.
- The extended combined diagnostics are archived as
  `gvclass_summary.extended.tar.gz` containing `gvclass_summary.extended.tsv`
  and `.csv`.

### Changed — standardized marker-panel columns in the final summary table
- Marker panels are now reported uniformly as `{panel}_completeness` (`n/N`:
  distinct models present out of the panel size) plus `{panel}_dup` (duplication
  factor = total marker hits / distinct markers present), matching the existing
  `vp_completeness`/`mirus_completeness` convention.
- `gvog4_unique` → `gvog4_completeness` (n/4) + `gvog4_dup`;
  `gvog8_unique`/`gvog8_total` → `gvog8_completeness` (n/8) (`gvog8_dup` kept);
  `mrya_unique`/`mrya_total` → `mrya_completeness` (n/6) + `mrya_dup`;
  `phage_unique`/`phage_total` → `phage_completeness` (n/20) + `phage_dup`.
- The combined `cellular_unique`/`cellular_total` columns are split into
  `busco_completeness` (n/255) + `busco_dup` and `cog_completeness` (n/56) +
  `cog_dup`.
- Removed `mcp_total` (all-MCP count). `ncldv_mcp_total` is retained, and
  `capsid_group` continues to report capsid typing.
- `vp_completeness`/`mirus_completeness` denominators are now derived from the
  number of core categories (still `n/4`) so they cannot go stale if the
  category definitions change.
- All four `species_tree_*` columns report `nd` (not determined) when no species
  tree was built, replacing the empty cells and the prior
  `no-species-tree-calculated` sentinel on `species_tree_nn_taxonomy`.
- Internal cellular/phage marker features consumed by the trained contamination
  model are unchanged; only the summary-table presentation changed.

### Added — IQ-TREE as an opt-in tree builder
- **VeryFastTree remains the default** (`tree_method: veryfasttree`). IQ-TREE
  (`--fast`, model `Q.pfam+R10+F`) is now selectable via `--tree-method iqtree`
  (or `pipeline.tree_method: iqtree`) for more accurate placement at a much higher
  runtime cost. IQ-TREE is pinned to 3.1.2 (latest on bioconda; 3.1.3 is GitHub-only).
  `fasttree` is accepted as an alias for `veryfasttree`.
- **New `--iqtree-mode` / `pipeline.iqtree_mode`** for the species tree:
  `fast` (default, `--fast`) or `ufboot` (ultrafast bootstrap `-B 1000 -bnni`, which
  adds branch-support values and writes a `.contree` consensus tree). UFBoot is much
  slower, so it is scoped to the species tree; per-marker gene trees always use
  `--fast`.

### Fixed — database auto-update on a configured version bump
- gvclass now re-downloads the database when the installed `DB_VERSION` is older
  than the version pinned in `config/gvclass_config.yaml`, even when the configured
  source is not Zenodo (e.g. the gvclass-dev tunnel bundle). Previously a config
  bump (v1.7.0 -> v1.7.1) on a non-Zenodo source left testers on the stale bundle.
  Interactive runs prompt; non-interactive runs update automatically.

### Changed — final summary schema overhaul (breaking)

The `gvclass_summary.tsv` columns changed:
- Removed `taxonomy_strict`. Promoted `species_tree_nn_taxonomy` next to
  `taxonomy_majority` (all four `species_tree_*` columns are `nd` without
  `--species-tree`).
- Renamed `estimated_completeness_quality` → `completeness_model_reliability`;
  dropped `estimated_completeness_advisory` and `estimated_completeness_r2_holdout`.
- Replaced the Bellas-only `capscan_group` with a unified `capsid_group`
  (`label:count` across the NCLDV/Mirus phyla + caps groups, e.g.
  `Nucleocytoviricota:4,Gossevirus:1`).
- Retired `vp_df` / `mirus_df`; dropped `cellular_dup`. (The marker-panel columns
  were then standardized — see the marker-panel entry above.)
- Moved the per-contig contamination diagnostics (`cellular_coherent_*`,
  `cellular_lineage_purity_median`, `cellular_hit_identity_median`,
  `viral_bearing_contig_count`, `contig_attribution_mode`) to a new always-on
  supplementary table `gvclass_summary.extended.tsv`.

Stale per-query `.final_summary.tsv` resume files lose renamed/removed columns
gracefully on re-read (name-keyed; no crash).

### Fixed — spurious single-gene PPV hits on NCLDV genomes
- Taxonomy now casts **one vote per physical query protein** (its nearest tree
  placement across markers), so a shared core gene — e.g. the A32 packaging
  ATPase, marker `PLV_PC_054` — can no longer cast multiple/cross-domain votes.
- Re-curated the PPV reference set: 254 protist-derived proteins (mostly
  NCLDV-like A32 carved out of EUK-co-labeled P10K/GCA/EP assemblies and
  mislabeled PPV) were reassigned to their true domain by blast against
  non-protist viral references. Ships with the new resource bundle.

### Added
- **`--species-tree`: opt-in NCLDV GVOG8 supermatrix species tree + tree-placement
  taxonomy.** For genomes gvclass calls NCLDV, the run gathers each query's nearest
  NCLDV references per GVOG8 marker (dedicated NCLDV-only search → per-group gene
  tree), concatenates one representative protein per marker into a supermatrix
  (pyfamsa alignment → `witchi` chi-squared column pruning, with pytrimal/unpruned
  fallbacks → VeryFastTree), and assigns each query a **root-invariant**
  nearest-reference lineage. By default it builds **one species tree per query**
  (`NEIGHBORS_PER_QUERY_TREE=30` reference genomes/marker), writing
  `out/species_tree/<query>/` (`<query>.treefile`, `<query>.partitions.txt`,
  `species_tree_taxonomy.tsv`) and adding `species_tree_nn_taxonomy`,
  `species_tree_nn_genome`, `species_tree_nn_distance`, and `species_tree_clade_id`
  to `gvclass_summary.tsv` (the per-query placement owns these columns).
  `--species-tree-combined` additionally builds one combined tree over all queries
  (`NEIGHBORS_PER_COMBINED_TREE=20`), writing `out/species_tree/combined.*` (or
  `out/species_tree/_combined/<panel>/` for multi-domain batches) as an extra
  artifact. Both neighbor counts are adjustable in
  `src/core/species_tree/config.py`. The supermatrix column trimmer is selectable via
  `--species-tree-trim {witchi,pytrimal,none}` (default `witchi`, sped up to 20
  permutations / 10 columns per pruning step). The route is fully isolated from
  standard per-marker classification — a run without the flag, or a non-NCLDV genome,
  is byte-for-byte unchanged. PPV and MIRUS domains are registered as additional
  panels. Adds the pure-python `witchi` dependency.

### Changed
- **Preplasmiviricota domain renamed `PLV`/`VP` → `PPV`.** The unified
  Preplasmiviricota domain (Polinton-like viruses + virophages, merged in the
  v1.7.0 bundle) now reports as `PPV` in the `domain` column and lineage strings,
  replacing the class-level misnomer `PLV`. Marker names (`PLV_MCP_*`, `VP_MCP_*`,
  `plv_mcp_caps_*`, …) and the `plv` / `vp_*` summary columns are unchanged. This is
  a label-only rename (tree topology and vote counts are identical); it requires the
  matching `PPV` resource bundle.
- **Same-family marker models are consolidated into grouped trees.** Models that
  share reference proteins (e.g. the MCP / PolB / ATPase families) are merged into a
  single alignment → tree → nearest-neighbour vote, removing the per-model
  vote over-counting that previously fragmented classification. Raw per-marker count
  columns (`*_total`, `*_unique`) are unaffected; contamination scoring is group-aware.

### Added
- **Bellas & Sommaruga capscan MCP markers are now first-class tree markers.** The
  caps major-capsid-protein HMMs feed the `mcp_plv` / `mcp_ncldv` group trees, with
  external (figshare `annotated_PLVs`) reference proteins labelled by Bellas group so
  nearest-neighbour placement can resolve Polinton-like / virophage / NCV-like capsid
  groups. References are curated anti-circularly (external, not GVClass's own homologs);
  the PgVV group is held out for validation (held-out PgVV proteins place to the PgVV
  group 43/43).
- **`capsid_group` summary column** (superseding the earlier Bellas-only `capscan_group`;
  see *Changed* above) reports a unified capsid-type tally from the MCP-tree
  nearest-neighbour placements — caps families plus the Nucleocytoviricota/Mirusviricota
  phyla — as `label:count` (e.g. `Nucleocytoviricota:4,Gossevirus:1`).

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
