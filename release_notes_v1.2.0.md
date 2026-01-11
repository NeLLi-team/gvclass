# GVClass v1.2.0

## Highlights
- Updated reference database v1.2.0 with refined GA thresholds and model annotations.
- Improved HMM model accuracy through updated gathering thresholds.
- Added functional annotations to HMM models.
- Bug fixes and documentation improvements.

## Database updates
- New database tarball: `resources_v1_2_0.tar.gz`.
- Expanded reference coverage:
  - Nucleocytoviricota and Mirusviricota genomes from Vasquez et al. (2025) bioRxiv.
  - Mirusviricota models and genomes from Medvedeva et al. (2026) Nature Microbiology.
  - Models for Polinton-like viruses (PLV) and virophages (PV) from Roux et al. (2023) Biomolecules.
  - Extended VP, PLV and phage reference set from MetaVR database (Fiamenghi et al. 2025).

## New output columns
- **VP (Virophage) metrics**: `vp_completeness`, `vp_mcp`, `vp_df`.
- **PLV metrics**: `plv`.
- **Mirus metrics**: `mirus_completeness`, `mirus_df`.
- **NCLDV MCP**: `ncldv_mcp_total` (includes OG1352, OG484).
- Removed old columns `mirus_unique`, `mirus_total`, `mirus_dup` (replaced by completeness metrics).

## CLI improvements
- `--no-mode-fast` renamed to `--extended` / `-e` (fast mode remains default).

## Bug fixes
- Removed `time.sleep()` usage in `prefect_flow.py`.
- Fixed inconsistent bytes encoding in `genetic_code_optimizer.py`.

## Documentation
- Added full CLI reference table and expanded `--contigs` documentation.
- Documented genetic code selection logic (9 codes tested, selection criteria).
- Mermaid diagram now correctly shows code 0 (meta mode with pretrained models).
