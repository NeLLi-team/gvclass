# GVClass Marker Sets Configuration

This directory contains the configuration for marker sets used in the GVClass pipeline.

Every model name in a marker set must exist in the reference database. A name that is
absent never matches a hit, so its panel reports a wrong size and wrong completeness.
`tests/test_marker_panels_match_database.py` checks all sets against
`docs/reference/marker-annotations.tsv`.

## Marker Sets

### MIRUS_MODELS
- Mirusviricota-specific markers
- 7 markers total, in 4 categories (MCP, ATPase, Portal, Triplex)
- Source: doi: 10.1101/2024.01.18.576163

### BUSCO_MODELS
- Benchmarking Universal Single-Copy Orthologs
- 69 markers total (the BUSCO models the database ships)
- Source: https://v10-1.orthodb.org

### PHAGE_MODELS
- Phage-specific markers from genomad (`GENOMAD.<number>.VV`)
- 20 markers total
- Used for detecting phage contamination

### GVOG4M_MODELS
- Giant Virus Orthologous Groups - 4 marker set
- Core NCLDV markers for rapid analysis

### GVOG8M_MODELS
- Giant Virus Orthologous Groups - 8 marker set
- Extended NCLDV markers for more detailed analysis

### UNI56_MODELS
- Universal 56 gene set
- 56 markers total
- Source: DOI: 10.7554/eLife.26580

### MCP_MODELS
- Major Capsid Protein markers
- 6 markers total
- Key structural proteins for NCLDV identification

### MRYA_MODELS
- Metagenomic Russian Yokohama-Asfarviridae markers
- 6 markers total
- Source: https://doi.org/10.1016/j.virol.2014.06.032

## Marker groups

`marker_groups.tsv` maps gene-family groups to their member models. The loader raises an
error when the table is missing, has a short row, or holds no group. A group counts as
cellular (or phage) only when all its members are in that set.

## Usage

To use these marker sets in scripts:

```python
from src.config import MIRUS_MODELS, BUSCO_MODELS, PHAGE_MODELS
# etc.
```

`missing_panel_models(available_models)` returns, per panel, the panel models that are not
in a set of database model names.

## Updating Marker Sets

To update marker sets:
1. Edit `marker_sets.py`
2. Update the corresponding list
3. Ensure all markers exist in the database (`resources/database/`)
4. Run `tests/test_marker_panels_match_database.py`, then test the pipeline with the updated markers
