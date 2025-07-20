# GVClass Marker Sets Configuration

This directory contains the configuration for marker sets used in the GVClass pipeline.

## Marker Sets

### MIRUS_MODELS
- Mirusviricota-specific markers
- 5 markers total
- Source: doi: 10.1101/2024.01.18.576163

### BUSCO_MODELS
- Benchmarking Universal Single-Copy Orthologs
- 255 markers total
- Source: https://v10-1.orthodb.org

### PHAGE_MODELS
- Phage-specific markers from genomad
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
- 5 markers total
- Key structural proteins for NCLDV identification

### MRYA_MODELS
- Metagenomic Russian Yokohama-Asfarviridae markers
- 6 markers total
- Source: https://doi.org/10.1016/j.virol.2014.06.032

## Usage

To use these marker sets in scripts:

```python
from config import MIRUS_MODELS, BUSCO_MODELS, PHAGE_MODELS
# etc.
```

## Updating Marker Sets

To update marker sets:
1. Edit `marker_sets.py`
2. Update the corresponding list
3. Ensure all markers exist in the database (`resources/database/`)
4. Test the pipeline with updated markers