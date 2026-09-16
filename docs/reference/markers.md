# Marker panels and genetic codes

Marker panels provide presence, duplication, and capsid-type counts in the [summary table](output.md). Nucleotide inputs are evaluated under nine genetic codes.

## Marker panels

| Panel | Column prefix | Size | Purpose |
| --- | --- | --- | --- |
| GVOG4 | `gvog4` | n/4 | Core NCLDV single-copy orthologs |
| GVOG8 | `gvog8` | n/8 | Core NCLDV single-copy orthologs |
| BUSCO eukaryotic | `busco` | n/255 | Eukaryotic carry-over flag |
| Universal COG (UNI56) | `cog` | n/56 | Universal cellular carry-over flag |
| Mryavirus | `mrya` | n/6 | Mryavirus markers |
| Phage (geNomad) | `phage` | n/20 | Phage contamination flag |
| Virophage core | `vp` | n/4 | MCP, Penton, ATPase, Protease |
| Mirusviricota core | `mirus` | n/4 | MCP, ATPase, Portal, Triplex |
| Capsid typing | `capsid_group`, `ncldv_mcp_total` | count | Capsid (MCP) type tally |
| PPV flag | `plv` | count | A32 (`PLV_PC_054`) proteins placing with PPV references |

`{panel}_completeness` reports markers present over panel size, such as `8/8`. For panels with `{panel}_dup`, duplication is total hits divided by distinct markers present.

Virophage and Mirusviricota completeness counts marker categories. The virophage panel also reports `vp_mcp`; neither panel reports a duplication column. `capsid_group` reports `label:count` values for Nucleocytoviricota, Mirusviricota, and Bellas & Sommaruga capsid groups. The PPV (Preplasmiviricota) group includes Polinton-like viruses and virophages; `plv` is `0` for ordinary NCLDV.

Functional annotations are listed in [marker annotations](marker-annotations.md).

The 576 order-level marker groups are searched in both modes. Their trees are built only with `-e`/`--extended`. See [speed and tree settings](../how-to/tune-speed-and-accuracy.md).

See [quality metrics](../explanation/quality-metrics.md) for interpretation.

## Genetic codes

Nine genetic codes are tested during gene calling:

| Code | Translation table |
| --- | --- |
| 0 | Pyrodigal meta mode (pretrained models) |
| 1 | NCBI standard |
| 4 | NCBI mold, protozoan, coelenterate mitochondrial; Mycoplasma, Spiroplasma |
| 6 | NCBI ciliate, dasycladacean, hexamita nuclear |
| 11 | NCBI bacterial, archaeal, plant plastid |
| 15 | NCBI Blepharisma nuclear |
| 29 | NCBI Mesodinium nuclear |
| 106 | Genetic code similar to code 6, found in some novel giant virus genomes |
| 129 | Genetic code similar to code 29, found in some novel giant virus genomes |

Codes are ranked by complete marker hits, average best-hit score, coding density, then code preference. If the top-ranked candidate is not code `0` and code `0` is available, that candidate must exceed code `0` by at least one of these margins:

- 2 complete marker hits;
- 10% in average best-hit score;
- 5% in coding density (relative increase).

Otherwise, code `0` is retained. The code panel is fixed.

The selected code appears in `ttable`; code `0` is reported as `codemeta`. Protein (`.faa`) input skips gene calling, reports `ttable=no_fna`, and has `GCperc=0.00` and `CODINGperc=0.00`.

See [how it works](../explanation/how-it-works.md) for the gene-calling and marker-detection workflow.
