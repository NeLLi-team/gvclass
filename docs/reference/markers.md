# Marker panels and genetic codes

Every query is scored against a fixed set of HMM marker panels. Each panel contributes a completeness column and, where applicable, a duplication column to the summary table. Gene calling for nucleotide input is run under nine genetic codes, and the selected code is recorded per query.

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

For each panel, `{panel}_completeness` reports distinct marker models present over panel size (for example `8/8`), and `{panel}_dup` reports total hits over distinct models present (a duplication factor).

!!! note
    Two panels deviate from the `{panel}_dup` convention. The virophage panel emits `vp_completeness` and `vp_mcp` (count of VP MCP hits), not `vp_dup`. The Mirusviricota panel emits `mirus_completeness` only. `capsid_group` is a `label:count` tally across the Nucleocytoviricota and Mirusviricota phyla and the Bellas & Sommaruga capsid groups; `ncldv_mcp_total` is the NCLDV-specific MCP count; `plv` counts A32 proteins (`PLV_PC_054`) that place with PPV references and flags Polinton-like viruses and virophages within the PPV (Preplasmiviricota) domain. It is `0` for ordinary NCLDV.

Each marker model carries a majority functional annotation derived from its member proteins. The full per-model table is in [Marker annotations](marker-annotations.md).

Order-level markers are a separate panel of 576 order-conserved orthologous groups, built only when fast mode is off (`-e`/`--extended`). The default fast mode skips them. See [Tune speed and accuracy](../how-to/tune-speed-and-accuracy.md) for the speed and resolution trade-off.

Per-column definitions for the full summary table are in [Output reference](output.md); interpretation of completeness, contamination, and duplication is in [Quality metrics](../explanation/quality-metrics.md).

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

Every code above is tested for every nucleotide query. This is not configurable; see [Configuration reference](configuration.md).

Selection rule, applied per query:

- Start from meta (code 0).
- Replace it when another code yields at least 2 more complete marker hits.
- Or an average best-hit score at least 10 percent higher.
- Or a coding density at least 5 percent higher.

Any one of the three conditions is sufficient.

The selected code is written to the `ttable` column. When the pyrodigal meta model wins, `ttable` reads `codemeta`. Protein (`.faa`) input skips gene calling and reports `ttable` `no_fna`, and its nucleotide statistics (`GCperc`, `CODINGperc`) are `0.00`.

For the full gene-calling and marker-detection sequence, see [How it works](../explanation/how-it-works.md).
