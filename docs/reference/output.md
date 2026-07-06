# Output files and columns

GVClass writes a run-level summary plus per-query artifacts to the output directory (`<query_dir>_results` by default, or the path passed to `-o`). The flags that add or remove files are listed in the [CLI reference](cli.md).

## Files

### Run-level files

Written to the root of the output directory.

| File | Contents |
| --- | --- |
| `gvclass_summary.tsv` | Main summary table, one row per query, tab-separated. |
| `gvclass_summary.csv` | Same content as `gvclass_summary.tsv`, comma-separated. |
| `gvclass_summary.extended.tar.gz` | Archive containing `gvclass_summary.extended.tsv` and `gvclass_summary.extended.csv`, the per-contig contamination diagnostics. |
| `gvclass_failed_queries.tsv` | Written only when one or more queries fail; one row per failed query with its error. Tab-separated. |
| `gvclass_failed_queries.csv` | Same content as the failed-queries TSV, comma-separated. |
| `run_status.json` | Run manifest used by `--resume`; records software version, database version/path, settings, per-query status, timestamps, artifact names, sizes, and checksums for newly completed queries. |
| `run.log` | Human-readable chronological run log with run start, query start, query completion/failure, and run completion lines. |

### Per-query files

Written for each input query. The top-level file is the archive; the summary and diagnostic files listed below are members inside `<query>.tar.gz`.

| File | Contents |
| --- | --- |
| `<query>.tar.gz` | Bundled per-query artifacts. |
| `<query>/<query>.final_summary.tsv` | Single-query row in the 44-column main schema; `--resume` rebuilds the run summary from these rows. |
| `<query>/<query>.summary.tab` | Per-query summary table, tab-delimited. |
| `<query>/stats/<query>.contamination_candidates.tsv` | Suspicious contigs; written only when `estimated_contamination` is at least 10, `contamination_type` is interpretable, and suspicious contigs are found. |

### Species-tree files

Written only with `--species-tree`. See [Build a species tree](../how-to/build-a-species-tree.md) and the [species tree explanation](../explanation/species-tree.md).

| Path | Contents |
| --- | --- |
| `species_tree/<query>/<query>.treefile` | Concatenated-marker species tree for the query. |
| `species_tree/<query>/<query>.partitions.txt` | Per-marker partition definitions for the query supermatrix. |
| `species_tree/<query>/species_tree_taxonomy.tsv` | Taxonomy and full-precision nearest-reference distances for the query species tree. |
| `species_tree/combined.*` | Combined tree across all queries; written only with `--species-tree-combined`. |
| `species_tree/_combined/<panel>/` | Per-panel combined trees for multi-domain batches; written only with `--species-tree-combined`. |

!!! note
    Without `--species-tree`, no `species_tree/` directory is produced and the four `species_tree_*` summary columns are `nd`.

## Summary columns

`gvclass_summary.tsv` and `gvclass_summary.csv` share the same 44 columns in this order.

| Column | Description |
| --- | --- |
| `query` | Input filename for the query. |
| `taxonomy_majority` | Full lineage from the per-marker single-gene tree nearest-neighbor majority vote. |
| `species_tree_nn_taxonomy` | Lineage of the nearest reference in the concatenated-marker species tree; `nd` unless `--species-tree` was set. |
| `taxonomy_confidence` | `high` when every emitted rank cleared its distinct-marker threshold, otherwise one or more of `low_support`, `reduced_fastmode`, `no_support`. |
| `capsid_group` | Combined capsid-type tally as `label:count` across the Nucleocytoviricota, Mirusviricota, and Bellas & Sommaruga capsid groups. |
| `species` | Species-rank call with per-taxon counts. |
| `genus` | Genus-rank call with per-taxon counts. |
| `family` | Family-rank call with per-taxon counts. |
| `order` | Order-rank call with per-taxon counts. |
| `class` | Class-rank call with per-taxon counts. |
| `phylum` | Phylum-rank call with per-taxon counts. |
| `domain` | Domain-rank call with per-taxon counts. |
| `avgdist` | Average tree distance to the reference neighbors. |
| `order_dup` | Average copy number of expected order-level markers; elevated values indicate duplicated, chimeric, or mixed bins. |
| `estimated_completeness` | Estimated percent of the expected genome recovered for the assigned lineage; the only completeness field in the main table. |
| `completeness_model_reliability` | `advisory_only`, `moderate`, or `high`, from the per-order model hold-out R^2; a property of the model, not the genome. |
| `estimated_contamination` | Output of the trained `extra_trees_v1` contamination model and the primary contamination estimate. |
| `contamination_type` | Likely contamination source (`clean`, `cellular`, `mixed_viral`, `phage`, `duplication`, `uncertain`), populated when `estimated_contamination` is at least 10. |
| `gvog4_completeness` | Distinct core NCLDV GVOG4 markers present, as `n/4`. |
| `gvog4_dup` | GVOG4 duplication factor (total marker hits / distinct markers present). |
| `gvog8_completeness` | Distinct core NCLDV GVOG8 markers present, as `n/8`. |
| `gvog8_dup` | GVOG8 duplication factor. |
| `busco_completeness` | Eukaryotic BUSCO markers present, as `n/255`. |
| `busco_dup` | BUSCO duplication factor; elevated values indicate cellular (eukaryote) carry-over. |
| `cog_completeness` | Universal COG (UNI56) markers present, as `n/56`. |
| `cog_dup` | COG duplication factor; elevated values indicate cellular carry-over. |
| `mrya_completeness` | Mryavirus markers present, as `n/6`. |
| `mrya_dup` | Mryavirus duplication factor. |
| `phage_completeness` | Phage (geNomad) markers present, as `n/20`. |
| `phage_dup` | Phage duplication factor. |
| `ncldv_mcp_total` | Count of NCLDV-specific major capsid protein (MCP) markers. |
| `vp_completeness` | Virophage core markers present (MCP, Penton, ATPase, Protease), as `n/4`. |
| `vp_mcp` | Count of virophage MCP hits. |
| `plv` | Count of A32 (`PLV_PC_054`) proteins placing with PPV references; 0 for ordinary NCLDV. |
| `mirus_completeness` | Mirusviricota core markers present (MCP, ATPase, Portal, Triplex), as `n/4`. |
| `contigs` | Number of contigs in the query. |
| `LENbp` | Total length in base pairs. |
| `GCperc` | GC content as a percentage. |
| `genecount` | Number of predicted genes. |
| `CODINGperc` | Coding density as a percentage. |
| `ttable` | Genetic code used for gene calling; `no_fna` for protein inputs. |
| `species_tree_nn_genome` | Nearest-reference genome id from the species tree; `nd` without `--species-tree`. |
| `species_tree_nn_distance` | Distance to the nearest reference in the species tree; `nd` without `--species-tree`. |
| `species_tree_clade_id` | Clade id of the nearest reference in the species tree; `nd` without `--species-tree`. |

The per-contig diagnostic columns (`cellular_coherent_*`, `cellular_lineage_purity_median`, `viral_bearing_contig_count`, `contig_attribution_mode`) appear only in the extended table inside `gvclass_summary.extended.tar.gz`.

## Taxonomy label namespaces

The selected database bundle defines the taxonomy label namespace that appears in `taxonomy_majority` and the per-rank count columns. Bundles with putative endogenous viral element references can emit `d_EUK-pEVE` at domain rank and lower-rank names ending in `-pEVE`, such as `p_Discosea-pEVE` or `g_Vannella-pEVE`.

`EUK-pEVE` labels mark eukaryotic source-lineage references with pEVE evidence. They remain separate from ordinary `EUK` labels at every rank and are not NCLDV, PPV, or MIRUS assignments. With `--species-tree`, these references can appear as auxiliary nearest-reference leaves inside the NCLDV, PPV, or MIRUS species-tree panels if they pass that panel's marker threshold; the reported taxonomy remains `EUK-pEVE`. See [Taxonomy and classification](../explanation/taxonomy.md#putative-eve-references) for interpretation.

For how to read completeness, contamination, and duplication, see [Quality metrics](../explanation/quality-metrics.md). For the marker panels behind the completeness columns, see [Markers](markers.md).
