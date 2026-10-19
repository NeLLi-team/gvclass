# Output files and columns

Results are written to `-o OUTPUT_DIR`, or `<query_dir>_results` by default. The main result is `gvclass_summary.tsv`, with one row per query.

## Files

### Run-level files

| File | Contents |
| --- | --- |
| `gvclass_summary.tsv` | Main summary table, one row per query, tab-separated. |
| `gvclass_summary.csv` | Same content as `gvclass_summary.tsv`, comma-separated. |
| `gvclass_summary.extended.tar.gz` | Archive containing `gvclass_summary.extended.tsv` and `gvclass_summary.extended.csv`, the per-contig contamination diagnostics. |
| `gvclass_failed_queries.tsv` | Written only when one or more queries fail; one row per failed query with its error. Tab-separated. A failed required search, tree, or summary step fails the query; GVClass does not fall back silently. |
| `gvclass_failed_queries.csv` | Same content as the failed-queries TSV, comma-separated. |
| `run_status.json` | Resume manifest: software and database versions, database path, settings, query status, timestamps, and artifact names, sizes, and checksums. |
| `run.log` | Run and query start, completion, and failure messages. |

### Per-query files

Each `<query>.tar.gz` archive contains the summary and diagnostic files below.

| File | Contents |
| --- | --- |
| `<query>.tar.gz` | Bundled per-query artifacts. |
| `<query>/<query>.final_summary.tsv` | Single-query row with the main columns plus extended quality fields. `--resume` rebuilds both run summaries from these rows. |
| `<query>/<query>.summary.tab` | Per-query summary table, tab-delimited. |
| `<query>/stats/<query>.contamination_candidates.tsv` | Suspicious contigs; written when contamination reaches the reporting threshold (at least 10), its type is interpretable, and candidates are found. |
| `<query>/stats/<query>.contamination_review.tsv` | Contigs with mixed viral and coherent cellular evidence; written whenever this evidence exists, even for a low or unavailable estimate. |

Both contamination tables include the contig ID, `candidate_type`, `reason`, marker counts, `recommended_action`, `length_bp`, `length_basis`, and `ambiguous_eve_hit_count`. `length_basis` identifies a measured DNA length or a protein-based estimate. A review row is not a deletion instruction.

### Species-tree files

Produced with `--species-tree` or `--species-tree-combined` when eligible queries are present. See [build a species tree](../how-to/build-a-species-tree.md).

| Path | Contents |
| --- | --- |
| `species_tree/<query>/<query>.treefile` | Concatenated-marker species tree for the query. |
| `species_tree/<query>/<query>.supermatrix.faa` | Concatenated protein alignment used for the tree. |
| `species_tree/<query>/<query>.partitions.txt` | Per-marker partition definitions for the query supermatrix. |
| `species_tree/<query>/species_tree_taxonomy.tsv` | Query placement and nearest-reference distances to six decimal places. |
| `species_tree/combined.treefile` | Combined tree of eligible queries and references when eligible queries use one panel. Requires `--species-tree-combined`; alignment and partition files use the same `combined` prefix. |
| `species_tree/species_tree_taxonomy.tsv` | Placements from the single-panel combined tree. |
| `species_tree/_combined/<panel>/` | Combined tree, alignment, partitions, and taxonomy for each panel in batches containing multiple panels. Requires `--species-tree-combined`. |

The four `species_tree_*` summary columns describe per-query trees. Combined-tree placements are recorded separately. `nd` means no value was determined.

## Summary columns

`gvclass_summary.tsv` and `gvclass_summary.csv` share the columns below in this order. `NaN` means an unavailable numeric estimate, not zero.

Resume requires saved quality-support fields. Summaries that lack those fields are reprocessed so that old quality estimates do not enter the current output unchanged.
When a query runs again, any nonempty working directory is retained as a sibling
directory, `<query>.previous.<unique>/`. The new attempt starts in an empty directory so
old search results and sequence files cannot affect its scores.

| Column | Description |
| --- | --- |
| `query` | Input filename for the query. |
| `taxonomy_majority` | Accepted lineage from per-marker nearest-neighbor votes. Unsupported ranks are blank. |
| `species_tree_nn_taxonomy` | Nearest-reference lineage from the query's concatenated-marker tree; `nd` if unavailable. |
| `taxonomy_confidence` | `high`, or flags: `low_support` for insufficient marker support; `mixed_cellular` when eukaryotic markers win the domain vote without 80% of the cellular-domain markers; `reduced_fastmode` for an order accepted with two markers in fast mode; `no_support` for absent taxonomy evidence. |
| `capsid_group` | Combined capsid-type tally as `label:count` across the Nucleocytoviricota, Mirusviricota, and Bellas & Sommaruga capsid groups. |
| `species` | Species-rank call with per-taxon counts. |
| `genus` | Genus-rank call with per-taxon counts. |
| `family` | Family-rank call with per-taxon counts. |
| `order` | Order-rank call with per-taxon counts. |
| `class` | Class-rank call with per-taxon counts. |
| `phylum` | Phylum-rank call with per-taxon counts. |
| `domain` | Domain-rank call with per-taxon counts. |
| `avgdist` | Average tree distance to the reference neighbors. |
| `order_dup` | Average copy number of expected order-level markers; elevated values can indicate duplication, chimerism, or mixed bins. |
| `estimated_completeness` | Marker-recovery estimate (%). It does not certify absolute genome completeness. Read its strategy, support, and reliability. |
| `completeness_query_support` | Query evidence status, separate from model validation. Statuses include `supported`, `limited_support`, `insufficient_evidence`, `unassigned`, `unsupported_target`, and `unavailable`. |
| `estimated_completeness_strategy` | Calculation used: fixed-panel ratio, novelty-aware regression, tier fallback, or `unavailable`. |
| `completeness_model_group` | Resource group used for completeness; `MIRUS` is the Mirusviricota resource group, not an order assignment. `unavailable` means no usable group. |
| `completeness_model_reliability` | Validation status. Models report `unvalidated_legacy`; unavailable estimates can report `unavailable`. Resources validated by marker dropout only can report `advisory_only`. |
| `estimated_contamination` | Selected contamination amount estimate. The method is `extra_trees_v1`. Unsupported inputs or missing evidence produce `NaN`. |
| `contamination_query_support` | Query evidence status. Supported evidence does not establish model calibration. |
| `estimated_contamination_strategy` | `extra_trees_v1` or `unavailable`. |
| `contamination_model_reliability` | `unvalidated_legacy` for the model and its corrected feature definitions. |
| `contamination_type` | `clean` below the reporting threshold (at least 10); otherwise `cellular`, `mixed_viral`, `phage`, `duplication`, or `uncertain`. `clean` is a threshold result, not proof of absence. Missing estimates and evidence limited to mixed contigs are `uncertain`. |
| `gvog4_completeness` | Distinct core NCLDV GVOG4 markers present, as `n/4`. |
| `gvog4_dup` | GVOG4 duplication factor (total marker hits / distinct markers present). |
| `gvog8_completeness` | Distinct core NCLDV GVOG8 markers present, as `n/8`. |
| `gvog8_dup` | GVOG8 duplication factor. |
| `busco_completeness` | Eukaryotic BUSCO markers present, as `n/255`. |
| `busco_dup` | BUSCO duplication factor; elevated values can reflect eukaryotic sequence contamination. |
| `cog_completeness` | Universal COG (UNI56) markers present, as `n/56`. |
| `cog_dup` | COG duplication factor; elevated values can reflect cellular sequence contamination. |
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
| `species_tree_nn_genome` | Nearest-reference genome identifier from the query's species tree; `nd` if unavailable. |
| `species_tree_nn_distance` | Tree distance to that reference, rounded to two decimal places; `nd` if unavailable. |
| `species_tree_clade_id` | Query-only clade identifier; `nd` when no such clade is present. Combined-tree clades are reported in their separate taxonomy table. |

## Contig screen files

`gvclass screen` writes its files to its own output directory.

| File | Columns |
| --- | --- |
| `<bin>.cellular.tsv` | `contig_id`, `source_contig`, `source_start`, `length_bp`, `s_n_genes`, `c_gc`, `s_coding_density`, `p_cellular` (stage 1), `p_gbm`, `set_logit`, `p_context` (stage 2), `called_cellular` |
| `<bin>.giant_virus.tsv` | `contig`, `length_bp`, `n_markers`, `n_dup_sc`, `dup_weighted`, `pair_shared`, `d_bulk`, `p_foreign`, `called_foreign` |
| `<bin>.summary.txt` | The printed summaries of both screens. |
| `contig_screen_summary.tsv` | `bin`, `cellular_status`, `giant_virus_status` (`ok` or `failed`). |

`called_cellular` and `called_foreign` mark windows or contigs at or above the
thresholds stored in the resource archive. A call is a review flag.

## Extended quality fields

The table inside `gvclass_summary.extended.tar.gz` contains the query ID and these diagnostic groups:

| Fields | Meaning |
| --- | --- |
| `cellular_coherent_contig_count`, `cellular_coherent_protein_fraction`, `cellular_coherent_bp_fraction` | Coherent cellular evidence, including mixed contigs. Fractions are percentages. These are evidence fractions, not validated foreign-DNA estimates. |
| `cellular_lineage_purity_median`, `cellular_hit_identity_median` | Agreement among cellular lineage hits and their sequence identity. |
| `viral_bearing_contig_count`, `mixed_evidence_contig_count`, `mixed_evidence_status` | Viral and mixed cellular/viral contig evidence. Mixed contigs require review. |
| `ambiguous_eve_contig_count`, `ambiguous_eve_protein_count`, `ambiguous_eve_evidence_status` | Putative endogenous viral element evidence, kept separate from ordinary cellular evidence. |
| `contig_length_basis`, `protein_estimated_contig_count`, `contig_attribution_mode` | Whether lengths come from DNA or protein estimates, and whether protein headers support contig assignment. |
| `contamination_model_applicability`, `contamination_feature_version`, `contamination_feature_contract` | Target and search applicability, and the runtime model-input definition. |
| `contamination_prediction_range_status`, `contamination_model_prediction_upper_bound` | Model range status and fitted upper bound. The bundled ensemble has an upper bound of `58.200481`. |

## Taxonomy label namespaces

Taxonomy labels come from the installed database. Putative endogenous viral element references use `d_EUK-pEVE` and lower-rank names ending in `-pEVE`, such as `p_Discosea-pEVE` or `g_Vannella-pEVE`.

`EUK-pEVE` identifies the eukaryotic source lineage and is distinct from `EUK`, NCLDV, PPV, and MIRUS assignments. These references can enter species trees if they meet the panel's marker threshold; their taxonomy remains `EUK-pEVE`. See [taxonomy and classification](../explanation/taxonomy.md).

See [quality metrics](../explanation/quality-metrics.md) for interpretation and [markers](markers.md) for panel definitions.
