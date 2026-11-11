# Output files and columns

Results are written to `-o OUTPUT_DIR`, or `<name>_results` in the current directory by default (`<name>` is the input directory name, or the file name without its suffix). The main result is `gvclass_summary.tsv`, with one row per query.

## Files

### Run-level files

| File | Contents |
| --- | --- |
| `gvclass_summary.tsv` | Main summary table, one row per query, tab-separated. |
| `gvclass_summary.csv` | Same content as `gvclass_summary.tsv`, comma-separated. |
| `gvclass_summary.extended.tar.gz` | Archive containing `gvclass_summary.extended.tsv` and `gvclass_summary.extended.csv`, the per-contig contamination diagnostics. |
| `gvclass_failed_queries.tsv` | Written only when one or more queries fail; one row per failed query with its error, including input files rejected by validation (status `invalid_input`). Tab-separated. A failed required search, tree, or summary step fails the query; GVClass does not fall back silently. |
| `gvclass_failed_queries.csv` | Same content as the failed-queries TSV, comma-separated. |
| `results.jsonl` | The record of finished queries. It is append-only: one JSON line per query that finished, written and flushed to disk when the query ends. The last line of a query counts. |
| `run_status.json` | Description of the run, written when the run starts (`running`) and when it ends (`complete` or `failed`): software and database versions, database path, all settings, and each query's status, finish time and error. |
| `run.log` | Run and query start, completion, and failure messages. |
| `run_config.json` | The resolved settings of the run, which the command writes and the runner reads. |
| `gvclass.log` | Log of the runner process. `--verbose` adds debug messages. |

A line of `results.jsonl` holds `query`, `status` (`complete`, `failed` or
`invalid_input`), `error`, `input_name` (the input file name with its suffix),
`input_sha256`, `settings` (the settings that change
results), `database_version`, `software_version`, `finished_at`, `archive` (file name,
size and SHA-256) and `summary` (the final-schema row of the query, including the
extended columns). GVClass computes `input_sha256` before it processes the query. If the
input file has another SHA-256 when the query ends, the query is recorded as `failed`.
Output validation reads this file for the queries of the current run (the queries in
`run_status.json`); a record of an input that this run no longer has does not count.
It reports a query whose last record is not `complete`, an archive that is missing or
has another size than recorded, and an archive whose summary or final summary lacks
a data row for the query.

`--resume` reads `results.jsonl` once and skips a query only when its last record is
`complete`, and all of these equal the current run: the name and the SHA-256 of the
input file (a rename from `.faa` to `.fna` runs the query again), the
result settings (`tree_method`, `iqtree_mode`, `mode_fast`, `completeness_mode`,
`sensitive_mode`, `species_tree`, `species_tree_trim`), the database version and the
software version. The archive must exist with the recorded size. The log names the
field that differed for every query that runs again. Thread counts, `--screen` and
`--contigs` do not change results and do not stop a skip. The combined summary tables
are built from the queries processed now and from the `summary` of the skipped ones.
The screen also covers the skipped queries and sets the five contamination columns of the
combined tables from its result, so a run with the screen and a run without it differ in
those columns.

### Per-query files

Each `<query>.tar.gz` archive contains the summary and diagnostic files below.

| File | Contents |
| --- | --- |
| `<query>.tar.gz` | Bundled per-query artifacts. |
| `<query>/<query>.final_summary.tsv` | Single-query row with the main columns plus extended quality fields. The same row is the `summary` of the query's `results.jsonl` record; `--resume` rebuilds both run summaries from that field. |
| `<query>/<query>.summary.tab` | Per-query summary table, tab-delimited. |
| `<query>/stats/<query>.contamination_candidates.tsv` | Suspicious contigs; written when contamination reaches the reporting threshold (at least 10), its type is interpretable, and candidates are found. |
| `<query>/stats/<query>.contamination_review.tsv` | Contigs with mixed viral and coherent cellular evidence; written whenever this evidence exists, even for a low or unavailable estimate. |

The per-query files are written before the screen runs. Their contamination columns hold the marker-model values, also for a bin whose row in `gvclass_summary.tsv` carries the screen estimate. This applies to `<query>.final_summary.tsv`, `<query>.summary.tab`, the two contamination tables in `stats/`, and the `summary` field of `results.jsonl`.

A running query keeps its working files in `.gvclass_work/<query>/`, never in a folder
named like the query, so a query name cannot collide with `contig_screen/`,
`species_tree/` or the folder of another query. A query that fails leaves its folder
there; `.gvclass_work/` itself goes when no query left a folder. A new attempt deletes a leftover folder before it starts, so old search results
and sequence files cannot affect its scores. A new `<query>.tar.gz` replaces the
previous one only when it is complete; a failed attempt leaves the previous archive in
place.

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

| Column | Description |
| --- | --- |
| `query` | Input filename for the query. |
| `taxonomy_majority` | Accepted lineage from per-marker nearest-neighbor votes. Tied and unsupported ranks are blank, with all lower ranks. The domain tie exception selects exactly one viral domain tied with one or more cellular domains if it meets the support floor; lower ranks are judged as usual within that domain. A protein whose nearest references tie votes only at the ranks those references share. |
| `species_tree_nn_taxonomy` | Nearest-reference lineage from the query's concatenated-marker tree; `nd` if unavailable. |
| `taxonomy_confidence` | `high`, or comma-separated flags in this order: `no_support` for absent taxonomy evidence; `low_support` for insufficient marker support; `tie:<rank>` when two or more taxa share the most marker votes at that rank, retained when a viral domain wins the domain tie exception; `mixed_cellular` when eukaryotic markers win the domain vote without 80% of the cellular-domain markers (viral winners are exempt); `plurality:<rank>` when an assigned rank is held by half or fewer of the markers, including a viral domain that wins a tie; `reduced_fastmode` for an order accepted with two markers in fast mode. `<rank>` is the shallowest rank the flag applies to, for example `plurality:genus`. See [Taxonomy](../explanation/taxonomy.md#confidence-flags). |
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
| `estimated_contamination` | Contamination estimate (%). A bin that the cellular screen scored takes `screen_cellular_pct`. Every other bin takes the prediction of the marker model (`extra_trees_v1`), and unsupported inputs or missing evidence produce `NaN`. |
| `contamination_query_support` | Query evidence status. `supported` for a screen estimate. For the marker model, supported evidence does not establish model calibration. |
| `estimated_contamination_strategy` | `contig_screen_cellular_v1` when the cellular screen scored the bin. Otherwise `extra_trees_v1` or `unavailable`. |
| `contamination_model_reliability` | `simulated_benchmark` for a screen estimate: its error was measured on simulated bins. `unvalidated_legacy` for the marker model and its corrected feature definitions. |
| `contamination_type` | Screen estimate: `cellular` if the cellular share is at least 10, `mixed_viral` if it is below 10 and `screen_giant_virus_pct` is at least 10, otherwise `clean`. Marker model: `clean` below the reporting threshold (at least 10); otherwise `cellular`, `mixed_viral`, `phage`, `duplication`, or `uncertain`. `clean` is a threshold result, not proof of absence. Marker-model estimates that are missing, or that rest on mixed contigs alone, are `uncertain`. |
| `gvog4_completeness` | Distinct core NCLDV GVOG4 markers present, as `n/4`. |
| `gvog4_dup` | GVOG4 duplication factor (total marker hits / distinct markers present). |
| `gvog8_completeness` | Distinct core NCLDV GVOG8 markers present, as `n/8`. |
| `gvog8_dup` | GVOG8 duplication factor. |
| `busco_completeness` | Eukaryotic BUSCO markers present, as `n/69`. |
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
| `ttable` | Genetic code used for gene calling; `codemeta` for metagenomic mode (code 0); `no_fna` for protein inputs. |
| `species_tree_nn_genome` | Nearest-reference genome identifier from the query's species tree; `nd` if unavailable. |
| `species_tree_nn_distance` | Tree distance to that reference, rounded to two decimal places; `nd` if unavailable. |
| `species_tree_clade_id` | Query-only clade identifier; `nd` when no such clade is present. Combined-tree clades are reported in their separate taxonomy table. |

## Screen columns

The last four columns of `gvclass_summary.tsv` report the contig screen. They are
percentages of the input base pairs with two decimals, and `NaN` when the screen did
not produce the number.

| Column | Meaning |
| --- | --- |
| `screen_status` | `ok`: both screens scored the bin. `partial`: a screen could not score the bin (`not_scored`), or one screen failed while the other scored. `failed`: a screen failed and none scored. `disabled`: `--no-screen` or `pipeline.contig_screen: false`. `not_installed`: the add-on is missing. `not_applicable`: protein input or `--contigs` mode. |
| `screen_cellular_pct` | Percent of base pairs called by the cellular screen. |
| `screen_giant_virus_pct` | Percent of base pairs in contigs of 5 kb or more that the giant-virus screen calls foreign. |
| `screen_flagged_pct` | Percent of base pairs called by either screen. `NaN` unless both screens scored the bin. |

`screen_cellular_pct` is also the contamination estimate. A bin with `screen_status` `ok`, or `partial` with a number in `screen_cellular_pct`, takes it as `estimated_contamination` and also sets `contamination_query_support`, `estimated_contamination_strategy`, `contamination_model_reliability` and `contamination_type`. For every other status, these five columns keep the marker-model values. A `NaN` share is never read as 0.

`screen_giant_virus_pct` and `screen_flagged_pct` are review flags, not estimates. See [How the contig screen works](../explanation/contig-screen.md#what-the-summary-columns-report).
A failed query has no row in the summary, so it has no screen columns.

## Contig screen files

A classification run writes the screen files of each bin to
`<output>/contig_screen/<bin>/`. `gvclass screen` writes them to its own output
directory. The folder is outside the query folders, so the per-query archives do not
hold it.

| File | Columns |
| --- | --- |
| `<bin>.cellular.tsv` | `contig_id`, `source_contig`, `source_start`, `length_bp`, `s_n_genes`, `c_gc`, `s_coding_density`, `p_cellular` (stage 1), `p_gbm`, `set_logit`, `p_context` (stage 2), `called_cellular` |
| `<bin>.giant_virus.tsv` | `contig`, `length_bp`, `n_markers`, `n_dup_sc`, `dup_weighted`, `pair_shared`, `d_bulk`, `p_foreign`, `called_foreign` |
| `<bin>.summary.txt` | The printed summaries of both screens. |
| `screen.log` | Output of the screen process, written in a classification run only. The warning on a failed screen shows the end of it. |
| `contig_screen_summary.tsv` | `bin`, `cellular_status`, `giant_virus_status` (`ok`, `failed` or `not_scored`), `input_bp`, `cellular_bp`, `cellular_pct`, `giant_virus_bp`, `giant_virus_pct`, `flagged_bp`, `flagged_pct`, `manifest_sha256` (SHA-256 of the add-on manifest that the screens used), `input_sha256` (SHA-256 of the bin file). A number of a screen that did not score the bin is `NaN`. |

`called_cellular` and `called_foreign` mark windows or contigs at or above the
thresholds stored in the add-on. A call is a review flag.

## Extended quality fields

The table inside `gvclass_summary.extended.tar.gz` contains the query ID and these diagnostic groups:

| Fields | Meaning |
| --- | --- |
| `cellular_coherent_contig_count`, `cellular_coherent_protein_fraction`, `cellular_coherent_bp_fraction` | Coherent cellular evidence, including mixed contigs. Fractions are percentages. These are evidence fractions, not validated foreign-DNA estimates. |
| `cellular_lineage_purity_median`, `cellular_hit_identity_median` | Agreement among cellular lineage hits and their sequence identity. |
| `viral_bearing_contig_count`, `mixed_evidence_contig_count`, `mixed_evidence_status` | Viral and mixed cellular/viral contig evidence. Mixed contigs require review. |
| `ambiguous_eve_contig_count`, `ambiguous_eve_protein_count`, `ambiguous_eve_evidence_status` | Putative endogenous viral element evidence, kept separate from ordinary cellular evidence. |
| `contig_length_basis`, `protein_estimated_contig_count`, `contig_attribution_mode` | Whether lengths come from DNA or protein estimates, and whether protein headers support contig assignment. |
| `contamination_model_applicability`, `contamination_feature_version`, `contamination_feature_contract` | Target and search applicability, and the runtime model-input definition of the marker model. |
| `contamination_prediction_range_status`, `contamination_model_prediction_upper_bound` | Range status and fitted upper bound of the marker model. The bundled ensemble has an upper bound of `58.200481`. |

The contamination fields in this table describe the marker model, also for a bin whose `estimated_contamination` comes from the screen.

## Taxonomy label namespaces

Taxonomy labels come from the installed database. Putative endogenous viral element references use `d_EUK-pEVE` and lower-rank names ending in `-pEVE`, such as `p_Discosea-pEVE` or `g_Vannella-pEVE`.

`EUK-pEVE` identifies the eukaryotic source lineage and is distinct from `EUK`, NCLDV, PPV, and MIRUS assignments. These references can enter species trees if they meet the panel's marker threshold; their taxonomy remains `EUK-pEVE`. See [taxonomy and classification](../explanation/taxonomy.md).

See [quality metrics](../explanation/quality-metrics.md) for interpretation and [markers](markers.md) for panel definitions.
