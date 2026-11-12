# Completeness and contamination

GVClass estimates genome quality from marker recovery, marker duplication, taxonomy, reference matches, and contig-level evidence. The main fields are `estimated_completeness`, `completeness_model_reliability`, `estimated_contamination`, and `contamination_type`. These estimates support manual review. They do not prove that a genome is complete or free of contamination.

Read each estimate together with its query-support and model-reliability fields. `NaN` means that an estimate is unavailable. It does not mean zero.

For commands to inspect these fields, see [Assess genome quality](../how-to/assess-genome-quality.md).

## Completeness

`estimated_completeness` is a marker-based estimate on a 0–100 scale. It compares recovered markers with the expected complement for a reference lineage. A value of 100 means that the marker estimate reached its upper limit. It does not demonstrate that the full nucleotide sequence has been recovered.

The default `novelty-aware` mode compares the query's markers with order or family expectations for its accepted `taxonomy_majority`. It uses lineage-specific marker tiers, reference marker counts, and, where available, an order-specific regression model. If a model's recorded hold-out R² is below 0.5, GVClass uses the marker-tier estimate. `--completeness-mode legacy` uses the order-marker recovery ratio scaled to a reference baseline. `estimated_completeness_strategy` records the calculation used. A reference baseline that is not a finite number cannot scale a ratio. The order-marker ratio is then unavailable (`NaN`), not 100. If only the baseline of the marker-tier model is not finite, the order-marker ratio replaces that model and `completeness_query_support` reports `limited_support`.

Mirusviricota uses the `MIRUS` resource group, which is separate from the order in the taxonomy output. `completeness_model_group` records the group used for scoring.

The training target of the bundled regression model is simulated marker retention, not retained genome bases. Missing genomic regions can remove several linked markers or none. Gene calling can also miss partial genes at contig ends.

### Reliability and query support

`completeness_model_reliability` describes the validation evidence for the model:

| Value | Interpretation |
| --- | --- |
| `unvalidated_legacy` | Bundled model. Its historical hold-out scores do not establish reliability. |
| `advisory_only` | Resource metadata marks the model as validated by marker dropout only. |
| `unavailable` | No estimate is available. |

`completeness_query_support` describes the evidence in this query. A model with good validation still needs enough suitable markers to score a new query. An unassigned lineage or a query with no suitable markers cannot support a numeric estimate. The marker tiers of the novelty-aware model can lack informative markers for a query that still has order markers. The novelty-aware estimate is then unavailable. With `--completeness-mode legacy` the order-marker ratio stays available and `completeness_query_support` reports `limited_support`. The `order_completeness_v2*` fields always describe the novelty-aware model. They hold their unavailable values when its resources are missing.

The training of the bundled models used marker dropout, and some reference-derived statistics and repeated subsamples crossed validation boundaries. Their historical hold-out scores do not establish independent nucleotide-completeness accuracy. The runtime marks this calibration as unvalidated. Reference coverage, taxonomic distance, and low marker recovery add query-specific limits.

A low-support adjustment can reduce a prediction but cannot increase one. Neither the adjustment nor clipping to 0–100 creates evidence for a complete genome.

Strategy and support fields are recorded in `<query>/<query>.summary.tab` inside each query's `.tar.gz` archive.

## Contamination

`estimated_contamination` has two sources. When the cellular screen of the [contig screen](contig-screen.md) add-on scored the bin (`screen_status` is `ok`, or `partial` with a number in `screen_cellular_pct`), the estimate is `screen_cellular_pct`, the percent of the input base pairs that the cellular screen calls cellular. The bin then reports `estimated_contamination_strategy` `contig_screen_cellular_v1`, `contamination_query_support` `supported` and `contamination_model_reliability` `simulated_benchmark`. Every other bin takes the prediction of the marker model described [below](#the-marker-model).

On 3,010 simulated bins (bacterial, archaeal, eukaryotic and second-giant-virus additions; 164 new target genomes plus the development panel of the screen), the cellular share had a mean absolute error of 3.1 to 4.4 percentage points in the five groups of more than 100 bins, and 6.6 in one group of 12 bins. The marker model had 14.1 to 16.1 points on the same bins. The cellular share does not count a second giant virus. `screen_giant_virus_pct` and the `mixed_viral` type report it. The `phage` and `duplication` types occur only with the marker model.

For a screen estimate, `contamination_type` is `cellular` when the share is at least 10, `mixed_viral` when the share is below 10 and `screen_giant_virus_pct` is at least 10, and `clean` otherwise. A call is a review flag, not a rule to delete a contig; see [Screen bins for contaminating contigs](../how-to/screen-contigs.md).

The archived per-query summary (`<query>/<query>.final_summary.tsv` inside `<query>.tar.gz`) is written before the screen and keeps the marker-model values. The same holds for `<query>.summary.tab` and the contamination tables in `stats/`.

### The marker model

The marker model is the bundled ExtraTrees regressor (`extra_trees_v1`). The model uses 21 features. These include cellular and phage marker counts, duplicated viral markers, disagreement among viral taxonomic hits, suspicious contig lengths, completeness, and contig-level cellular evidence. `estimated_contamination_strategy` records the model identifier.

Giant viruses can carry genes acquired from cellular organisms. A cellular-like protein or a BUSCO/UNI56 hit alone is insufficient to identify a contaminating contig. Six marker groups hold both cellular and viral models, for example `grp_COG0085`, which holds the cellular RNA polymerase subunit and a giant-virus homolog. The category of a protein in such a group is the category of its best-scoring model in the group, and a tie between categories counts as viral. BUSCO and UNI56 hits contribute to the cellular signal but are not its only source. Reference matches from other marker searches can also contribute. GVClass combines these signals with contig-level lineage evidence, but it cannot guarantee that it separates all horizontally transferred genes from contamination. Inspect the implicated contigs before removing DNA from a bin.

The contig classifier checks whether several cellular matches agree on a lineage and have high sequence identity. A protein whose nearest tree references tie takes the lineage they share, and it has no domain when they differ in domain. A protein whose references have no order stays in the count of cellular proteins but cannot agree with the others. A few cellular-like genes do not establish a coherent cellular contig. A contig with no tree-placed protein cannot be assessed this way. It is judged by the cellular marker and cellular hit rules that also apply when the classifier is inactive. When a contig has both coherent cellular and viral evidence, the classifier keeps both signals. Such a contig can be a chimera, an integrated viral region, or a biological mixture. Review it before sequence removal.

`EUK-pEVE` reference hits identify putative endogenous viral element evidence in eukaryotic source sequences. They are ambiguous for contamination assessment. Do not read them as ordinary cellular contamination or as proof of a free viral genome.

For nucleotide input, contig fractions use the supplied DNA lengths. For protein-only input, length and contig assignment can require estimates from protein headers. Check `contig_attribution_mode` before interpreting those fractions.

For the marker model, `contamination_type` is `clean` below the reporting threshold, normally 10. At or above that threshold it reports a likely source:

| Value | Evidence |
| --- | --- |
| `cellular` | Cellular markers or coherent cellular-lineage assignments on contigs. |
| `mixed_viral` | Conflicting viral order or family assignments. |
| `phage` | Phage markers or matches to phage, PPV, PLV or virophage references. |
| `duplication` | Elevated marker copy numbers. |
| `uncertain` | No source could be resolved from the available evidence. |

`clean` means that the estimate is below the threshold, not that contamination has been excluded. A `mixed_viral` signal is reported as `uncertain` when at least three contigs carry viral evidence and none has a coherent cellular lineage. Sparse viral references can produce this pattern as well as mixtures. An unavailable estimate is also `uncertain`. So is an estimate that rests on mixed contigs alone, that is, when every contig with coherent cellular evidence also carries viral evidence and the phage, duplication and viral-mixture signals do not name the source. When one of those signals names the source, the type reports it.

### Limits of the bundled contamination model

For the marker model, `contamination_model_reliability` and `contamination_query_support` expose the limits of the bundled model. Its training record describes NCLDV examples with sensitive-mode features. This record does not establish calibration for Mirusviricota, PPV, other search modes, or the corrected feature definitions. An exploratory model value is not a validated contamination percentage.

The mean of the maximum leaf value in each bundled tree is 58.20%. The ensemble cannot predict above that value, even if most of the input is contamination. Clipping the output at 100% does not remove this limit.

Accuracy depends on the test genomes, donors, fragmentation, reference overlap, and search settings. An accuracy claim needs the model hash, input provenance, source-genome split, feature definition, and per-scenario results. A clean-genome test alone cannot establish sensitivity to contamination.

## Marker duplication

`order_dup` and `gvog8_dup` report the average copies per detected marker in their respective panels. A value near 1 indicates mostly single-copy markers. Higher values can reflect mixed populations, assembly duplication, or gene duplication and need inspection. A value of 0 can indicate no detected markers. Neither field gives a direct foreign-DNA percentage.

Two genomes in one bin can recover most expected markers while inflating copy counts. Review duplication together with completeness, contamination, and taxonomic support. The [output reference](../reference/output.md) lists the supporting fields.
