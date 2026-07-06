# Completeness and contamination

Standard genome quality checks assume a tidy single-copy core and treat anything host-like as contamination. This does not apply to giant virus genomes. Nucleocytoviricota genomes are gene-rich, and they routinely carry genes that look eukaryotic, so a checker built for bacterial or eukaryotic MAGs would flag almost every real giant virus as contaminated. GVClass measures quality with metrics tuned to that biology. Each bin gets two headline numbers, `estimated_completeness` and `estimated_contamination`, plus supporting columns that say how far to trust them.

## How complete is the genome

`estimated_completeness` is the percentage of the expected gene complement recovered for the lineage the bin was assigned to. By default it comes from a novelty-aware, per-lineage model. Rather than scoring every bin against one fixed marker list, the model conditions the expectation on the assigned lineage and on how novel that lineage is relative to the reference set. The `--completeness-mode legacy` flag selects a fixed-panel calculation instead; `novelty-aware` is the default. It is the only completeness field in the main table.

The companion column `completeness_model_reliability` is easy to misread. It takes three values, `advisory_only`, `moderate`, and `high`, and it reflects the hold-out R^2 of the completeness model for the bin's assigned order. It describes the model, not the genome. A bin can score `estimated_completeness` 100.00 and still carry `advisory_only` reliability when the model for its order validated poorly in hold-out. The bundled example shows this directly. PkV-RF01 returns `estimated_completeness` 100.00 with `8/8` core GVOG markers, yet its reliability is `advisory_only`, because the completeness model for its order, Imitervirales, validates only at the advisory tier. The 10-contig Pimascovirales bin GVMAG-S-1096109-37 returns 82.86 at `moderate` reliability. Read the percentage through its tier: a high number under `advisory_only` is a plausible estimate, not a measured fact.

!!! warning

    `completeness_model_reliability` is a property of the order model, not of the genome in front of you. An `advisory_only` value flags an order where the model lacks hold-out support, so treat even a 100.00 completeness in that order as an estimate to confirm, not a guarantee.

## What contamination means for a giant virus

Here the giant virus case departs sharply from cellular genome QC. Some Nucleocytoviricota genomes carry hundreds of eukaryote-like genes acquired by horizontal gene transfer from their hosts and prey, spanning metabolism, cytoskeletal proteins, nutrient transport, and DNA repair. A contamination checker built for cellular MAGs would read those host-like genes as foreign DNA and condemn the genome. GVClass does not count them. The HGT-derived gene content is a genuine feature of giant virus biology, and penalizing it would reject nearly every legitimate genome.

To separate real cellular contamination from this HGT background, GVClass restricts the cellular signal to conserved cellular marker genes that giant viruses obligately lack: the eukaryotic BUSCO set (255 markers) and the universal cellular UNI56 set (56 markers). Giant viruses lack a complete cellular housekeeping repertoire and depend on the host for it, so a clean giant virus bin carries essentially none of these markers. When they do appear, they point to cellular DNA that co-binned with the virus rather than to viral HGT. That is the contamination GVClass is built to catch. The [marker reference](../reference/markers.md) lists each panel and its size.

!!! note

    GVClass does not penalize the eukaryote-like genes that giant viruses acquire by HGT. Only the conserved cellular markers (BUSCO and UNI56) that a virus cannot legitimately carry feed `estimated_contamination`.

## The contamination estimate

`estimated_contamination` is the headline contamination figure, the output of a trained model rather than a marker ratio. The model is an `ExtraTreesRegressor` (identifier `extra_trees_v1`) shipped in the resource bundle at `resources/contamination/model.joblib`. It is trained on features computed in sensitive mode, the default search setting. On the real-contig benchmark it reaches a mean absolute error of 3.92%, and on clean bins its mean predicted contamination is 0.14%, close to zero. The bundled clean examples land where that benchmark predicts: PkV-RF01 and AC3300027503___Ga0255182_1000024 both return `estimated_contamination` 0.00, and the 10-contig GVMAG-S-1096109-37 returns 0.08.

## Naming the source of contamination

When `estimated_contamination` reaches 10 or higher, GVClass fills in `contamination_type` to suggest where the extra DNA came from. The categories are `clean`, `cellular`, `mixed_viral`, `phage`, `duplication`, and `uncertain`. Two of them carry the conceptual weight of this page. `cellular` means the bin picked up host or co-occurring microbial DNA, evidenced by those mainly single-copy BUSCO or UNI56 markers a virus should not have. `mixed_viral` means two or more viral orders are mixed on one bin; it is multi-order viral mixing, not host gene content, and not the HGT genes described above. The remaining two are narrower: `phage` flags bacteriophage markers (the geNomad panel) on the bin, and `duplication` flags inflated marker copy numbers.

One refinement protects novel lineages from a false `mixed_viral` call. When a bin shows no coherent cellular contigs, three or more viral-bearing contigs, and a `viral_mixture` signature, that pattern can equally mean a genuinely novel virus whose markers do not yet match a single reference order. GVClass downgrades that case from `mixed_viral` to `uncertain`, which routes it to manual review rather than rejection. Notably, the per-contig evidence behind these calls (`cellular_coherent_*`, `cellular_lineage_purity_median`, `viral_bearing_contig_count`, `contig_attribution_mode`) is written to the extended summary tables, keeping the main table readable.

## Duplication as complementary QC

Two duplication columns add a second, orthogonal view that does not depend on the contamination model. `order_dup` is the average copy number of the order-level markers expected for the lineage, and `gvog8_dup` is the duplication factor across the eight core NCLDV markers. A single clean genome should carry one copy of each, so both sit near 1. Values above roughly 2 point to multiple populations, an assembly chimera, or a mixed bin; values below about 1.5 are typically clean. Duplication catches a failure mode the percentages can miss: a bin can look complete and score low contamination yet still hold two near-complete genomes of the same lineage, which inflates the duplication factors before it surfaces anywhere else.

## Reading the metrics together

Together these columns describe a bin from complementary angles. High `estimated_completeness` at a reliability tier you trust, low `estimated_contamination`, `contamination_type` `clean`, and duplication factors near 1 is the signature of a high-quality GVMAG. Low completeness, contamination at or above 10, a non-`clean` type, or duplication above 2 each flag a bin for manual curation. The [assess genome quality](../how-to/assess-genome-quality.md) how-to turns these thresholds into a step-by-step triage you can run across a batch of bins, and [the output reference](../reference/output.md) documents every column these metrics live in.
