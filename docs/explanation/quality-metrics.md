# Completeness and contamination

GVClass estimates genome quality from marker recovery, marker duplication, reference matches, and contig-level evidence. The main fields are `estimated_completeness`, `completeness_model_reliability`, `estimated_contamination`, and `contamination_type`.

For commands to inspect these fields, see [Assess genome quality](../how-to/assess-genome-quality.md).

## Completeness

`estimated_completeness` is a marker-based estimate on a 0–100 scale. It compares recovered markers with the expected complement for a reference lineage. A value of 100 does not demonstrate that the full nucleotide sequence has been recovered.

The default `novelty-aware` mode uses lineage-specific marker tiers and, where available, an order-specific prediction model. If a model's recorded hold-out R² is below 0.5, GVClass uses the marker-tier estimate. `--completeness-mode legacy` uses the order-marker recovery ratio scaled to a reference baseline.

`completeness_model_reliability` describes the available model and its recorded hold-out R²:

| Value | Interpretation |
| --- | --- |
| `advisory_only` | Model or validation information is missing, or R² is below 0.5. |
| `moderate` | R² is at least 0.5 and below 0.7. |
| `high` | R² is at least 0.7. |

These tiers describe model validation, not confidence intervals for individual genomes. A high completeness estimate with `advisory_only` reliability needs support from marker counts and other assembly evidence. Strategy and support fields are recorded in `<query>/<query>.summary.tab` inside each query's `.tar.gz` archive.

## Contamination

`estimated_contamination` is the prediction of the bundled ExtraTrees regressor. Its inputs include cellular and phage marker signals, marker duplication, taxonomic disagreement, and the distribution of reference matches across contigs. The archived per-query summary records the model identifier in `estimated_contamination_strategy`.

Giant viruses can carry genes acquired from cellular organisms. A cellular-like protein or a BUSCO/UNI56 hit alone is insufficient to identify a contaminating contig. GVClass combines these signals with contig-level lineage evidence. Inspect the implicated contigs before removing DNA from a bin.

`contamination_type` is `clean` below the reporting threshold, normally 10. At or above that threshold it reports a likely source:

| Value | Evidence |
| --- | --- |
| `cellular` | Cellular markers or coherent cellular-lineage assignments on contigs. |
| `mixed_viral` | Conflicting viral order or family assignments. |
| `phage` | Phage markers or matches to phage, PPV, PLV or virophage references. |
| `duplication` | Elevated marker copy numbers. |
| `uncertain` | No source could be resolved from the available evidence. |

`clean` means that the estimate is below the threshold, not that contamination has been excluded. A `mixed_viral` signal is reported as `uncertain` when at least three contigs carry viral evidence and none has a coherent cellular lineage. Sparse viral references can produce this pattern as well as mixtures.

## Marker duplication

`order_dup` and `gvog8_dup` report the average copies per detected marker in their respective panels. A value near 1 indicates mostly single-copy markers. Higher values can reflect mixed populations, assembly duplication, or gene duplication and need inspection. A value of 0 can indicate no detected markers.

Two genomes in one bin can recover most expected markers while inflating copy counts. Review duplication together with completeness, contamination, and taxonomic support. The [output reference](../reference/output.md) lists the supporting fields.
