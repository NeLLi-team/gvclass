# Assess genome quality

Review the taxonomy, marker counts, and quality estimates together. Completeness
and contamination estimates depend on marker evidence and reference data.

The commands below use `example_results/` from the
[getting-started tutorial](../tutorials/getting-started.md). Replace that directory
with your results directory.

## 1. Read the assignment and estimates together

```bash
cut -f1,2,4,15-18 example_results/gvclass_summary.tsv
```

The selected columns show the query, lineage, taxonomy confidence, completeness,
model reliability, contamination, and contamination type. Check the lineage and
its marker support first; the estimates depend on the assignment.

`estimated_completeness` is a marker-based estimate on a 0–100 scale. A value of
100 does not establish that every part of the genome has been assembled.
`completeness_model_reliability` describes the per-order model, not the query:

| Value | Interpretation |
| --- | --- |
| `advisory_only` | Model validation is missing or insufficient, or a fallback estimate is used. |
| `moderate` | The model meets the intermediate hold-out fit threshold. |
| `high` | The model meets the higher hold-out fit threshold. |

See [Quality metrics](../explanation/quality-metrics.md) for the estimators and
reliability thresholds.

## 2. Review the contamination evidence

`estimated_contamination` is the trained model's estimate. The type label
summarises the detected signal:

| `contamination_type` | Signal to review |
| --- | --- |
| `clean` | Estimate below the reporting threshold; not proof that contamination is absent. |
| `cellular` | Cellular marker or contig evidence. |
| `mixed_viral` | Marker placements suggest more than one viral lineage. |
| `phage` | Phage markers or matches to phage, PPV, PLV or virophage references. |
| `duplication` | Excess marker copies. |
| `uncertain` | Unresolved signal or an unavailable numeric estimate. |

The reporting threshold is at least 10 and can be higher if the model specifies
a higher threshold. Interpret a missing estimate as unavailable, not zero.

Inspect marker counts and duplication alongside the estimate:

```bash
cut -f1,14,19-26 example_results/gvclass_summary.tsv
```

This includes `order_dup`, GVOG4 and GVOG8 counts and duplication, and the BUSCO
and UNI56 cellular panels. A duplication factor above 1 means multiple copies
were found for at least some detected markers. Check their contig locations and
annotations before deciding whether they represent mixed genomes, duplicated
sequence, or genuine gene copies.

Review cellular marker hits in their contig context. Gene transfer and integrated
viral sequences can produce such hits without a separate contaminating genome.

## 3. Inspect the archived diagnostics

Extract the extended table:

```bash
mkdir -p example_results/diagnostics
tar -xzf example_results/gvclass_summary.extended.tar.gz \
    -C example_results/diagnostics
head -n 4 example_results/diagnostics/gvclass_summary.extended.tsv
```

It includes contig-level fields such as `cellular_coherent_contig_count`,
`cellular_lineage_purity_median`, `viral_bearing_contig_count`, and
`contig_attribution_mode`.

List the files for one query:

```bash
tar -tzf example_results/GVMAG-S-1096109-37.tar.gz
```

For flagged queries, the archive can contain
`<query>/stats/<query>.contamination_candidates.tsv`. This file is written when
the estimate reaches the reporting threshold, the type is interpretable, and
candidate contigs are found. Review those contigs against gene annotations,
coverage, and assembly context before removing them. The absence of a candidate
file does not prove that a bin is clean.

All columns and archive paths are listed in the
[output reference](../reference/output.md).
