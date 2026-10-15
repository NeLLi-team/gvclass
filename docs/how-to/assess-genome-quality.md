# Assess genome quality

Review the taxonomy, marker counts, and quality estimates together. Completeness
and contamination estimates depend on marker evidence and reference data. One
score cannot certify genome quality.

The commands below use `example_results/` from the
[getting-started tutorial](../tutorials/getting-started.md). Replace that directory
with your results directory.

## 1. Read the assignment and estimates together

```bash
cut -f1,3,5,16,17,20,21,22,25 example_results/gvclass_summary.tsv
```

The selected columns show the query, lineage, taxonomy confidence, completeness
with its query support and model reliability, contamination with its query
support, and contamination type. Check the lineage and its marker support first;
the estimates depend on the assignment. Confirm that the accepted lineage fits
the viral target you expect.

`estimated_completeness` is a marker-based estimate on a 0–100 scale. A value of
100 does not establish that every part of the genome has been assembled.
Compare it with the marker counts and duplication fields. Inspect low-support or
sparsely represented lineages separately, and compare assembly length, missing
genes, contig ends, and read coverage when those data are available.

Two fields describe the support for an estimate:

| Field | What it describes |
| --- | --- |
| `completeness_query_support`, `contamination_query_support` | The evidence in this query. Treat an unavailable estimate (`NaN`) as unknown, not zero. |
| `completeness_model_reliability`, `contamination_model_reliability` | The validation of the model. Bundled models report `unvalidated_legacy`, which marks an exploratory estimate. A historical model score is not a claim of independent accuracy. |

Use `completeness_model_group`, `estimated_completeness_strategy`, and
`estimated_contamination_strategy` to identify the reference group and the
calculation. See [Quality metrics](../explanation/quality-metrics.md) for the
estimators and their limits.

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
a higher threshold. Interpret a missing estimate as unavailable, not zero. Do not
use a low numeric score to dismiss strong conflicting evidence: the bundled
contamination model cannot predict above about 58.20%, so it can underestimate
heavily contaminated inputs.

Inspect marker counts and duplication alongside the estimate:

```bash
cut -f1,15,26-33 example_results/gvclass_summary.tsv
```

This includes `order_dup`, GVOG4 and GVOG8 counts and duplication, and the BUSCO
and UNI56 cellular panels. A duplication factor above 1 means multiple copies
were found for at least some detected markers. Check their contig locations and
annotations before deciding whether they represent mixed genomes, duplicated
sequence, or genuine gene copies.

Review cellular marker hits in their contig context. Gene transfer and integrated
viral sequences, including putative endogenous viral elements, can produce such
hits without a separate contaminating genome. High duplication can identify a
mixed bin that has a low contamination estimate. Confirm the cause with read
coverage, gene context, and contig links when they are available.

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
`contig_attribution_mode`. Protein-only contig assignments and lengths can be
estimates.

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

The archive can also contain `<query>/stats/<query>.contamination_review.tsv`.
This file lists contigs with mixed viral and coherent cellular evidence. It is
written whenever that evidence exists, even when the estimate is low or
unavailable. Inspect the contig boundaries and the `length_basis`. Do not remove
a whole contig solely because it has a review or candidate label.

Keep the scores, support fields, marker evidence, and the reason for each
curation decision. If the evidence is insufficient, record the quality as
unresolved.

All columns and archive paths are listed in the
[output reference](../reference/output.md).
