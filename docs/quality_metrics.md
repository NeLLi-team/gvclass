# Quality Metrics

This document is the developer-facing reference for GVClass completeness and contamination metrics.
It complements the user-facing summaries in the README by describing the actual runtime formulas, the model-selection logic, and the offline training provenance for the shipped resources.

## Source Of Truth

- Runtime completeness logic:
  [src/core/summarize_full.py](/home/fschulz/dev/software/gvclass/src/core/summarize_full.py),
  [src/core/weighted_completeness.py](/home/fschulz/dev/software/gvclass/src/core/weighted_completeness.py),
  [src/core/novelty_completeness.py](/home/fschulz/dev/software/gvclass/src/core/novelty_completeness.py)
- Runtime contamination logic:
  [src/core/contamination_scoring.py](/home/fschulz/dev/software/gvclass/src/core/contamination_scoring.py)
- Offline completeness resource builder:
  [src/bin/build_novelty_completeness_resources.py](/home/fschulz/dev/software/gvclass/src/bin/build_novelty_completeness_resources.py)
- Offline contamination benchmark/export workflow:
  [src/bin/benchmark_contamination_methods_faa.py](/home/fschulz/dev/software/gvclass/src/bin/benchmark_contamination_methods_faa.py)

## Completeness Metrics

### Legacy Order Completeness

For an assigned order, GVClass loads the expected orthogroup list from `order_completeness.tab`.

Definitions:

- `present_ogs` = number of expected orthogroups with count `> 0`
- `total_hits` = sum of counts across the expected orthogroups
- `n_expected` = number of expected orthogroups for the order

Formulas:

```text
order_completeness_raw = (present_ogs / n_expected) * 100
order_dup = total_hits / present_ogs
order_completeness = clamp((order_completeness_raw / baseline_mean) * 100, 0, 100)
```

Fields:

- `order_completeness_raw`: direct marker recovery fraction
- `order_dup`: duplication factor for the expected order marker set
- `order_completeness`: baseline-normalized legacy completeness
- `order_completeness_baseline_mean` and `order_completeness_baseline_std`: the reference recovery baseline used for normalization
- `order_completeness_reference_order`: the order whose baseline was used
- `order_completeness_strategy`: currently `order_baseline_ratio_v1`

### Weighted Order Completeness

The weighted completeness path uses `marker_stats.tsv` to assign per-marker weights.
The shipped calculator uses:

- weight strategy: `information_theoretic`
- adaptive scaling: enabled
- outlier detection: `IsolationForest`

Base marker weights:

```text
conservation_rate = percent_genomes_with_marker / 100
base_weight = -log(1 - conservation_rate + 1e-10)
```

Weighted completeness:

```text
weighted_score = sum(weight_i * presence_i) / sum(weight_i)
weighted_order_completeness_raw = weighted_score * 100
```

Confidence score:

```text
coverage_factor = n_present / n_expected
weight_diversity = clamp(1 - std(weights) / mean(weights), 0, 1)
order_confidence_score = ((0.7 * coverage_factor) + (0.3 * weight_diversity)) * 100
```

`order_weighted_completeness` is the weighted raw score normalized through the same order baseline used for the legacy completeness path.

If weighted scoring fails, GVClass falls back to:

- `order_weighted_completeness_raw = order_completeness_raw`
- `order_confidence_score = 50.0`

### Novelty-Aware Completeness

The novelty-aware scorer is runtime-only.
Its tier definitions, baselines, and strategy-3 model bundle are built offline and shipped in `resources/`.

#### Group Selection

- Prefer `family` tiers when a family-specific tier set exists for the assigned order.
- Otherwise fall back to the `order` tier set.
- If neither exists, the novelty-aware path returns default zero/`unavailable` values and `estimated_completeness` falls back according to `--completeness-mode`.

#### Strategy 2 Tier Score

Each group has three marker tiers:

- `core`
- `shared`
- `accessory`

Presence fractions:

```text
core_fraction = present(core) / len(core)
shared_fraction = present(shared) / len(shared)
accessory_fraction = present(accessory) / len(accessory)
```

The shipped resource builder uses:

```text
weight_core = 0.7
weight_shared = 0.2
weight_accessory = 0.1
```

Strategy-2 raw score:

```text
strategy2_raw =
  100 * (
    weight_core * core_fraction +
    weight_shared * shared_fraction +
    weight_accessory * accessory_fraction
  ) / (sum of weights that are actually present as tiers)
```

Normalization:

```text
strategy2_normalized = min(100, (strategy2_raw / baseline_mean) * 100)
```

#### Support And OOD Diagnostics

Informative markers are `core + shared`.

Support score:

```text
informative_present = round(core_fraction * n_core) + round(shared_fraction * n_shared)
informative_fraction = informative_present / informative_marker_count
ref_support = min(1, baseline_n_refs / 10)
support_score = ((0.85 * informative_fraction) + (0.15 * ref_support)) * 100
```

OOD flag rules:

- `unassigned`
- `no_informative_markers`
- `zero_informative_hits`
- `low_support`
- `low_reference_support`
- `ok`

These appear in:

- `order_completeness_v2_support_score`
- `order_completeness_v2_ood_flag`
- `order_completeness_v2_reference_group`
- `order_completeness_v2_validation_mode`
- `order_completeness_v2_informative_fraction`

#### Strategy 3 Model And OOD Shrinkage

If a model exists for the assigned order, GVClass predicts a calibrated completeness score from these inputs:

- `raw_score`
- `strategy1_score`
- `strategy2_raw_score`
- `strategy2_score`
- `strategy2_baseline_mean`
- `strategy2_baseline_n_refs`
- `core_fraction`
- `shared_fraction`
- `accessory_fraction`
- `n_expected_markers`
- `n_core_markers`
- `n_shared_markers`
- `n_accessory_markers`
- `group_is_family`

The raw strategy-3 prediction is then bounded to `0..100` and optionally shrunk toward the conservative floor:

```text
conservative_floor = max(strategy1_score, strategy2_score)
```

The shrinkage strength depends on:

- `ood_flag`
- `support_score`
- `baseline_n_refs`

This final bounded/shrunk score becomes `order_completeness_v2`.

#### Surfaced Primary Completeness Estimate

`estimated_completeness` is not an independent metric.
It is a selector:

- `legacy` mode -> `estimated_completeness = order_completeness`
- `novelty-aware` mode -> `estimated_completeness = order_completeness_v2`

`estimated_completeness_strategy` reports which branch was surfaced.

## Contamination Metrics

### Rule-Based Diagnostic Score

The rule-based path is still computed at runtime because it contributes useful interpretable diagnostics and forms part of the feature set for the trained contamination model.
It is no longer the production contamination estimate.

The rule-based path combines four components:

- `cellular_signal`
- `phage_signal`
- `duplication_signal`
- `viral_mixture_signal`

All component scores are clipped to `0..100`.

#### Blast/Tree Inputs

Only strong best hits are considered for the blast-derived features:

- identity `>= 60`
- score `>= 80`

Tree-derived features use majority/secondary fractions across the assigned taxonomic counters.

#### Cellular Signal

```text
cellular_signal =
  max(0, (cellular_unique - 1) * 9.0)
  + max(0, cellular_total - cellular_unique) * 1.5
  + strong_cellular_hit_count * 3.0
  + nonviral_best_hit_fraction * 0.20
  + dominant_nonviral_lineage_fraction * 0.10    if strong_cellular_hit_count >= 3
```

#### Phage Signal

```text
phage_signal =
  (phage_unique - 1) * 10.0                      if phage_unique >= 2
  + max(0, phage_total - phage_unique) * 2.0
  + strong_phage_hit_count * 4.0
```

#### Duplication Signal

```text
duplication_signal =
  max(0, order_dup - 1.20) * 45.0
  + max(0, gvog8_dup - 1.30) * 30.0
```

#### Viral-Mixture Signal

```text
viral_mixture_signal =
  max(0, order_secondary_fraction - 10.0) * 1.3
  + max(0, family_secondary_fraction - 15.0) * 0.8
  + foreign_viral_order_hit_count * 5.0
  + foreign_viral_family_hit_count * 2.5
  + dominant_foreign_viral_order_fraction * 0.15   if foreign_viral_order_hit_count >= 3
```

#### Final Rule-Based Score

```text
contamination_score_v1 =
  0.38 * cellular_signal
  + 0.12 * phage_signal
  + 0.20 * duplication_signal
  + 0.30 * viral_mixture_signal
```

Severity bins:

- `clean` if `< 8`
- `low` if `>= 8`
- `moderate` if `>= 20`
- `high` if `>= 45`

`contamination_source_v1` is the largest component unless the final score is `< 8`, in which case it is forced to `none`.

Fields:

- `contamination_score_v1`
- `contamination_flag_v1`
- `contamination_source_v1`
- `contamination_cellular_signal_v1`
- `contamination_phage_signal_v1`
- `contamination_duplication_signal_v1`
- `contamination_viral_mixture_signal_v1`
- `contamination_nonviral_hit_fraction_v1`

### Contig-Aware Suspicious Sequence Burden

`suspicious_bp_fraction_v2` and `suspicious_contig_count_v2` are feature-extraction metrics, not the final production contamination estimate.

Contig lengths come from:

- the real nucleotide contig lengths when an `.fna` is available
- otherwise a coding-length proxy from `.faa` sequence lengths times `3`

A contig is labeled suspicious if any of these triggers fire:

- `cellular_marker_count >= 2`
- `len(cellular_hits) >= 3` and `nonviral_fraction >= 0.5` and `viral_marker_count <= 1`
- `phage_marker_count >= 2` and `viral_marker_count <= 1`
- `len(phage_hits) >= 3` and `nonviral_fraction >= 0.5` and `viral_marker_count <= 1`
- `len(foreign_viral) >= 2` and `foreign_viral_fraction >= 0.6` and `viral_marker_count >= 2`

The summary fields report:

- the number of suspicious contigs
- the fraction of total assembly bp assigned to suspicious contigs

### Production Contamination Estimate

The shipped production estimate is the trained model in `src/bundled_models/contamination_model.joblib`.

The shipped bundle currently reports:

- `model_name = hist_gbm`
- runtime strategy label = `hist_gbm_v1`

Model input features:

- `contamination_score_v1`
- `contamination_cellular_signal_v1`
- `contamination_phage_signal_v1`
- `contamination_duplication_signal_v1`
- `contamination_viral_mixture_signal_v1`
- `contamination_nonviral_hit_fraction_v1`
- `order_dup`
- `gvog8_dup`
- `cellular_unique`
- `cellular_total`
- `phage_unique`
- `phage_total`
- `contigs`
- `suspicious_bp_fraction_v2`
- `suspicious_contig_count_v2`
- `estimated_completeness`

The model output is clipped to `0..100` and written to:

- `estimated_contamination`
- `estimated_contamination_strategy`

## Training Provenance

### Completeness Resources

Offline completeness resources are built with:

- [src/bin/build_novelty_completeness_resources.py](/home/fschulz/dev/software/gvclass/src/bin/build_novelty_completeness_resources.py)

The builder:

- derives order and family marker profiles from the reference FAA database
- constructs `core`, `shared`, and `accessory` tiers
- computes family/order baselines
- trains order-specific strategy-3 regressors
- uses leave-one-family-out validation where possible

The shipped resource builder config is:

- `core_prevalence = 0.9`
- `shared_prevalence = 0.65`
- `accessory_prevalence = 0.5`
- `core_max_copy = 1.2`
- `shared_max_copy = 1.5`
- `family_min_ref_genomes = 3`
- `family_min_informative_markers = 2`
- `weight_core = 0.7`
- `weight_shared = 0.2`
- `weight_accessory = 0.1`

### Contamination Model

Offline contamination model training/export is handled by:

- [src/bin/benchmark_contamination_methods_faa.py](/home/fschulz/dev/software/gvclass/src/bin/benchmark_contamination_methods_faa.py)

The benchmark design:

- one-contig isolate references from 7 families
- train/test split by genome within family
- synthetic scenarios:
  - `clean`
  - `hgt_cellular_1gene`
  - `hgt_cellular_5genes`
  - `cellular_contig_10`
  - `cellular_contig_25`
  - `viral_divergent_10`
  - `viral_divergent_25`
  - `viral_related_10`
  - `viral_related_25`
- HGT-like scenarios are labeled as `0` true contamination

Model family selection:

- candidate models:
  - `random_forest`
  - `extra_trees`
  - `hist_gbm`
- selection objective:
  - grouped cross-validation MAE on held-out genomes
- threshold optimization:
  - best F1 over thresholds `1.0..50.0` in `0.5` increments

The trained production bundle is exported from the selected winning family on the full FAA feature table.

## Interpreting The Shipped Examples

The bundled example outputs are useful for checking which branch is active:

- in novelty-aware runs, `estimated_completeness_strategy = novelty_aware_v1`
- in legacy runs, `estimated_completeness_strategy = order_baseline_ratio_v1`
- in all current production runs, `estimated_contamination_strategy = hist_gbm_v1`

Current example values from the novelty-aware run:

- `AC3300027503___Ga0255182_1000024`
  - `estimated_completeness = 79.62`
  - `estimated_contamination = 21.66`
  - `contamination_score_v1 = 7.98`
- `GVMAG-S-1096109-37`
  - `estimated_completeness = 89.88`
  - `estimated_contamination = 3.43`
  - `contamination_score_v1 = 0.80`
- `PkV-RF01`
  - `estimated_completeness = 95.15`
  - `estimated_contamination = 0.56`
  - `contamination_score_v1 = 26.22`

These examples illustrate an important point:

- `contamination_score_v1` is not the primary contamination estimate
- `estimated_contamination` can diverge substantially from the rule-based precursor because the trained model uses the full feature set rather than a single linear threshold on the rule-based score
