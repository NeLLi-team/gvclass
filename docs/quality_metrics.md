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
  [src/bin/benchmark_real_contig_contamination.py](/home/fschulz/dev/software/gvclass/src/bin/benchmark_real_contig_contamination.py)

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
See the companion model card at
[src/bundled_models/contamination_model.yaml](/home/fschulz/dev/software/gvclass/src/bundled_models/contamination_model.yaml)
for the full provenance, integrity digest, and threshold.

The shipped bundle currently reports (v1.4.3):

- `model_name = extra_trees`
- runtime strategy label = `extra_trees_v1`
- training profile = `sensitive_mode_features` (features generated with
  `sensitive_mode=true`, E=1e-5, GA/TC/NC cutoffs disabled)
- held-out MAE on the real-contig benchmark test set = `3.92`
- mean predicted contamination on clean bins = `0.14%`
- integrity gate: SHA-256 of the `.joblib` file is verified against
  `CONTAMINATION_MODEL_SHA256` in
  [src/core/contamination_scoring.py](/home/fschulz/dev/software/gvclass/src/core/contamination_scoring.py)
  before `joblib.load` is allowed to run

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

- [src/bin/benchmark_real_contig_contamination.py](/home/fschulz/dev/software/gvclass/src/bin/benchmark_real_contig_contamination.py)

The benchmark design:

- 30 real giant-virus base genomes with complete taxonomy coverage
- per-base one-contig mixtures with cellular or foreign-viral contigs
  drawn from independent donor genomes; donor bp targets span the
  `clean`, `viral_pair`, `viral_triple`, `cellular_single_*`, and
  `cellular_multi*` scenarios listed in the run manifest
- train/test split by `base_accession` (no genome appears in both
  splits) with a GroupKFold over base genomes for candidate selection
- HGT-like single-gene cellular contributions are NOT labeled as
  contamination — the model learns to predict cellular CONTIG burden,
  not the ubiquitous eukaryote-like gene content that giant viruses
  naturally carry via horizontal gene transfer (see the explanation of
  `contamination_type = mixed_viral` below)

Feature-generation mode for v1.4.3:

- `sensitive_mode=true` is active through the HMM step, matching the
  runtime default. Training and inference therefore see the same
  feature distribution.

Model family selection:

- candidate models:
  - `random_forest` (`n_estimators=400, min_samples_leaf=2`)
  - `extra_trees` (`n_estimators=500, min_samples_leaf=2`)
  - `hist_gbm` (`max_depth=6, min_samples_leaf=5`)
- selection objective:
  - grouped cross-validation MAE on held-out base genomes
- threshold optimization:
  - best F1 over thresholds `1.0..50.0` in `0.5` increments

The trained production bundle is exported from the winning family on
the full feature table. v1.4.3 selected `extra_trees`.

### Per-contig taxonomic-purity classifier (v1.4.3 Phase 2)

Added in v1.4.3 to distinguish two bin shapes that the rule-based
`viral_mixture_signal` collapsed into one label:

1. **Novel-virus sparse-reference case** — a new giant-virus lineage
   whose markers tree-resolve to scattered eukaryotic HGT-recipients
   because no close NCLDV reference is available. EUK-leaning
   placements sit on the same contigs as NCLDV markers.
2. **Real host-contig contamination** — a eukaryotic contig
   co-assembled with the viral bin. That contig carries ribosomal /
   translation markers whose tree placements all agree on ONE
   cellular lineage, and those contigs have NO giant-virus markers.

The classifier runs per-contig, aggregating each protein's dominant
tree-NN lineage (drawn from `stats/<query>.tree_nn`) across every
marker that protein hit. A contig is classified as:

- **`viral_bearing`** if ≥ 1 protein's dominant NN lineage is in
  `{NCLDV, MIRUS, MRYA, VP, PLV}` — these contigs are genuinely
  viral and the classifier refuses to call them cellular even if
  they carry scattered EUK markers alongside.
- **`cellular_coherent`** if the contig has 0 viral proteins AND
  ≥ `MIN_CELLULAR_MARKERS = 3` cellular proteins AND
  ≥ `PURITY_THRESHOLD = 60 %` of those cellular proteins agree on
  the same top-vote `(domain, order)` AND the median BLAST identity
  of the cellular hits on that contig is ≥ `MIN_CELLULAR_IDENTITY
  = 70 %`. The identity floor is the key guard against novel-virus
  HGT contigs whose cellular-resolving genes sit at 25-55 % identity
  to their donors; only genuine host contigs land at ≥ 70 % across
  multiple ribosomal / translation markers.
- **`ambiguous`** otherwise (fragmentary contigs, scattered
  lineages, low identity, mixed HGT signals).

Six new summary columns expose the classifier state:

- `cellular_coherent_contig_count` — count of contigs that passed all
  four gates. Primary host-contig-contamination indicator.
- `cellular_coherent_protein_fraction` — fraction of
  marker-bearing proteins on `cellular_coherent` contigs. Primary ML
  feature for the size of the cellular burden; invariant across
  `.fna` and `.faa` inputs.
- `cellular_coherent_bp_fraction` — fraction of total assembly bp on
  `cellular_coherent` contigs. Diagnostic only; omitted from the ML
  feature set because the `len(seq) * 3` proxy used in `.faa` runs
  would mis-scale it.
- `cellular_lineage_purity_median` — median purity across the
  `cellular_coherent` contigs (how sharp the cellular call is
  when it fires).
- `cellular_hit_identity_median` — median BLAST identity of the
  cellular-hit proteins on those contigs (genuine host contigs
  carry high-identity hits).
- `viral_bearing_contig_count` — dilution signal; a bin with many
  viral-bearing contigs and zero cellular-coherent contigs is a
  novel-virus candidate rather than a contaminated bin.

#### Fallback contract

The classifier has two degenerate regimes in which it emits all six
features at their neutral values (`0`) and the legacy
`cellular_marker_count >= 2` rule becomes the sole cellular trigger:

- **Low-data regime** — fewer than
  `MIN_MARKER_PROTEINS_FOR_FRACTION = 5` marker-bearing proteins in
  the bin. Per-contig aggregation is not informative below this floor.
- **Weak-attribution regime** — `.faa`-only input whose protein IDs
  do not encode the source contig, measured as
  `n_distinct_contigs / n_proteins >= 0.8`. In this regime every
  protein maps to its own pseudo-contig and the per-contig classifier
  is structurally meaningless. A new `contig_attribution_mode` column
  reports which regime fired (`fna_gene_calling`,
  `faa_contig_encoded`, or `faa_weak_per_protein`) and an INFO log
  is emitted when the fallback activates.

### Interpreting `contamination_type`

The per-query `contamination_type` in the summary is a categorical
label that names the dominant rule-based signal when the final
`estimated_contamination` clears the reporting threshold:

- `clean` — below threshold
- `cellular` — the `cellular_signal` component dominates; at least one
  eukaryotic / bacterial / archaeal / plastid / mitochondrial CONTIG
  appears to carry conserved housekeeping markers (BUSCO odb10
  eukaryotic set + UNI56 prokaryotic translation machinery). Giant
  viruses are not expected to carry these markers at all — they are
  ribosomal proteins, tRNA synthetases, and other core translation
  genes that viruses obligately scavenge from their hosts. Two or more
  such markers on a single contig flag the contig as cellular.
- `phage` — phage-marker burden on contigs otherwise lacking giant
  virus markers
- `duplication` — elevated per-order or GVOG8 marker duplication,
  typically a chimeric assembly rather than a foreign-contig problem
- `mixed_viral` — the `viral_mixture_signal` component dominates;
  multiple distinct giant-virus orders or families appear to be mixed
  in the bin. This is NOT about "has EUK-like genes": giant viruses
  routinely carry hundreds of eukaryote-derived genes via HGT and
  those genes do not count as contamination here. `mixed_viral` is
  triggered when, for example, the assembly carries proteins whose
  tree neighbors land in both Pimascovirales and Imitervirales, or
  where the second-place order/family in the per-marker taxonomy vote
  exceeds a support threshold. See the `viral_mixture_signal` formula
  in [Rule-Based Diagnostic Score](#rule-based-diagnostic-score).
- `uncertain` — threshold cleared but no single signal dominates

Why giant-virus eukaryote-like genes do not inflate the cellular
signal: the `CELLULAR_MODELS` HMM set is deliberately limited to
translation-machinery genes (ribosomal proteins, tRNA synthetases)
that giant viruses do not encode. A Mimivirus translation-factor
homolog, for example, is carried in the giant-virus-specific GVOG
panel rather than in BUSCO / UNI56, and it therefore does not
contribute to `cellular_signal`. Only when a eukaryotic CONTIG
carrying at least two ribosomal / translation markers is co-assembled
with the viral bin does the cellular path fire — a signal that
reliably points at contaminating host sequence rather than legitimate
HGT acquisition.

## Interpreting The Shipped Examples

The bundled example outputs are useful for checking which branch is active:

- in novelty-aware runs, `estimated_completeness_strategy = novelty_aware_v1`
- in legacy runs, `estimated_completeness_strategy = order_baseline_ratio_v1`
- in v1.4.3 production runs, `estimated_contamination_strategy = extra_trees_v1`

Current example values from the sensitive-mode novelty-aware run
(v1.4.3, bundled `example/` directory):

- `AC3300027503___Ga0255182_1000024`
  - `estimated_completeness = 81.73` (quality = advisory_only)
  - `estimated_contamination = 0.00`
  - `contamination_type = clean`
- `GVMAG-S-1096109-37`
  - `estimated_completeness = 82.86` (quality = advisory_only)
  - `estimated_contamination = 26.93`
  - `contamination_type = uncertain` (downgraded from `mixed_viral` by
    the Phase 2 novel-virus-candidate rule — see below)
  - `cellular_coherent_contig_count = 0` (no host-contig evidence)
  - `viral_bearing_contig_count = 3` (three contigs carry NCLDV markers)
  - rule-based signal breakdown (from the `.summary.tab`):
    - `cellular_signal = 3.00` (low — `cellular_unique = 1`,
      `cellular_total = 3`; no CONTIG carries ≥2 translation markers)
    - `phage_signal = 0.00`
    - `duplication_signal = 18.04` (`order_dup = 1.52`,
      `gvog8_dup = 1.43` — elevated but not extreme)
    - `viral_mixture_signal = 69.50` — **dominant**
    - `suspicious_contig_count_v2 = 0` — no contig-level cellular or
      foreign-viral contig was flagged as suspicious
  - `viral_mixture_signal` is driven here by the per-marker tax vote:
    the order-level counter splits roughly NCLDV__Pimascovirales 50%,
    EUK__NA 45%, with a PLV trace. The formula `max(0,
    order_secondary_fraction - 10) * 1.3` contributes ~46 points on
    its own, plus a similar term for family. Note that the
    `order_secondary_fraction` term is agnostic to whether the
    secondary order is viral or cellular; it is therefore better read
    as "low purity of the per-marker tax call" than as
    "multiple giant-virus orders present".
  - **Phase 2 per-contig classifier resolves the ambiguity**: the
    per-contig taxonomic-purity classifier (see the dedicated section
    below) examines each contig in this bin and finds **zero**
    `cellular_coherent` contigs — no contig has the combination of
    coherent cellular lineage, no viral markers, and high-identity
    cellular BLAST hits that would indicate real host contamination.
    The three `viral_bearing` contigs (≥1 NCLDV marker each) dominate
    the bin. With `cellular_coherent_contig_count = 0`,
    `viral_bearing_contig_count = 3`, and
    `contamination_source_v1 = viral_mixture`,
    `_classify_contamination_type` recognises this as the
    novel-virus fingerprint and downgrades the label from
    `mixed_viral` to `uncertain`, so curators see the evidence
    for-and-against contamination rather than a firm mixed-viral
    call. The absolute ML score (26.93 %) is still surfaced because
    the trained model has not yet seen novel-virus examples with
    truth label 0 in its benchmark, so it over-estimates here — a
    future release will expand the benchmark to close that gap.
- `PkV-RF01`
  - `estimated_completeness = 95.43` (quality = advisory_only)
  - `estimated_contamination = 0.09`
  - `contamination_type = clean`

These examples illustrate several important points:

- Reference NCLDV genomes (`PkV-RF01`) predict near-zero contamination,
  as expected for a clean viral reference.
- `estimated_contamination` is NOT driven by the simple presence of
  eukaryote-like genes. Giant viruses carry hundreds of HGT-derived
  eukaryotic genes as part of their legitimate gene content, and the
  contamination pipeline is explicitly designed not to flag those.
- `GVMAG-S-1096109-37` is called `mixed_viral` because the per-marker
  tree neighbors spread across multiple giant-virus orders AND show
  substantial eukaryotic tree placement (see the per-level columns
  `domain`, `phylum`, `class`, `order`, `family` in the detailed
  summary: roughly 50% NCLDV, 45% EUK neighbors). A bin that
  simultaneously lands in two distinct giant-virus orders or carries
  proteins that tree-resolve primarily against eukaryotic references
  (rather than the giant-virus references where viral HGT-origin
  genes should fall) is an indicator of mixed-source sequence, not
  of normal giant-virus HGT content. The `mixed_viral` label makes
  the dominant signal visible without claiming the ML prediction
  distinguishes mixed-virus from cellular-carryover — that
  finer-grained call is still out of scope for the current model.
