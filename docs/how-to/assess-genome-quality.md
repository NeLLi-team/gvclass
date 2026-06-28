# Assess genome quality

Open `gvclass_summary.tsv` (or the `.csv`) after a run and read the columns below to decide whether a classified giant virus MAG (GVMAG) is high quality or needs manual curation. You need a finished classification first; see [Classify bins](classify-bins.md). Full column definitions live in [the output reference](../reference/output.md).

## Read completeness

`estimated_completeness` is the percentage of the expected genome recovered for the assigned lineage. A value near 100 means little is missing. In the bundled example, `PkV-RF01` and `AC3300027503___Ga0255182_1000024` both score `100.00`, while the 10-contig `GVMAG-S-1096109-37` scores `82.86`.

Read `completeness_model_reliability` next to it. The value (`advisory_only`, `moderate`, or `high`) describes the per-order model that produced the estimate, not your genome. A `100.00` tagged `advisory_only` means the assigned order has too few references to calibrate the model, so treat that number as a rough guide. The `GVMAG-S-1096109-37` estimate is tagged `moderate`, which carries more weight.

!!! note
    `estimated_completeness` is the only completeness field in the main table. Switch the underlying estimator with `--completeness-mode legacy` for a fixed-panel calculation instead (see [the CLI reference](../reference/cli.md)).

## Read contamination

`estimated_contamination` is the trained-model estimate. On clean bins it sits near zero; the three example genomes report `0.00`, `0.08`, and `0.00`.

`contamination_type` names the likely source. It reads `clean` below the threshold and carries a specific label once `estimated_contamination` reaches 10 or higher.

| `contamination_type` | likely source |
| --- | --- |
| `clean` | contamination below the threshold; nothing to remove |
| `cellular` | eukaryotic or prokaryotic cellular sequence |
| `mixed_viral` | two or more viral orders mixed in one bin |
| `phage` | bacteriophage sequence |
| `duplication` | duplicated content, often an assembly chimera |
| `uncertain` | ambiguous signal flagged for triage (can be a novel-virus signature, not true contamination) |

## Read duplication

Duplication catches bins that merge multiple populations or chimeric assemblies. Check `order_dup` and `gvog8_dup`:

- above ~2: multiple populations or assembly chimeras are likely; inspect the bin
- below ~1.5: typically clean

Each `{panel}_dup` column is a duplication factor, total marker hits divided by distinct markers present. A single-copy marker found twice pushes the factor above 1.

## Check for cellular carry-over

Giant viruses obligately lack the conserved cellular markers that cellular genomes carry, so any of those markers signals sequence from a host or a co-binned cell. Watch two panels and their duplication factors:

- `busco_completeness` (`n/255`) with `busco_dup`: eukaryotic BUSCO markers
- `cog_completeness` (`n/56`) with `cog_dup`: universal cellular COG (UNI56) markers

A non-zero `busco_completeness` or `cog_completeness`, especially with an elevated `_dup`, points to cellular carry-over. See [the markers reference](../reference/markers.md) for every panel.

!!! warning
    Giant viruses carry many eukaryote-like genes acquired by horizontal gene transfer. GVClass does NOT count those as contamination; the cellular signal is restricted to conserved cellular HMMs (BUSCO eukaryotic plus UNI56 universal cellular) that giant viruses lack. `mixed_viral` likewise means several viral orders mixed in one bin, not host genes. See [quality metrics](../explanation/quality-metrics.md) for the reasoning.

## Look deeper

Two extra outputs explain a flagged genome.

- `gvclass_summary.extended.tsv` is always written. It holds the per-contig diagnostics (`cellular_coherent_*`, `cellular_lineage_purity_median`, `viral_bearing_contig_count`, `contig_attribution_mode`) that the main table omits to stay readable.
- `<query>.contamination_candidates.tsv` is written per query when `estimated_contamination` reaches 10, the type is interpretable, and suspicious contigs are found. It names the specific contigs to drop or check.

## The verdict

Apply one combined rule:

- high `estimated_completeness`, low duplication (`order_dup` and `gvog8_dup` below ~1.5), and `contamination_type` `clean`: a high-quality GVMAG
- low completeness, high duplication, or any non-clean type: send the bin for manual curation

!!! tip
    All three bundled example genomes pass: high completeness, low duplication, and `clean`. Use them as a reference point when reading your own [summary table](../reference/output.md).
