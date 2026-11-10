"""Summary file writers for query and pipeline outputs.

A query's summary row is formatted once, by ``final_summary_record``. The
per-query ``final_summary.tsv`` in the archive, the ``summary`` of its
``results.jsonl`` record, and the combined tables all use that formatting, so a
``--resume`` reproduces the combined tables from the records without a second
code path.
"""

import csv
import math
import os
import tarfile
from pathlib import Path
from typing import Any

LEGACY_SUMMARY_HEADERS: list[str] = [
    "query",
    "taxonomy_majority",
    "taxonomy_confidence",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "domain",
    "avgdist",
    "order_dup",
    "order_completeness",
    "order_completeness_raw",
    "order_completeness_baseline_mean",
    "order_completeness_baseline_std",
    "order_completeness_reference_order",
    "order_completeness_strategy",
    "order_completeness_v2",
    "order_completeness_v2_strategy",
    "order_completeness_v2_strategy2_raw",
    "order_completeness_v2_strategy2_normalized",
    "order_completeness_v2_support_score",
    "order_completeness_v2_ood_flag",
    "order_completeness_v2_reference_group",
    "order_completeness_v2_validation_mode",
    "order_completeness_v2_informative_fraction",
    "estimated_completeness",
    "estimated_completeness_strategy",
    "completeness_query_support",
    "completeness_model_group",
    "completeness_model_reliability",
    "contamination_score_v1",
    "contamination_flag_v1",
    "contamination_source_v1",
    "contamination_cellular_signal_v1",
    "contamination_phage_signal_v1",
    "contamination_duplication_signal_v1",
    "contamination_viral_mixture_signal_v1",
    "contamination_nonviral_hit_fraction_v1",
    "estimated_contamination",
    "estimated_contamination_strategy",
    "contamination_query_support",
    "contamination_model_reliability",
    "contamination_type",
    "suspicious_bp_fraction_v2",
    "suspicious_contig_count_v2",
    "gvog4_completeness",
    "gvog4_dup",
    "gvog8_completeness",
    "gvog8_dup",
    "ncldv_mcp_total",
    "vp_completeness",
    "vp_mcp",
    "plv",
    "vp_df",
    "mirus_completeness",
    "mirus_df",
    "mrya_completeness",
    "mrya_dup",
    "phage_completeness",
    "phage_dup",
    "busco_completeness",
    "busco_dup",
    "cog_completeness",
    "cog_dup",
    "contigs",
    "LENbp",
    "GCperc",
    "genecount",
    "CODINGperc",
    "ttable",
    "weighted_order_completeness",
    "weighted_order_completeness_raw",
]

FINAL_SUMMARY_COLUMNS: list[str] = [
    "query",
    "taxonomy_majority",
    # Concatenated-marker species-tree placement — a strong phylogenomic signal,
    # so it sits up front next to the single-gene-consensus taxonomy_majority.
    # Value is "nd" (not determined) when the run built no species tree.
    "species_tree_nn_taxonomy",
    "taxonomy_confidence",
    # Unified capsid-type tally across all MCP panels (NCLDV/Mirus/PPV caps
    # groups), e.g. "Nucleocytoviricota:4,Gossevirus:1".
    "capsid_group",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "domain",
    "avgdist",
    "order_dup",
    "estimated_completeness",
    "completeness_query_support",
    "estimated_completeness_strategy",
    "completeness_model_group",
    "completeness_model_reliability",
    "estimated_contamination",
    "contamination_query_support",
    "estimated_contamination_strategy",
    "contamination_model_reliability",
    "contamination_type",
    "gvog4_completeness",
    "gvog4_dup",
    "gvog8_completeness",
    "gvog8_dup",
    "busco_completeness",
    "busco_dup",
    "cog_completeness",
    "cog_dup",
    "mrya_completeness",
    "mrya_dup",
    "phage_completeness",
    "phage_dup",
    "ncldv_mcp_total",
    "vp_completeness",
    "vp_mcp",
    "plv",
    "mirus_completeness",
    "contigs",
    "LENbp",
    "GCperc",
    "genecount",
    "CODINGperc",
    "ttable",
    # --species-tree placement detail ("nd" unless a species tree was built).
    "species_tree_nn_genome",
    "species_tree_nn_distance",
    "species_tree_clade_id",
    # Contig screen (src/pipeline/contig_screen_step.py), filled after the queries
    # finish; the per-query summary files are written before it runs.
    "screen_status",
    "screen_cellular_pct",
    "screen_giant_virus_pct",
    "screen_flagged_pct",
]


# Always-on supplementary table archived in gvclass_summary.extended.tar.gz:
# per-contig taxonomic-purity / contamination diagnostics that are 0 for clean
# single-virus genomes and only meaningful under cellular contamination.
EXTENDED_SUMMARY_COLUMNS: list[str] = [
    "query",
    "cellular_coherent_contig_count",
    "cellular_coherent_protein_fraction",
    "cellular_coherent_bp_fraction",
    "cellular_lineage_purity_median",
    "cellular_hit_identity_median",
    "viral_bearing_contig_count",
    "mixed_evidence_contig_count",
    "ambiguous_eve_contig_count",
    "ambiguous_eve_protein_count",
    "mixed_evidence_status",
    "ambiguous_eve_evidence_status",
    "contig_length_basis",
    "protein_estimated_contig_count",
    "contig_attribution_mode",
    "contamination_model_applicability",
    "contamination_feature_version",
    "contamination_feature_contract",
    "contamination_prediction_range_status",
    "contamination_model_prediction_upper_bound",
]

FAILED_QUERY_COLUMNS: list[str] = ["query", "status", "error"]


#: Decimals of a floating-point column. A float in any other column has none.
DECIMALS: dict[str, int] = {
    **dict.fromkeys(
        (
            "avgdist",
            "order_dup",
            "order_completeness",
            "order_completeness_raw",
            "order_completeness_baseline_mean",
            "order_completeness_baseline_std",
            "order_completeness_v2",
            "order_completeness_v2_strategy2_raw",
            "order_completeness_v2_strategy2_normalized",
            "order_completeness_v2_support_score",
            "order_completeness_v2_informative_fraction",
            "estimated_completeness",
            "estimated_contamination",
            "contamination_score_v1",
            "contamination_cellular_signal_v1",
            "contamination_phage_signal_v1",
            "contamination_duplication_signal_v1",
            "contamination_viral_mixture_signal_v1",
            "contamination_nonviral_hit_fraction_v1",
            "suspicious_bp_fraction_v2",
            "gvog4_dup",
            "gvog8_dup",
            "busco_dup",
            "cog_dup",
            "mrya_dup",
            "phage_dup",
            "vp_df",
            "mirus_df",
            "GCperc",
            "CODINGperc",
            "weighted_order_completeness",
            "weighted_order_completeness_raw",
            "cellular_coherent_protein_fraction",
            "cellular_coherent_bp_fraction",
            "cellular_lineage_purity_median",
            "cellular_hit_identity_median",
            # Full precision lives in species_tree_taxonomy.tsv.
            "species_tree_nn_distance",
            "screen_cellular_pct",
            "screen_giant_virus_pct",
            "screen_flagged_pct",
        ),
        2,
    ),
    "contamination_model_prediction_upper_bound": 6,
}

#: Legacy ``summary.tab`` headers whose value sits under another summary key.
_LEGACY_ALIASES: dict[str, str] = {
    "weighted_order_completeness": "order_weighted_completeness",
    "weighted_order_completeness_raw": "order_weighted_completeness_raw",
}

#: Species-tree placement columns. A run without a species tree has no value for
#: them; the table says ``nd`` (not determined) rather than leaving a blank.
_SPECIES_TREE_ND_COLUMNS = {
    "species_tree_nn_taxonomy",
    "species_tree_nn_genome",
    "species_tree_nn_distance",
    "species_tree_clade_id",
}

INDIVIDUAL_FINAL_SUMMARY_COLUMNS: list[str] = FINAL_SUMMARY_COLUMNS + [
    col for col in EXTENDED_SUMMARY_COLUMNS if col not in FINAL_SUMMARY_COLUMNS
]


def write_individual_summary_file(
    query_output_dir: Path, query_name: str, summary_data: dict[str, Any]
) -> Path:
    """Write the legacy per-query ``<query>.summary.tab``."""
    summary_file = query_output_dir / f"{query_name}.summary.tab"
    row = [
        format_value(header, summary_data.get(_LEGACY_ALIASES.get(header, header), ""))
        for header in LEGACY_SUMMARY_HEADERS
    ]
    with summary_file.open("w") as handle:
        handle.write("\t".join(LEGACY_SUMMARY_HEADERS) + "\n")
        handle.write("\t".join(row) + "\n")
    return summary_file


def final_summary_record(summary_data: dict[str, Any]) -> dict[str, str]:
    """Return the formatted final-schema row of a query, extended columns included."""
    return dict(
        zip(
            INDIVIDUAL_FINAL_SUMMARY_COLUMNS,
            _format_row(summary_data, INDIVIDUAL_FINAL_SUMMARY_COLUMNS),
            strict=True,
        )
    )


def write_individual_final_summary_file(
    query_output_dir: Path, query_name: str, summary_data: dict[str, Any]
) -> Path:
    """Write ``<query>.final_summary.tsv``: the row that the run record keeps."""
    summary_file = query_output_dir / f"{query_name}.final_summary.tsv"
    record = final_summary_record(summary_data)
    with summary_file.open("w", newline="") as handle:
        handle.write("\t".join(record) + "\n")
        handle.write("\t".join(record.values()) + "\n")
    return summary_file


def write_final_summary_files(results: list[dict[str, Any]], output_dir: Path) -> Path:
    """Write the final TSV and CSV summary tables."""
    summary_tsv = output_dir / "gvclass_summary.tsv"
    summary_csv = output_dir / "gvclass_summary.csv"
    summary_results, failed_results = _partition_summary_results(results)
    _write_failed_query_report(failed_results, output_dir)
    _write_tables(summary_results, FINAL_SUMMARY_COLUMNS, summary_tsv, summary_csv)
    return summary_tsv


def write_final_summary_extended_files(
    results: list[dict[str, Any]], output_dir: Path
) -> Path:
    """Write the always-on supplementary diagnostics table.

    gvclass_summary.extended.tsv / .csv carries the per-contig contamination
    features that are 0 for clean single-virus genomes (kept out of the main
    table for readability).
    """
    extended_tsv = output_dir / "gvclass_summary.extended.tsv"
    extended_csv = output_dir / "gvclass_summary.extended.csv"
    summary_results, _ = _partition_summary_results(results)
    _write_tables(summary_results, EXTENDED_SUMMARY_COLUMNS, extended_tsv, extended_csv)
    return extended_tsv


def _write_tables(
    results: list[dict[str, Any]], columns: list[str], tsv_path: Path, csv_path: Path
) -> None:
    """Write one row per result, sorted by query, as a TSV and a CSV table."""
    with tsv_path.open("w") as tsv_file, csv_path.open("w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        tsv_file.write("\t".join(columns) + "\n")
        csv_writer.writerow(columns)
        for result in sorted(results, key=lambda item: item["query"]):
            row = _format_row(result["summary_data"], columns)
            tsv_file.write("\t".join(row) + "\n")
            csv_writer.writerow(row)


def archive_final_summary_extended_files(output_dir: Path) -> Path | None:
    """Archive extended combined tables and remove the loose copies."""
    extended_paths = [
        output_dir / "gvclass_summary.extended.tsv",
        output_dir / "gvclass_summary.extended.csv",
    ]
    existing_paths = [path for path in extended_paths if path.exists()]
    if not existing_paths:
        return None

    archive_path = output_dir / "gvclass_summary.extended.tar.gz"
    tmp_path = output_dir / "gvclass_summary.extended.tar.gz.part"
    tmp_path.unlink(missing_ok=True)
    with tarfile.open(tmp_path, "w:gz") as tar_handle:
        for path in existing_paths:
            tar_handle.add(path, arcname=path.name)
    with tmp_path.open("rb") as handle:
        os.fsync(handle.fileno())
    if not tarfile.is_tarfile(tmp_path):
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Archive verification failed for {tmp_path}")
    tmp_path.replace(archive_path)
    for path in existing_paths:
        path.unlink(missing_ok=True)
    return archive_path


def _partition_summary_results(
    results: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Split results into those with a summary row and failed-query report rows."""
    summary_results = [result for result in results if result["status"] == "complete"]
    failed_results = [
        {
            "query": result["query"],
            "status": result["status"],
            "error": result.get("error") or "",
        }
        for result in results
        if result["status"] != "complete"
    ]
    return summary_results, failed_results


def _write_failed_query_report(
    failed_results: list[dict[str, str]], output_dir: Path
) -> None:
    failed_tsv = output_dir / "gvclass_failed_queries.tsv"
    failed_csv = output_dir / "gvclass_failed_queries.csv"
    if not failed_results:
        failed_tsv.unlink(missing_ok=True)
        failed_csv.unlink(missing_ok=True)
        return

    with (
        failed_tsv.open("w", newline="") as tsv_handle,
        failed_csv.open("w", newline="") as csv_handle,
    ):
        tsv_writer = csv.DictWriter(
            tsv_handle, fieldnames=FAILED_QUERY_COLUMNS, delimiter="\t"
        )
        csv_writer = csv.DictWriter(csv_handle, fieldnames=FAILED_QUERY_COLUMNS)
        tsv_writer.writeheader()
        csv_writer.writeheader()
        for row in sorted(failed_results, key=lambda item: item["query"]):
            tsv_writer.writerow(row)
            csv_writer.writerow(row)


def _format_row(summary_data: dict[str, Any], columns: list[str]) -> list[str]:
    row = []
    for column in columns:
        value = summary_data.get(column, "")
        if column in _SPECIES_TREE_ND_COLUMNS and value in ("", None):
            value = "nd"
        row.append(format_value(column, value))
    return row


def format_value(column: str, value: object) -> str:
    """Return ``value`` as table text; a float gets the decimals of its column.

    NaN stays ``NaN`` so a reader can tell "not computed" (for example the
    sensitive-mode contamination gate) from a genuine 0.0.
    """
    if not isinstance(value, float):
        return str(value)
    if math.isnan(value):
        return "NaN"
    return f"{value:.{DECIMALS.get(column, 0)}f}"
