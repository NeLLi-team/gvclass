"""Summary file writers for query and pipeline outputs."""

import csv
import io
import math
import os
import tarfile
import traceback
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

LEGACY_SUMMARY_HEADERS: List[str] = [
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

LEGACY_SUMMARY_KEY_MAPPING: Dict[str, str] = {
    "query": "query",
    "taxonomy_majority": "taxonomy_majority",
    "taxonomy_confidence": "taxonomy_confidence",
    "species": "species",
    "genus": "genus",
    "family": "family",
    "order": "order",
    "class": "class",
    "phylum": "phylum",
    "domain": "domain",
    "avgdist": "avgdist",
    "order_dup": "order_dup",
    "order_completeness": "order_completeness",
    "order_completeness_raw": "order_completeness_raw",
    "order_completeness_baseline_mean": "order_completeness_baseline_mean",
    "order_completeness_baseline_std": "order_completeness_baseline_std",
    "order_completeness_reference_order": "order_completeness_reference_order",
    "order_completeness_strategy": "order_completeness_strategy",
    "order_completeness_v2": "order_completeness_v2",
    "order_completeness_v2_strategy": "order_completeness_v2_strategy",
    "order_completeness_v2_strategy2_raw": "order_completeness_v2_strategy2_raw",
    "order_completeness_v2_strategy2_normalized": "order_completeness_v2_strategy2_normalized",
    "order_completeness_v2_support_score": "order_completeness_v2_support_score",
    "order_completeness_v2_ood_flag": "order_completeness_v2_ood_flag",
    "order_completeness_v2_reference_group": "order_completeness_v2_reference_group",
    "order_completeness_v2_validation_mode": "order_completeness_v2_validation_mode",
    "order_completeness_v2_informative_fraction": "order_completeness_v2_informative_fraction",
    "estimated_completeness": "estimated_completeness",
    "estimated_completeness_strategy": "estimated_completeness_strategy",
    "contamination_score_v1": "contamination_score_v1",
    "contamination_flag_v1": "contamination_flag_v1",
    "contamination_source_v1": "contamination_source_v1",
    "contamination_cellular_signal_v1": "contamination_cellular_signal_v1",
    "contamination_phage_signal_v1": "contamination_phage_signal_v1",
    "contamination_duplication_signal_v1": "contamination_duplication_signal_v1",
    "contamination_viral_mixture_signal_v1": "contamination_viral_mixture_signal_v1",
    "contamination_nonviral_hit_fraction_v1": "contamination_nonviral_hit_fraction_v1",
    "estimated_contamination": "estimated_contamination",
    "estimated_contamination_strategy": "estimated_contamination_strategy",
    "suspicious_bp_fraction_v2": "suspicious_bp_fraction_v2",
    "suspicious_contig_count_v2": "suspicious_contig_count_v2",
    "gvog4_completeness": "gvog4_completeness",
    "gvog4_dup": "gvog4_dup",
    "gvog8_completeness": "gvog8_completeness",
    "gvog8_dup": "gvog8_dup",
    "ncldv_mcp_total": "ncldv_mcp_total",
    "vp_completeness": "vp_completeness",
    "vp_mcp": "vp_mcp",
    "plv": "plv",
    "vp_df": "vp_df",
    "mirus_completeness": "mirus_completeness",
    "mirus_df": "mirus_df",
    "mrya_completeness": "mrya_completeness",
    "mrya_dup": "mrya_dup",
    "phage_completeness": "phage_completeness",
    "phage_dup": "phage_dup",
    "busco_completeness": "busco_completeness",
    "busco_dup": "busco_dup",
    "cog_completeness": "cog_completeness",
    "cog_dup": "cog_dup",
    "contigs": "contigs",
    "LENbp": "LENbp",
    "GCperc": "GCperc",
    "genecount": "genecount",
    "CODINGperc": "CODINGperc",
    "ttable": "ttable",
    "order_weighted_completeness": "weighted_order_completeness",
    "order_weighted_completeness_raw": "weighted_order_completeness_raw",
}

FINAL_SUMMARY_COLUMNS: List[str] = [
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
    "completeness_model_reliability",
    "estimated_contamination",
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
]


# Always-on supplementary table archived in gvclass_summary.extended.tar.gz:
# per-contig taxonomic-purity / contamination diagnostics that are 0 for clean
# single-virus genomes and only meaningful under cellular contamination.
EXTENDED_SUMMARY_COLUMNS: List[str] = [
    "query",
    "cellular_coherent_contig_count",
    "cellular_coherent_protein_fraction",
    "cellular_coherent_bp_fraction",
    "cellular_lineage_purity_median",
    "cellular_hit_identity_median",
    "viral_bearing_contig_count",
    "contig_attribution_mode",
]

FAILED_QUERY_COLUMNS: List[str] = ["query", "status", "error"]

TWO_DECIMAL_COLUMNS = {
    "avgdist",
    "order_dup",
    "gvog4_dup",
    "gvog8_dup",
    "busco_dup",
    "cog_dup",
    "mrya_dup",
    "phage_dup",
    "estimated_completeness",
    "estimated_contamination",
    "GCperc",
    "CODINGperc",
    # Per-contig taxonomic-purity features (v1.4.3 Phase 2).
    "cellular_coherent_protein_fraction",
    "cellular_coherent_bp_fraction",
    "cellular_lineage_purity_median",
    "cellular_hit_identity_median",
    # Species-tree placement distance (full precision lives in species_tree_taxonomy.tsv).
    "species_tree_nn_distance",
}

# Quality metrics that must keep higher precision than the default .0f fallback.
FOUR_DECIMAL_COLUMNS: set = set()

LEGACY_TWO_DECIMAL_COLUMNS = {
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
    "contamination_score_v1",
    "contamination_cellular_signal_v1",
    "contamination_phage_signal_v1",
    "contamination_duplication_signal_v1",
    "contamination_viral_mixture_signal_v1",
    "contamination_nonviral_hit_fraction_v1",
    "estimated_contamination",
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
}


def write_individual_summary_file(
    query_output_dir: Path, query_name: str, summary_data: Dict[str, Any], logger
) -> Path:
    """Write the legacy per-query summary file."""
    summary_file = query_output_dir / f"{query_name}.summary.tab"
    try:
        row_data = _build_legacy_summary_row(summary_data)
        with open(summary_file, "w") as handle:
            handle.write("\t".join(LEGACY_SUMMARY_HEADERS) + "\n")
            handle.write("\t".join(row_data) + "\n")
        logger.info(f"Summary file written successfully: {summary_file}")
    except Exception as exc:
        logger.error(f"Failed to write summary file for {query_name}: {exc}")
        logger.error(traceback.format_exc())
    return summary_file


def write_individual_final_summary_file(
    query_output_dir: Path, query_name: str, summary_data: Dict[str, Any], logger
) -> Path:
    """Write a per-query final-schema summary row for resume reconstruction."""
    summary_file = query_output_dir / f"{query_name}.final_summary.tsv"
    try:
        row_data = _build_final_summary_row(
            {
                "query": query_name,
                "status": "complete",
                "summary_data": summary_data,
            }
        )
        with open(summary_file, "w", newline="") as handle:
            handle.write("\t".join(FINAL_SUMMARY_COLUMNS) + "\n")
            handle.write("\t".join(row_data) + "\n")
        logger.info(f"Final-schema summary file written successfully: {summary_file}")
    except Exception as exc:
        logger.error(
            f"Failed to write final-schema summary file for {query_name}: {exc}"
        )
        logger.error(traceback.format_exc())
    return summary_file


def _build_legacy_summary_row(summary_data: Dict[str, Any]) -> List[str]:
    row_data: List[str] = []
    for header in LEGACY_SUMMARY_HEADERS:
        row_data.append(_get_legacy_summary_value(summary_data, header))
    return row_data


def _get_legacy_summary_value(summary_data: Dict[str, Any], header: str) -> str:
    for data_key, header_name in LEGACY_SUMMARY_KEY_MAPPING.items():
        if header_name == header:
            return _format_legacy_summary_value(header, summary_data.get(data_key, ""))
    return _format_legacy_summary_value(header, summary_data.get(header, ""))


def _format_legacy_summary_value(header: str, value: Any) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            # Preserve NaN verbatim so downstream readers can distinguish
            # "skipped / not computed" (e.g. sensitive-mode contamination
            # gate) from a genuine 0.0 score.
            return "NaN"
        if header in LEGACY_TWO_DECIMAL_COLUMNS:
            return f"{value:.2f}"
        return f"{value:.0f}"
    return str(value)


def write_final_summary_files(results: List[Dict[str, Any]], output_dir: Path) -> Path:
    """Write the final TSV and CSV summary tables."""
    summary_tsv = output_dir / "gvclass_summary.tsv"
    summary_csv = output_dir / "gvclass_summary.csv"
    summary_results, failed_results = _partition_summary_results(results)
    _write_failed_query_report(failed_results, output_dir)
    with (
        open(summary_tsv, "w") as tsv_file,
        open(summary_csv, "w", newline="") as csv_file,
    ):
        csv_writer = csv.writer(csv_file)
        tsv_file.write("\t".join(FINAL_SUMMARY_COLUMNS) + "\n")
        csv_writer.writerow(FINAL_SUMMARY_COLUMNS)
        for result in sorted(summary_results, key=lambda item: item["query"]):
            row = _build_final_summary_row(result)
            tsv_file.write("\t".join(row) + "\n")
            csv_writer.writerow(row)
    return summary_tsv


def write_final_summary_extended_files(
    results: List[Dict[str, Any]], output_dir: Path
) -> Path:
    """Write the always-on supplementary diagnostics table.

    gvclass_summary.extended.tsv / .csv carries the per-contig contamination
    features that are 0 for clean single-virus genomes (kept out of the main
    table for readability).
    """
    extended_tsv = output_dir / "gvclass_summary.extended.tsv"
    extended_csv = output_dir / "gvclass_summary.extended.csv"
    summary_results, _ = _partition_summary_results(results)
    with (
        open(extended_tsv, "w") as tsv_file,
        open(extended_csv, "w", newline="") as csv_file,
    ):
        csv_writer = csv.writer(csv_file)
        tsv_file.write("\t".join(EXTENDED_SUMMARY_COLUMNS) + "\n")
        csv_writer.writerow(EXTENDED_SUMMARY_COLUMNS)
        for result in sorted(summary_results, key=lambda item: item["query"]):
            row = [
                _format_final_summary_value(
                    column, result["summary_data"].get(column, "")
                )
                for column in EXTENDED_SUMMARY_COLUMNS
            ]
            tsv_file.write("\t".join(row) + "\n")
            csv_writer.writerow(row)
    return extended_tsv


def archive_final_summary_extended_files(output_dir: Path) -> Optional[Path]:
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
    with open(tmp_path, "rb") as handle:
        os.fsync(handle.fileno())
    if not tarfile.is_tarfile(tmp_path):
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Archive verification failed for {tmp_path}")
    os.replace(tmp_path, archive_path)
    for path in existing_paths:
        path.unlink(missing_ok=True)
    return archive_path


def _partition_summary_results(
    results: List[Dict[str, Any]],
) -> tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
    summary_results = []
    failed_results = []
    for result in results:
        if result.get("status") == "complete" and "summary_data" in result:
            summary_results.append(result)
            continue
        failed_results.append(_build_failed_query_row(result))
    return summary_results, failed_results


def _build_failed_query_row(result: Dict[str, Any]) -> Dict[str, str]:
    status = result.get("status") or "missing_summary"
    error = result.get("error") or ""
    if status == "complete" and "summary_data" not in result:
        error = "complete result missing summary_data"
    return {
        "query": str(result.get("query", "")),
        "status": str(status),
        "error": str(error),
    }


def _write_failed_query_report(
    failed_results: List[Dict[str, str]], output_dir: Path
) -> None:
    failed_tsv = output_dir / "gvclass_failed_queries.tsv"
    failed_csv = output_dir / "gvclass_failed_queries.csv"
    if not failed_results:
        failed_tsv.unlink(missing_ok=True)
        failed_csv.unlink(missing_ok=True)
        return

    with (
        open(failed_tsv, "w", newline="") as tsv_handle,
        open(failed_csv, "w", newline="") as csv_handle,
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


def read_individual_summary_results(
    output_dir: Path, query_names: Optional[Iterable[str]] = None
) -> List[Dict[str, Any]]:
    """Read existing per-query summary files into final-summary result records."""
    wanted_names = (
        set(query_names)
        if query_names is not None
        else _discover_individual_summary_names(output_dir)
    )
    results: List[Dict[str, Any]] = []
    for query_name in sorted(wanted_names):
        result = _read_best_individual_summary_result(output_dir, query_name)
        if result is not None:
            results.append(result)
    return results


def _discover_individual_summary_names(output_dir: Path) -> set[str]:
    names = {
        path.name.removesuffix(".summary.tab")
        for path in output_dir.glob("*.summary.tab")
    }
    names.update(
        path.name.removesuffix(".final_summary.tsv")
        for path in output_dir.glob("*.final_summary.tsv")
    )
    names.update(
        path.name.removesuffix(".tar.gz")
        for path in output_dir.glob("*.tar.gz")
        if not path.name.startswith("gvclass_summary.")
    )
    return names


def _read_best_individual_summary_result(
    output_dir: Path, query_name: str
) -> Optional[Dict[str, Any]]:
    final_summary = output_dir / f"{query_name}.final_summary.tsv"
    if final_summary.exists():
        return _read_individual_summary_result(final_summary, query_name)
    legacy_summary = output_dir / f"{query_name}.summary.tab"
    if legacy_summary.exists():
        return _read_individual_summary_result(legacy_summary, query_name)
    for member_name in (
        f"{query_name}/{query_name}.final_summary.tsv",
        f"{query_name}/{query_name}.summary.tab",
    ):
        result = _read_individual_summary_result_from_archive(
            output_dir / f"{query_name}.tar.gz",
            member_name,
            query_name,
        )
        if result is not None:
            return result
    return None


def _read_individual_summary_result(
    summary_file: Path, query_name: str
) -> Optional[Dict[str, Any]]:
    with open(summary_file, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
    if row is None:
        return None
    summary_data = {key: value for key, value in row.items() if key is not None}
    summary_data["query"] = query_name
    return {"query": query_name, "status": "complete", "summary_data": summary_data}


def _read_individual_summary_result_from_archive(
    archive_path: Path,
    member_name: str,
    query_name: str,
) -> Optional[Dict[str, Any]]:
    if not archive_path.exists():
        return None
    try:
        with tarfile.open(archive_path, "r:gz") as tar_handle:
            member = tar_handle.getmember(member_name)
            member_file = tar_handle.extractfile(member)
            if member_file is None:
                return None
            text = member_file.read().decode()
    except (KeyError, OSError, tarfile.TarError, UnicodeDecodeError):
        return None

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    row = next(reader, None)
    if row is None:
        return None
    summary_data = {key: value for key, value in row.items() if key is not None}
    summary_data["query"] = query_name
    return {"query": query_name, "status": "complete", "summary_data": summary_data}


# Species-tree placement columns. When a run builds no species tree these are
# absent from summary_data; emit the "nd" (not determined) sentinel rather than a
# blank so all four read consistently. ``_SPECIES_TREE_ND_VALUES`` also maps the
# retired ``no-species-tree-calculated`` literal to ``nd`` so a ``--resume`` over
# per-query files written by an older gvclass version is migrated, not preserved.
_SPECIES_TREE_ND_COLUMNS = {
    "species_tree_nn_taxonomy",
    "species_tree_nn_genome",
    "species_tree_nn_distance",
    "species_tree_clade_id",
}
_SPECIES_TREE_ND_VALUES = ("", None, "no-species-tree-calculated")


def _build_final_summary_row(result: Dict[str, Any]) -> List[str]:
    if result["status"] != "complete" or "summary_data" not in result:
        return [result["query"]] + [""] * (len(FINAL_SUMMARY_COLUMNS) - 1)

    summary_data = result["summary_data"]
    row: List[str] = []
    for column in FINAL_SUMMARY_COLUMNS:
        value = summary_data.get(column, "")
        if column in _SPECIES_TREE_ND_COLUMNS and value in _SPECIES_TREE_ND_VALUES:
            value = "nd"
        row.append(_format_final_summary_value(column, value))
    return row


def _format_final_summary_value(column: str, value: Any) -> str:
    if not isinstance(value, float):
        return str(value)
    if math.isnan(value):
        # See _format_legacy_summary_value: preserve NaN rather than collapse
        # to 0.00, so the sensitive-mode contamination skip is visible.
        return "NaN"
    if column in FOUR_DECIMAL_COLUMNS:
        return f"{value:.4f}"
    if column in TWO_DECIMAL_COLUMNS:
        return f"{value:.2f}"
    return f"{value:.0f}"
