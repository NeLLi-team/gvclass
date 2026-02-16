"""Summary file writers for query and pipeline outputs."""

from pathlib import Path
from typing import Any, Dict, List
import csv
import traceback


LEGACY_SUMMARY_HEADERS: List[str] = [
    "query",
    "taxonomy_majority",
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
    "gvog4_unique",
    "gvog8_unique",
    "gvog8_total",
    "gvog8_dup",
    "ncldv_mcp_total",
    "mcp_total",
    "vp_completeness",
    "vp_mcp",
    "plv",
    "vp_df",
    "mirus_completeness",
    "mirus_df",
    "mrya_unique",
    "mrya_total",
    "phage_unique",
    "phage_total",
    "cellular_unique",
    "cellular_total",
    "cellular_dup",
    "contigs",
    "LENbp",
    "GCperc",
    "genecount",
    "CODINGperc",
    "ttable",
    "weighted_order_completeness",
]

LEGACY_SUMMARY_KEY_MAPPING: Dict[str, str] = {
    "query": "query",
    "taxonomy_majority": "taxonomy_majority",
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
    "gvog4_unique": "gvog4_unique",
    "gvog8_unique": "gvog8_unique",
    "gvog8_total": "gvog8_total",
    "gvog8_dup": "gvog8_dup",
    "ncldv_mcp_total": "ncldv_mcp_total",
    "mcp_total": "mcp_total",
    "vp_completeness": "vp_completeness",
    "vp_mcp": "vp_mcp",
    "plv": "plv",
    "vp_df": "vp_df",
    "mirus_completeness": "mirus_completeness",
    "mirus_df": "mirus_df",
    "mrya_unique": "mrya_unique",
    "mrya_total": "mrya_total",
    "phage_unique": "phage_unique",
    "phage_total": "phage_total",
    "cellular_unique": "cellular_unique",
    "cellular_total": "cellular_total",
    "cellular_dup": "cellular_dup",
    "contigs": "contigs",
    "LENbp": "LENbp",
    "GCperc": "GCperc",
    "genecount": "genecount",
    "CODINGperc": "CODINGperc",
    "ttable": "ttable",
    "order_weighted_completeness": "weighted_order_completeness",
}

FINAL_SUMMARY_COLUMNS: List[str] = [
    "query",
    "taxonomy_majority",
    "taxonomy_strict",
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
    "order_weighted_completeness",
    "order_confidence_score",
    "gvog4_unique",
    "gvog8_unique",
    "gvog8_total",
    "gvog8_dup",
    "ncldv_mcp_total",
    "mcp_total",
    "vp_completeness",
    "vp_mcp",
    "plv",
    "vp_df",
    "mirus_completeness",
    "mirus_df",
    "mrya_unique",
    "mrya_total",
    "phage_unique",
    "phage_total",
    "cellular_unique",
    "cellular_total",
    "cellular_dup",
    "contigs",
    "LENbp",
    "GCperc",
    "genecount",
    "CODINGperc",
    "ttable",
]

TWO_DECIMAL_COLUMNS = {
    "avgdist",
    "order_dup",
    "gvog8_dup",
    "vp_df",
    "mirus_df",
    "cellular_dup",
    "order_completeness",
    "order_weighted_completeness",
    "order_confidence_score",
    "GCperc",
    "CODINGperc",
}

LEGACY_TWO_DECIMAL_COLUMNS = {
    "avgdist",
    "order_dup",
    "order_completeness",
    "gvog8_dup",
    "vp_df",
    "mirus_df",
    "cellular_dup",
    "GCperc",
    "CODINGperc",
    "weighted_order_completeness",
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
        if header in LEGACY_TWO_DECIMAL_COLUMNS:
            return f"{value:.2f}"
        return f"{value:.0f}"
    return str(value)


def write_final_summary_files(results: List[Dict[str, Any]], output_dir: Path) -> Path:
    """Write the final TSV and CSV summary tables."""
    summary_tsv = output_dir / "gvclass_summary.tsv"
    summary_csv = output_dir / "gvclass_summary.csv"
    with open(summary_tsv, "w") as tsv_file, open(
        summary_csv, "w", newline=""
    ) as csv_file:
        csv_writer = csv.writer(csv_file)
        tsv_file.write("\t".join(FINAL_SUMMARY_COLUMNS) + "\n")
        csv_writer.writerow(FINAL_SUMMARY_COLUMNS)
        for result in sorted(results, key=lambda item: item["query"]):
            row = _build_final_summary_row(result)
            tsv_file.write("\t".join(row) + "\n")
            csv_writer.writerow(row)
    return summary_tsv


def _build_final_summary_row(result: Dict[str, Any]) -> List[str]:
    if result["status"] != "complete" or "summary_data" not in result:
        return [result["query"]] + [""] * (len(FINAL_SUMMARY_COLUMNS) - 1)

    summary_data = result["summary_data"]
    return [
        _format_final_summary_value(column, summary_data.get(column, ""))
        for column in FINAL_SUMMARY_COLUMNS
    ]


def _format_final_summary_value(column: str, value: Any) -> str:
    if not isinstance(value, float):
        return str(value)
    if column in TWO_DECIMAL_COLUMNS:
        return f"{value:.2f}"
    return f"{value:.0f}"
