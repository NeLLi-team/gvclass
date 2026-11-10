from __future__ import annotations

import csv
import tarfile
from pathlib import Path

import pandas as pd
import pytest

from src.pipeline import summary_writer
from src.pipeline.summary_writer import (
    EXTENDED_SUMMARY_COLUMNS,
    FINAL_SUMMARY_COLUMNS,
    INDIVIDUAL_FINAL_SUMMARY_COLUMNS,
    LEGACY_SUMMARY_HEADERS,
    archive_final_summary_extended_files,
    final_summary_record,
    write_final_summary_extended_files,
    write_final_summary_files,
    write_individual_final_summary_file,
)


def test_failed_queries_are_reported_separately_from_final_summary(
    tmp_path: Path,
) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [
            {
                "query": "completed_query",
                "status": "complete",
                "summary_data": {
                    "query": "completed_query",
                    "taxonomy_majority": "d_NCLDV;p_Nucleocytoviricota",
                    "taxonomy_confidence": "high",
                },
            },
            {
                "query": "failed_query",
                "status": "failed",
                "error": "simulated HMM search failure",
            },
        ],
        tmp_path,
    )

    with open(summary_file, newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert [row["query"] for row in summary_rows] == ["completed_query"]

    with open(tmp_path / "gvclass_failed_queries.tsv", newline="") as handle:
        failed_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert failed_rows == [
        {
            "query": "failed_query",
            "status": "failed",
            "error": "simulated HMM search failure",
        }
    ]
    assert (tmp_path / "gvclass_failed_queries.csv").exists()


def test_species_tree_sentinel_when_no_tree(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [{"query": "q1", "status": "complete", "summary_data": {"query": "q1"}}],
        tmp_path,
    )
    with open(summary_file, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["species_tree_nn_taxonomy"] == "nd"
    assert row["species_tree_nn_genome"] == "nd"
    assert row["species_tree_nn_distance"] == "nd"
    assert row["species_tree_clade_id"] == "nd"


def test_species_tree_taxonomy_value_passthrough(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "species_tree_nn_taxonomy": "d_NCLDV;p_Nucleocytoviricota",
                },
            }
        ],
        tmp_path,
    )
    with open(summary_file, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["species_tree_nn_taxonomy"] == "d_NCLDV;p_Nucleocytoviricota"


def test_extended_summary_written(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_extended_files

    extended = write_final_summary_extended_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "cellular_coherent_contig_count": 0,
                    "contig_attribution_mode": "fna_gene_calling",
                    "mixed_evidence_status": "review_required",
                    "contamination_model_prediction_upper_bound": 58.20048,
                },
            }
        ],
        tmp_path,
    )
    with open(extended, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["contig_attribution_mode"] == "fna_gene_calling"
    assert row["cellular_coherent_contig_count"] == "0"
    assert row["mixed_evidence_status"] == "review_required"
    assert row["contamination_model_prediction_upper_bound"] == "58.200480"


def test_extended_summary_archive_removes_loose_files(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import (
        archive_final_summary_extended_files,
        write_final_summary_extended_files,
    )

    write_final_summary_extended_files(
        [{"query": "q1", "status": "complete", "summary_data": {"query": "q1"}}],
        tmp_path,
    )

    archive = archive_final_summary_extended_files(tmp_path)

    assert archive == tmp_path / "gvclass_summary.extended.tar.gz"
    assert archive.exists()
    assert not (tmp_path / "gvclass_summary.extended.tsv").exists()
    assert not (tmp_path / "gvclass_summary.extended.csv").exists()
    with tarfile.open(archive, "r:gz") as tar_handle:
        assert {
            "gvclass_summary.extended.tsv",
            "gvclass_summary.extended.csv",
        }.issubset(set(tar_handle.getnames()))


def test_quality_context_survives_the_final_summary_record() -> None:
    record = summary_writer.final_summary_record(
        {
            "query": "q1",
            "estimated_completeness": 42.5,
            "completeness_query_support": "limited_support",
            "estimated_completeness_strategy": "novelty_aware_v1",
            "completeness_model_group": "MIRUS",
            "completeness_model_reliability": "unvalidated_legacy",
            "estimated_contamination": 12.5,
            "contamination_query_support": "supported",
            "estimated_contamination_strategy": "trained_model_v1",
            "contamination_model_reliability": "unvalidated_legacy",
        }
    )

    assert list(record) == summary_writer.INDIVIDUAL_FINAL_SUMMARY_COLUMNS
    assert record["estimated_completeness"] == "42.50"
    assert record["completeness_query_support"] == "limited_support"
    assert record["estimated_completeness_strategy"] == "novelty_aware_v1"
    assert record["completeness_model_group"] == "MIRUS"
    assert record["completeness_model_reliability"] == "unvalidated_legacy"
    assert record["contamination_query_support"] == "supported"
    assert record["contamination_model_reliability"] == "unvalidated_legacy"


@pytest.mark.parametrize(
    "writer_name",
    ["write_individual_summary_file", "write_individual_final_summary_file"],
)
def test_required_summary_write_failure_raises(
    tmp_path: Path, writer_name: str
) -> None:
    writer = getattr(summary_writer, writer_name)

    with pytest.raises(FileNotFoundError):
        writer(tmp_path / "missing", "q1", {"query": "q1"})


def test_the_final_summary_file_holds_the_record(tmp_path: Path) -> None:
    summary_data = {"query": "q1", "avgdist": 1.23456, "taxonomy_majority": "d_NCLDV"}

    summary_file = summary_writer.write_individual_final_summary_file(
        tmp_path, "q1", summary_data
    )

    with summary_file.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row == summary_writer.final_summary_record(summary_data)
    assert row["avgdist"] == "1.23"


@pytest.mark.parametrize(
    ("column", "value", "expected"),
    [
        ("avgdist", 1.23456, "1.23"),
        ("estimated_completeness", 88.5, "88.50"),
        ("order_completeness_v2_support_score", 0.12345, "0.12"),
        ("contamination_model_prediction_upper_bound", 0.8765432, "0.876543"),
        ("ncldv_mcp_total", 4.6, "5"),
        ("estimated_completeness", float("nan"), "NaN"),
        ("estimated_completeness", "88.50", "88.50"),
        ("contigs", 3, "3"),
    ],
)
def test_format_value_takes_the_decimals_from_one_table(
    column: str, value: object, expected: str
) -> None:
    assert summary_writer.format_value(column, value) == expected


def test_legacy_summary_reads_the_two_aliased_weighted_completeness_keys(
    tmp_path: Path,
) -> None:
    summary_file = summary_writer.write_individual_summary_file(
        tmp_path,
        "q1",
        {
            "query": "q1",
            "order_weighted_completeness": 55.556,
            "order_weighted_completeness_raw": 66.666,
            "estimated_completeness": 88.5,
        },
    )

    header, values = summary_file.read_text().splitlines()
    row = dict(zip(header.split("\t"), values.split("\t"), strict=True))
    assert header.split("\t") == summary_writer.LEGACY_SUMMARY_HEADERS
    assert row["weighted_order_completeness"] == "55.56"
    assert row["weighted_order_completeness_raw"] == "66.67"
    assert row["estimated_completeness"] == "88.50"


SCREEN_COLUMNS = [
    "screen_status",
    "screen_cellular_pct",
    "screen_giant_virus_pct",
    "screen_flagged_pct",
]


def test_screen_columns_close_the_final_summary_without_moving_the_others() -> None:
    columns = summary_writer.FINAL_SUMMARY_COLUMNS

    assert columns[-4:] == SCREEN_COLUMNS
    assert columns[-5] == "species_tree_clade_id"
    assert columns.index("taxonomy_majority") == 1
    assert columns.index("estimated_completeness") == 14
    assert columns.index("estimated_contamination") == 19


def test_screen_percentages_are_written_with_two_decimals_and_nan_when_missing(
    tmp_path: Path,
) -> None:
    summary_file = summary_writer.write_final_summary_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "screen_status": "partial",
                    "screen_cellular_pct": 12.3456,
                    "screen_giant_virus_pct": float("nan"),
                    "screen_flagged_pct": float("nan"),
                },
            }
        ],
        tmp_path,
    )

    with summary_file.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))

    assert row["screen_status"] == "partial"
    assert row["screen_cellular_pct"] == "12.35"
    assert row["screen_giant_virus_pct"] == "NaN"
    assert row["screen_flagged_pct"] == "NaN"


def test_a_saved_record_has_blank_screen_columns() -> None:
    record = summary_writer.final_summary_record(
        {"query": "q1", "taxonomy_majority": "d_NCLDV"}
    )

    assert record["taxonomy_majority"] == "d_NCLDV"
    assert record["screen_status"] == ""


def test_taxonomy_confidence_column_present_in_summary_schemas() -> None:
    assert "taxonomy_confidence" in FINAL_SUMMARY_COLUMNS
    assert "taxonomy_confidence" in LEGACY_SUMMARY_HEADERS


def test_completeness_reliability_column_in_summary_schema() -> None:
    """The completeness QC column is surfaced as completeness_model_reliability;
    the raw advisory/r2-holdout diagnostics were dropped from the final summary."""

    assert "completeness_model_reliability" in FINAL_SUMMARY_COLUMNS
    for column in (
        "estimated_completeness_quality",
        "estimated_completeness_advisory",
        "estimated_completeness_r2_holdout",
    ):
        assert column not in FINAL_SUMMARY_COLUMNS


def test_summary_writer_preserves_extended_columns_across_resume(
    tmp_path: Path,
) -> None:
    """The final-schema row must keep every summary and extended column.

    <query>.final_summary.tsv and the ``summary`` of the ``results.jsonl`` record
    hold this row; --resume rebuilds both combined tables from the record.
    """
    assert (
        INDIVIDUAL_FINAL_SUMMARY_COLUMNS[: len(FINAL_SUMMARY_COLUMNS)]
        == FINAL_SUMMARY_COLUMNS
    )
    assert set(EXTENDED_SUMMARY_COLUMNS) <= set(INDIVIDUAL_FINAL_SUMMARY_COLUMNS)

    summary_data = {col: "0" for col in INDIVIDUAL_FINAL_SUMMARY_COLUMNS}
    summary_data["query"] = "sample_q"
    summary_data["cellular_coherent_bp_fraction"] = 12.34
    summary_data["contig_attribution_mode"] = "contig_aware"

    write_individual_final_summary_file(tmp_path, "sample_q", summary_data)
    record = final_summary_record(summary_data)
    saved_lines = (tmp_path / "sample_q.final_summary.tsv").read_text().splitlines()
    assert saved_lines[0].split("\t") == INDIVIDUAL_FINAL_SUMMARY_COLUMNS
    assert saved_lines[1].split("\t") == list(record.values())

    resumed_results = [
        {"query": "sample_q", "status": "complete", "summary_data": record}
    ]
    write_final_summary_files(resumed_results, tmp_path)
    write_final_summary_extended_files(resumed_results, tmp_path)
    archive_final_summary_extended_files(tmp_path)

    main_df = pd.read_csv(tmp_path / "gvclass_summary.tsv", sep="\t")
    assert list(main_df.columns) == FINAL_SUMMARY_COLUMNS

    tar_path = tmp_path / "gvclass_summary.extended.tar.gz"
    assert tar_path.exists()
    with tarfile.open(tar_path, "r:gz") as tar:
        member = tar.extractfile("gvclass_summary.extended.tsv")
        assert member is not None
        ext_df = pd.read_csv(member, sep="\t")
    assert list(ext_df.columns) == EXTENDED_SUMMARY_COLUMNS
    assert float(ext_df.loc[0, "cellular_coherent_bp_fraction"]) == 12.34
    assert ext_df.loc[0, "contig_attribution_mode"] == "contig_aware"
