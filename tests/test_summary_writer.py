from __future__ import annotations

import csv
from pathlib import Path


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


def test_r2_holdout_formats_four_decimals() -> None:
    from src.pipeline.summary_writer import _format_final_summary_value

    assert (
        _format_final_summary_value("estimated_completeness_r2_holdout", 0.8734)
        == "0.8734"
    )
    assert (
        _format_final_summary_value("estimated_completeness_r2_holdout", 0.0)
        == "0.0000"
    )


def test_r2_holdout_empty_and_string_passthrough() -> None:
    from src.pipeline.summary_writer import _format_final_summary_value

    assert _format_final_summary_value("estimated_completeness_r2_holdout", "") == ""
    # Resume reads the value back as a string; must be idempotent.
    assert (
        _format_final_summary_value("estimated_completeness_r2_holdout", "0.8734")
        == "0.8734"
    )


def test_write_final_summary_files_r2_end_to_end(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "estimated_completeness_r2_holdout": 0.8734,
                },
            }
        ],
        tmp_path,
    )

    with open(summary_file, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))

    # Must NOT be truncated to "1" by the .0f fallback.
    assert row["estimated_completeness_r2_holdout"] == "0.8734"
