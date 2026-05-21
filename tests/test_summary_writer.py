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
