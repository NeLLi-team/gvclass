from __future__ import annotations

import tarfile
from pathlib import Path


def test_output_validator_accepts_archived_query_summaries(tmp_path: Path) -> None:
    from src.utils.output_validator import validate_pipeline_outputs

    work = tmp_path / "q1"
    work.mkdir()
    (work / "q1.summary.tab").write_text("query\ttaxonomy_majority\nq1\ttax\n")
    (work / "q1.final_summary.tsv").write_text("query\ttaxonomy_majority\nq1\ttax\n")
    with tarfile.open(tmp_path / "q1.tar.gz", "w:gz") as tar_handle:
        tar_handle.add(work, arcname="q1")
    (tmp_path / "gvclass_summary.tsv").write_text("query\ttaxonomy_majority\nq1\ttax\n")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["success"] is True
    assert validation["issues"] == []
