from __future__ import annotations

import json
import logging
import tarfile
from pathlib import Path

import pytest

from src.utils.output_validator import validate_pipeline_outputs


@pytest.mark.parametrize(
    "shared_directory", [None, "contig_screen", "species_tree", ".species_tree_work"]
)
def test_output_validator_accepts_archived_query_summaries(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, shared_directory: str | None
) -> None:
    """Shared output is not a query with missing summary/archive files."""
    work = tmp_path / "q1"
    work.mkdir()
    (work / "q1.summary.tab").write_text("query\ttaxonomy_majority\nq1\ttax\n")
    (work / "q1.final_summary.tsv").write_text("query\ttaxonomy_majority\nq1\ttax\n")
    with tarfile.open(tmp_path / "q1.tar.gz", "w:gz") as tar_handle:
        tar_handle.add(work, arcname="q1")
    (tmp_path / "gvclass_summary.tsv").write_text("query\ttaxonomy_majority\nq1\ttax\n")
    (tmp_path / "run_status.json").write_text(
        json.dumps({"queries": {"q1": {}}}), encoding="utf-8"
    )
    if shared_directory:
        (tmp_path / shared_directory / "q1").mkdir(parents=True)

    with caplog.at_level(logging.DEBUG, logger="src.utils.output_validator"):
        validation = validate_pipeline_outputs(tmp_path, verbose=True)

    assert validation["success"] is True
    assert validation["issues"] == []
    assert validation["warnings"] == []
    assert "Validation successful" in caplog.messages
    assert all(symbol not in caplog.text for symbol in ("✅", "❌", "⚠", "✓", "\x1b"))


@pytest.mark.parametrize(
    "query", ["contig_screen", "species_tree", ".species_tree_work"]
)
@pytest.mark.parametrize("metadata", ["known-query", "absent", "malformed", "nonutf8"])
def test_shared_folder_names_do_not_hide_query_warnings(
    tmp_path: Path, query: str, metadata: str
) -> None:
    """Known queries and uncertain folders retain missing-file diagnostics."""
    work = tmp_path / "q1"
    work.mkdir()
    for name in ("q1.summary.tab", "q1.final_summary.tsv"):
        (work / name).write_text(
            "query\ttaxonomy_majority\nq1\ttax\n", encoding="utf-8"
        )
    with tarfile.open(tmp_path / "q1.tar.gz", "w:gz") as archive:
        archive.add(work, arcname="q1")
    (tmp_path / query).mkdir()
    if metadata == "known-query":
        (tmp_path / "run_status.json").write_text(
            json.dumps({"queries": {"q1": {}, query: {}}}), encoding="utf-8"
        )
    elif metadata == "malformed":
        (tmp_path / "run_status.json").write_text("broken{", encoding="utf-8")
    elif metadata == "nonutf8":
        (tmp_path / "run_status.json").write_bytes(b"\xff")

    validation = validate_pipeline_outputs(tmp_path, verbose=True)

    assert validation["warnings"] == [
        f"Expected file not found for {query}: {query}.summary.tab",
        f"Expected file not found for {query}: {query}.tar.gz",
    ]


@pytest.mark.parametrize("valid_summary", [False, True])
def test_output_validator_diagnostics_are_plain(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, valid_summary: bool
) -> None:
    """Missing archives and invalid summaries retain plain, useful diagnostics."""
    summary = "query\ttaxonomy_majority\nq1\ttax\n"
    (tmp_path / "q1.summary.tab").write_text(
        summary if valid_summary else "wrong\nq1\n", encoding="utf-8"
    )
    (tmp_path / "gvclass_summary.tsv").write_text(summary, encoding="utf-8")

    with caplog.at_level(logging.DEBUG, logger="src.utils.output_validator"):
        validation = validate_pipeline_outputs(tmp_path, verbose=True)

    assert validation["success"] is valid_summary
    assert len(validation["warnings"]) == 1
    assert "No compressed archives" in validation["warnings"][0]
    assert "1 warning(s) found" in caplog.messages
    if not valid_summary:
        assert validation["issues"]
        assert any("Validation found" in message for message in caplog.messages)
    assert all(symbol not in caplog.text for symbol in ("✅", "❌", "⚠", "✓", "\x1b"))
