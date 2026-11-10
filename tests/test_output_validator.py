"""Output validation reads ``results.jsonl`` and the archives that it names."""

from __future__ import annotations

import hashlib
import io
import json
import logging
import tarfile
from dataclasses import asdict
from pathlib import Path

import pytest

from src.pipeline.run_status import (
    RESULTS_FILENAME,
    RUN_STATUS_FILENAME,
    SOFTWARE_VERSION,
    QueryRecord,
)
from src.utils.output_validator import validate_pipeline_outputs

GOOD_TABLE = "query\ttaxonomy_majority\nq1\ttax\n"
LOGGER_NAME = "src.utils.output_validator"


def write_archive(output: Path, name: str, members: dict[str, str]) -> Path:
    """Write ``<name>.tar.gz`` with the given members, named under ``<name>/``."""
    archive = output / f"{name}.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        for member, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(f"{name}/{member}")
            info.size = len(data)
            handle.addfile(info, io.BytesIO(data))
    return archive


def add_to_run(output: Path, name: str) -> None:
    """List ``name`` among the queries of the current run in ``run_status.json``."""
    path = output / RUN_STATUS_FILENAME
    queries = (
        json.loads(path.read_text(encoding="utf-8"))["queries"] if path.exists() else {}
    )
    queries[name] = {"status": "complete"}
    path.write_text(json.dumps({"queries": queries}), encoding="utf-8")


def write_record(
    output: Path,
    name: str,
    *,
    status: str = "complete",
    error: str | None = None,
    current: bool = True,
) -> None:
    """Append the record of a query; a complete one names the archive on disk.

    A ``current`` query is also listed in ``run_status.json``.
    """
    archive = output / f"{name}.tar.gz"
    info = None
    if status == "complete":
        info = {
            "filename": archive.name,
            "bytes": archive.stat().st_size,
            "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        }
    record = QueryRecord(
        query=name,
        status=status,
        error=error,
        input_name=f"{name}.faa",
        input_sha256=None,
        settings={},
        database_version="v1",
        software_version=SOFTWARE_VERSION,
        finished_at="2026-01-01T00:00:00+00:00",
        archive=info,
        summary=None,
    )
    with (output / RESULTS_FILENAME).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(asdict(record)) + "\n")
    if current:
        add_to_run(output, name)


def finished_run(output: Path) -> None:
    """Leave a finished run of query ``q1`` with good summaries."""
    members = {"q1.summary.tab": GOOD_TABLE, "q1.final_summary.tsv": GOOD_TABLE}
    write_archive(output, "q1", members)
    write_record(output, "q1")
    (output / "gvclass_summary.tsv").write_text(GOOD_TABLE, encoding="utf-8")


@pytest.mark.parametrize(
    "shared_directory",
    [None, "contig_screen", "species_tree", ".species_tree_work", ".gvclass_work"],
)
def test_output_validator_accepts_a_finished_run(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, shared_directory: str | None
) -> None:
    """Shared output folders are not queries with missing files."""
    finished_run(tmp_path)
    if shared_directory:
        (tmp_path / shared_directory / "q1").mkdir(parents=True)

    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        validation = validate_pipeline_outputs(tmp_path, verbose=True)

    assert validation == {"success": True, "issues": [], "warnings": []}
    assert "Validation successful" in caplog.messages
    assert all(symbol not in caplog.text for symbol in ("✅", "❌", "⚠", "✓", "\x1b"))


def test_a_missing_summary_member_is_an_issue(tmp_path: Path) -> None:
    write_archive(tmp_path, "q1", {"q1.summary.tab": GOOD_TABLE})
    write_record(tmp_path, "q1")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["success"] is False
    assert validation["issues"] == ["Archive member q1/q1.final_summary.tsv is missing"]


@pytest.mark.parametrize("text", ["", "query\ttaxonomy_majority\n"])
def test_an_empty_summary_member_is_an_issue(tmp_path: Path, text: str) -> None:
    write_archive(
        tmp_path, "q1", {"q1.summary.tab": text, "q1.final_summary.tsv": GOOD_TABLE}
    )
    write_record(tmp_path, "q1")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["success"] is False
    assert validation["issues"] == [
        "Archive member q1/q1.summary.tab is empty or has no data row"
    ]


@pytest.mark.parametrize(
    "text", ["query\ttaxonomy_majority\n\n", "query\ttaxonomy_majority\n\t\n"]
)
def test_a_summary_member_with_only_a_blank_row_is_an_issue(
    tmp_path: Path, text: str
) -> None:
    write_archive(
        tmp_path, "q1", {"q1.summary.tab": text, "q1.final_summary.tsv": text}
    )
    write_record(tmp_path, "q1")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["success"] is False
    assert validation["issues"] == [
        "Archive member q1/q1.summary.tab is empty or has no data row",
        "Archive member q1/q1.final_summary.tsv is empty or has no data row",
    ]


def test_a_summary_member_with_the_row_of_another_query_is_an_issue(
    tmp_path: Path,
) -> None:
    write_archive(
        tmp_path,
        "q1",
        {
            "q1.summary.tab": "query\ttaxonomy_majority\nother\ttax\n",
            "q1.final_summary.tsv": GOOD_TABLE,
        },
    )
    write_record(tmp_path, "q1")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["issues"] == [
        "Archive member q1/q1.summary.tab has no row for query q1"
    ]


def test_a_summary_member_finds_its_row_after_a_blank_line(tmp_path: Path) -> None:
    table = "query\ttaxonomy_majority\n\nq1\ttax\n"
    write_archive(
        tmp_path, "q1", {"q1.summary.tab": table, "q1.final_summary.tsv": table}
    )
    write_record(tmp_path, "q1")

    assert validate_pipeline_outputs(tmp_path)["success"] is True


def test_a_summary_member_without_the_key_columns_is_an_issue(tmp_path: Path) -> None:
    write_archive(
        tmp_path,
        "q1",
        {"q1.summary.tab": "wrong\nq1\n", "q1.final_summary.tsv": GOOD_TABLE},
    )
    write_record(tmp_path, "q1")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["issues"] == [
        (
            "Archive member q1/q1.summary.tab lacks required columns: "
            "query, taxonomy_majority"
        )
    ]


def test_a_failed_last_record_is_an_issue(tmp_path: Path) -> None:
    finished_run(tmp_path)
    write_record(tmp_path, "q2", status="failed", error="HMM search failed")
    write_record(tmp_path, "q3", status="invalid_input", error="too short")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["issues"] == [
        "Query q2 is failed: HMM search failed",
        "Query q3 is invalid_input: too short",
    ]


def test_only_the_last_record_of_a_query_counts(tmp_path: Path) -> None:
    write_record(tmp_path, "q1", status="failed", error="first attempt failed")
    finished_run(tmp_path)

    assert validate_pipeline_outputs(tmp_path)["success"] is True


def test_a_record_of_an_input_that_the_run_no_longer_has_is_not_validated(
    tmp_path: Path,
) -> None:
    finished_run(tmp_path)
    write_record(
        tmp_path, "removed", status="invalid_input", error="too short", current=False
    )

    assert validate_pipeline_outputs(tmp_path) == {
        "success": True,
        "issues": [],
        "warnings": [],
    }


def test_a_query_of_the_run_without_a_record_is_an_issue(tmp_path: Path) -> None:
    finished_run(tmp_path)
    add_to_run(tmp_path, "q2")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["issues"] == ["Query q2 has no record in results.jsonl"]


def test_a_missing_or_resized_archive_is_an_issue(tmp_path: Path) -> None:
    finished_run(tmp_path)
    write_archive(tmp_path, "q2", {"q2.summary.tab": GOOD_TABLE})
    write_record(tmp_path, "q2")
    (tmp_path / "q1.tar.gz").unlink()
    (tmp_path / "q2.tar.gz").write_bytes(b"short")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["issues"] == [
        "Archive of query q1 is missing or has changed in size",
        "Archive of query q2 is missing or has changed in size",
    ]


def test_an_output_without_a_run_status_is_an_issue(tmp_path: Path) -> None:
    validation = validate_pipeline_outputs(tmp_path)

    assert validation["success"] is False
    assert validation["issues"] == ["run_status.json is missing or unreadable"]


def test_a_missing_output_directory_is_an_issue(tmp_path: Path) -> None:
    validation = validate_pipeline_outputs(tmp_path / "absent")

    assert validation["success"] is False
    assert "does not exist" in validation["issues"][0]


@pytest.mark.parametrize(
    ("table", "expected"),
    [
        ("", "Combined summary file is empty"),
        (
            "wrong\nq1\n",
            "Combined summary missing required columns: query, taxonomy_majority",
        ),
    ],
)
def test_a_bad_combined_summary_is_an_issue(
    tmp_path: Path, table: str, expected: str
) -> None:
    finished_run(tmp_path)
    (tmp_path / "gvclass_summary.tsv").write_text(table, encoding="utf-8")

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["issues"] == [expected]


def test_a_combined_summary_without_rows_is_a_warning(tmp_path: Path) -> None:
    finished_run(tmp_path)
    (tmp_path / "gvclass_summary.tsv").write_text(
        "query\ttaxonomy_majority\n", encoding="utf-8"
    )

    validation = validate_pipeline_outputs(tmp_path)

    assert validation["success"] is True
    assert validation["warnings"] == ["Combined summary has header but no data rows"]


def test_diagnostics_are_plain(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    write_record(tmp_path, "q1", status="failed", error="boom")
    (tmp_path / "gvclass_summary.tsv").write_text(
        "query\ttaxonomy_majority\n", encoding="utf-8"
    )

    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        validate_pipeline_outputs(tmp_path, verbose=True)

    assert "Validation found 1 issue(s)" in caplog.messages
    assert "1 warning(s) found" in caplog.messages
    assert all(symbol not in caplog.text for symbol in ("✅", "❌", "⚠", "✓", "\x1b"))
