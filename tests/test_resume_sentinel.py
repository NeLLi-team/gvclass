"""Regression tests for the v1.4.3 resume SUCCESS sentinel + atomic tar path.

Covers Phase 1.3 of the remediation plan:
* SUCCESS sentinel is written only after summary + tar are durable on disk.
* Resume accepts the sentinel without re-inspecting the archive.
* Resume falls back to the legacy summary+tar pair only when the sentinel
  is missing AND the tar is well-formed (``tarfile.is_tarfile``).
* A corrupted or truncated tar causes the legacy fallback to reject the
  query so it re-runs.
* ``_archive_query_output`` never leaves a truncated ``.tar.gz`` on failure.
"""

from __future__ import annotations

import json
import logging
import tarfile
from pathlib import Path
from unittest.mock import Mock

import pytest


LOGGER = logging.getLogger("test_resume_sentinel")


# ---------------------------------------------------------------------------
# SUCCESS sentinel writer.
# ---------------------------------------------------------------------------


def _seed_post_process_artefacts(output_base: Path, query_name: str) -> tuple[Path, Path]:
    output_base.mkdir(parents=True, exist_ok=True)
    summary_file = output_base / f"{query_name}.summary.tab"
    tar_file = output_base / f"{query_name}.tar.gz"
    summary_file.write_text("query\ttaxonomy_majority\nq\td_NCLDV\n")
    # Not a real tar, but exists for sentinel payload checksum.
    tar_file.write_bytes(b"fake-tar-bytes")
    return summary_file, tar_file


def test_write_success_sentinel_records_version_and_checksums(tmp_path: Path) -> None:
    from src.pipeline.query_processing_engine import (
        SUCCESS_SENTINEL_VERSION,
        _write_success_sentinel,
    )

    output_base = tmp_path / "out"
    summary_file, tar_file = _seed_post_process_artefacts(output_base, "q1")

    sentinel_path = _write_success_sentinel("q1", output_base, LOGGER)

    assert sentinel_path is not None
    assert sentinel_path.name == "q1.SUCCESS"
    payload = json.loads(sentinel_path.read_text())
    assert payload["sentinel_version"] == SUCCESS_SENTINEL_VERSION
    assert payload["query"] == "q1"
    assert payload["software_version"].startswith("v")
    assert payload["tar_filename"] == tar_file.name
    assert payload["summary_filename"] == summary_file.name
    # SHA-256 hex digest is 64 hex chars.
    assert len(payload["tar_sha256"]) == 64
    assert len(payload["summary_sha256"]) == 64
    assert payload["tar_bytes"] == tar_file.stat().st_size


def test_write_success_sentinel_skips_when_artefact_missing(tmp_path: Path) -> None:
    from src.pipeline.query_processing_engine import _write_success_sentinel

    output_base = tmp_path / "out"
    output_base.mkdir()
    (output_base / "q1.summary.tab").write_text("query\n")
    # No tar file on disk -> sentinel must not be written.

    logger = Mock()
    assert _write_success_sentinel("q1", output_base, logger) is None
    assert not (output_base / "q1.SUCCESS").exists()
    logger.warning.assert_called()


def test_clear_success_sentinel_removes_prior_marker(tmp_path: Path) -> None:
    from src.pipeline.query_processing_engine import _clear_success_sentinel

    output_base = tmp_path / "out"
    output_base.mkdir()
    sentinel = output_base / "q1.SUCCESS"
    sentinel.write_text("{}")
    _clear_success_sentinel("q1", output_base, LOGGER)
    assert not sentinel.exists()


# ---------------------------------------------------------------------------
# Atomic tar writer.
# ---------------------------------------------------------------------------


def test_archive_query_output_is_atomic_on_failure(tmp_path: Path, monkeypatch) -> None:
    """Simulate an exception mid-tar and assert no ``.tar.gz`` is published."""
    from src.pipeline.query_processing_engine import _archive_query_output

    output_base = tmp_path / "out"
    output_base.mkdir()
    query_output_dir = tmp_path / "work" / "q1"
    query_output_dir.mkdir(parents=True)
    (query_output_dir / "hello.txt").write_text("hello")

    original_open = tarfile.open

    def exploding_tarfile_open(*args, **kwargs):
        raise RuntimeError("simulated tar failure")

    monkeypatch.setattr(tarfile, "open", exploding_tarfile_open)
    with pytest.raises(RuntimeError, match="simulated tar failure"):
        _archive_query_output("q1", query_output_dir, output_base, LOGGER)
    monkeypatch.setattr(tarfile, "open", original_open)

    assert not (output_base / "q1.tar.gz").exists()
    assert not (output_base / "q1.tar.gz.part").exists()


def test_archive_query_output_publishes_valid_archive(tmp_path: Path) -> None:
    from src.pipeline.query_processing_engine import _archive_query_output

    output_base = tmp_path / "out"
    output_base.mkdir()
    query_output_dir = tmp_path / "work" / "q1"
    query_output_dir.mkdir(parents=True)
    (query_output_dir / "hello.txt").write_text("hello")

    _archive_query_output("q1", query_output_dir, output_base, LOGGER)

    final_tar = output_base / "q1.tar.gz"
    assert final_tar.exists()
    assert not (output_base / "q1.tar.gz.part").exists()
    assert tarfile.is_tarfile(final_tar)


# ---------------------------------------------------------------------------
# Resume filter semantics.
# ---------------------------------------------------------------------------


def _make_tar(path: Path, content: bytes = b"ok") -> None:
    with tarfile.open(path, "w:gz") as tar:
        payload = path.parent / "payload.txt"
        payload.write_bytes(content)
        tar.add(payload, arcname="payload.txt")
        payload.unlink()


def test_query_is_resume_complete_accepts_sentinel(tmp_path: Path) -> None:
    from src.pipeline.parallel_runner import _query_is_resume_complete

    output_path = tmp_path / "out"
    output_path.mkdir()
    (output_path / "q1.SUCCESS").write_text("{}")

    logger = Mock()
    assert _query_is_resume_complete(output_path, "q1", logger) is True
    logger.warning.assert_not_called()


def test_query_is_resume_complete_legacy_fallback_accepts_valid_pair(tmp_path: Path) -> None:
    from src.pipeline.parallel_runner import _query_is_resume_complete

    output_path = tmp_path / "out"
    output_path.mkdir()
    (output_path / "q1.summary.tab").write_text("query\n")
    _make_tar(output_path / "q1.tar.gz")

    logger = Mock()
    assert _query_is_resume_complete(output_path, "q1", logger) is True
    # We log INFO (not warning) when falling back.
    logger.info.assert_called()


def test_query_is_resume_complete_rejects_corrupt_tar(tmp_path: Path) -> None:
    from src.pipeline.parallel_runner import _query_is_resume_complete

    output_path = tmp_path / "out"
    output_path.mkdir()
    (output_path / "q1.summary.tab").write_text("query\n")
    # Not a tar file at all.
    (output_path / "q1.tar.gz").write_bytes(b"not a real tar.gz")

    logger = Mock()
    assert _query_is_resume_complete(output_path, "q1", logger) is False
    logger.warning.assert_called()


def test_query_is_resume_complete_rejects_missing_pair(tmp_path: Path) -> None:
    from src.pipeline.parallel_runner import _query_is_resume_complete

    output_path = tmp_path / "out"
    output_path.mkdir()
    # Summary only, no tar.
    (output_path / "q1.summary.tab").write_text("query\n")

    logger = Mock()
    assert _query_is_resume_complete(output_path, "q1", logger) is False


def test_count_resume_skips_counts_sentinel_and_legacy(tmp_path: Path) -> None:
    from src.bin.gvclass_cli import count_resume_skips

    output_path = tmp_path / "out"
    output_path.mkdir()

    # Query A has a sentinel (modern).
    (output_path / "qa.SUCCESS").write_text("{}")
    (output_path / "qa.summary.tab").write_text("")
    (output_path / "qa.tar.gz").write_bytes(b"")

    # Query B has the legacy pair only.
    (output_path / "qb.summary.tab").write_text("")
    (output_path / "qb.tar.gz").write_bytes(b"")

    # Query C has a partial output (summary only) and should NOT count.
    (output_path / "qc.summary.tab").write_text("")

    assert count_resume_skips(output_path) == 2


# ---------------------------------------------------------------------------
# Atomic DB install.
# ---------------------------------------------------------------------------


def test_atomic_swap_directory_rolls_back_on_failure(tmp_path: Path, monkeypatch) -> None:
    from src.utils.database_manager import DatabaseManager

    db_path = tmp_path / "resources"
    staging_path = tmp_path / "resources.new"
    db_path.mkdir()
    (db_path / "existing.txt").write_text("old")
    staging_path.mkdir()
    (staging_path / "new.txt").write_text("new")

    # Force the second os.replace (staging -> db_path) to fail.
    import os as os_module

    real_replace = os_module.replace
    call_state = {"count": 0}

    def flaky_replace(src, dst):
        call_state["count"] += 1
        if call_state["count"] == 2:
            raise OSError("simulated rename failure")
        real_replace(src, dst)

    monkeypatch.setattr(os_module, "replace", flaky_replace)
    with pytest.raises(OSError, match="simulated rename failure"):
        DatabaseManager._atomic_swap_directory(staging_path, db_path)

    # Original db_path must have been restored.
    assert db_path.exists()
    assert (db_path / "existing.txt").exists()
    assert (db_path / "existing.txt").read_text() == "old"


def test_atomic_swap_directory_replaces_previous_and_cleans_backup(tmp_path: Path) -> None:
    from src.utils.database_manager import DatabaseManager

    db_path = tmp_path / "resources"
    staging_path = tmp_path / "resources.new"
    db_path.mkdir()
    (db_path / "existing.txt").write_text("old")
    staging_path.mkdir()
    (staging_path / "new.txt").write_text("new")

    DatabaseManager._atomic_swap_directory(staging_path, db_path)

    assert db_path.exists()
    assert (db_path / "new.txt").read_text() == "new"
    assert not (db_path / "existing.txt").exists()
    assert not (tmp_path / "resources.old").exists()
    assert not staging_path.exists()


def test_clear_prior_outputs_removes_legacy_artefacts(tmp_path: Path) -> None:
    """Regression for Codex-audit finding: clearing only SUCCESS is not enough.

    A rerun that crashes after the sentinel is cleared but before the new
    tar is published would otherwise leave the previous run's
    ``.summary.tab`` + ``.tar.gz`` pair on disk. ``_query_is_resume_complete``
    would then still return True via the legacy fallback and the rerun
    would be silently skipped next time. ``_clear_prior_outputs`` must
    wipe every resumable marker up front so the post-crash state is
    obviously incomplete.
    """
    from src.pipeline.parallel_runner import _query_is_resume_complete
    from src.pipeline.query_processing_engine import _clear_prior_outputs

    output_base = tmp_path / "out"
    output_base.mkdir()
    (output_base / "q1.SUCCESS").write_text("{}")
    (output_base / "q1.summary.tab").write_text("query\nq1\n")
    _make_tar(output_base / "q1.tar.gz", content=b"stale-payload")
    # Also drop a dangling .part from a previous crash to make sure it is
    # wiped along with the final artefacts.
    (output_base / "q1.tar.gz.part").write_bytes(b"half-written")
    (output_base / "q1.contamination_candidates.tsv").write_text("stale\n")

    _clear_prior_outputs("q1", output_base, LOGGER)

    assert not (output_base / "q1.SUCCESS").exists()
    assert not (output_base / "q1.summary.tab").exists()
    assert not (output_base / "q1.tar.gz").exists()
    assert not (output_base / "q1.tar.gz.part").exists()
    assert not (output_base / "q1.contamination_candidates.tsv").exists()

    # After cleanup, resume must NOT treat the query as complete via
    # either the sentinel or the legacy fallback.
    assert _query_is_resume_complete(output_base, "q1", Mock()) is False


def test_resume_all_skipped_does_not_crash_on_zero_workers(tmp_path: Path) -> None:
    """Regression: if every query has a SUCCESS sentinel, the flow must not
    instantiate ThreadPoolExecutor(max_workers=0)."""
    from src.pipeline.parallel_runner import gvclass_flow

    query_dir = tmp_path / "queries"
    output_dir = tmp_path / "out"
    query_dir.mkdir()
    output_dir.mkdir()
    (query_dir / "q1.fna").write_text(">c1\nACGT\n")
    (query_dir / "q2.fna").write_text(">c2\nACGT\n")

    for name in ("q1", "q2"):
        (output_dir / f"{name}.SUCCESS").write_text("{}")
        (output_dir / f"{name}.summary.tab").write_text("query\n" + name + "\n")
        _make_tar(output_dir / f"{name}.tar.gz")

    # Run with resume=True; the flow should detect zero remaining queries
    # and short-circuit without spawning any workers. It still writes a
    # final summary from whatever is on disk (empty in this synthetic case).
    # The synthetic FNAs are far below the Phase 3.2 20 kb floor, so pass
    # allow_short=True to let validate_query_directory accept them.
    summary_file = gvclass_flow(
        query_dir=str(query_dir),
        output_dir=str(output_dir),
        database_path=None,
        total_threads=4,
        resume=True,
        allow_short=True,
    )
    assert Path(summary_file).exists()
