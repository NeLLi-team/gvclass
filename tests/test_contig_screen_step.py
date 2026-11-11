"""The contig-screen step of a run: one subprocess per bin, totals in the summary row.

The subprocess is the stand-in of ``tests/fake_screen.py``: it follows the command
line of ``python -m src.bin.gvclass_screen`` and writes ``contig_screen_summary.tsv``
the way it does. The mode of a bin comes from the first word of its name.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pytest

from src.pipeline import contig_screen_step as step
from tests.fake_screen import INPUT_DIGEST, install_fake_screen, seed_summary

MANIFEST_DIGEST = hashlib.sha256(b"{}").hexdigest()
LOGGER = logging.getLogger("test_contig_screen_step")
NAN = math.nan
#: What the marker model leaves in the summary row of every query before the screen.
MARKER_MODEL = {
    "estimated_contamination": 41.5,
    "estimated_contamination_strategy": "extra_trees_v1",
    "contamination_query_support": "limited_support",
    "contamination_model_reliability": "unvalidated_legacy",
    "contamination_type": "phage",
}


@pytest.fixture(autouse=True)
def fake_screen(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    install_fake_screen(tmp_path, monkeypatch)


@pytest.fixture
def database(tmp_path: Path) -> Path:
    """A database with the add-on manifest in place."""
    root = tmp_path / "db"
    (root / "contig_screen").mkdir(parents=True)
    (root / "contig_screen" / "manifest.json").write_text("{}", encoding="utf-8")
    return root


def options(**overrides: Any) -> step.ScreenOptions:
    values = {
        "enabled": True,
        "contigs_mode": False,
        "resume": False,
        "plan_workers": lambda n_bins: (2, 3),
    }
    return step.ScreenOptions(**{**values, **overrides})


def make_queries(
    tmp_path: Path, names: list[str]
) -> tuple[list[dict[str, Any]], dict[str, Path]]:
    """Return completed results and the query files; ``.faa`` names are proteins."""
    folder = tmp_path / "queries"
    folder.mkdir(exist_ok=True)
    results = []
    files = {}
    for name in names:
        suffix = ".faa" if name.startswith("prot") else ".fna"
        path = folder / f"{name}{suffix}"
        path.write_text(">c1\nACGT\n", encoding="utf-8")
        files[name] = path
        results.append(
            {
                "query": name,
                "status": "complete",
                "summary_data": {"query": name, **MARKER_MODEL},
            }
        )
    return results, files


def run_step(
    tmp_path: Path,
    database: Path,
    names: list[str],
    **overrides: Any,
) -> list[dict[str, Any]]:
    results, files = make_queries(tmp_path, names)
    step.add_screen_columns(
        results, files, tmp_path / "out", database, options(**overrides), LOGGER
    )
    return results


def screen_columns(result: dict[str, Any]) -> dict[str, Any]:
    data = result["summary_data"]
    return {key: data[key] for key in data if key.startswith("screen_")}


def contamination_columns(result: dict[str, Any]) -> dict[str, Any]:
    data = result["summary_data"]
    return {key: data[key] for key in MARKER_MODEL}


def test_a_completed_screen_adds_its_numbers_to_the_summary_row(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["ok_a"])

    assert screen_columns(results[0]) == {
        "screen_status": "ok",
        "screen_cellular_pct": 10.0,
        "screen_giant_virus_pct": 20.0,
        "screen_flagged_pct": 25.0,
    }


def test_the_bin_runs_in_its_own_folder_outside_the_query_folder(
    tmp_path: Path, database: Path
) -> None:
    results, files = make_queries(tmp_path, ["ok_a"])

    step.add_screen_columns(
        results, files, tmp_path / "out", database, options(), LOGGER
    )

    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    args = json.loads((folder / "args.json").read_text(encoding="utf-8"))
    assert args == [
        str(files["ok_a"]),
        "-o",
        str(folder),
        "-d",
        str(database),
        "-t",
        "3",
    ]
    assert not (tmp_path / "out" / "ok_a").exists()


def test_a_screened_bin_reports_its_stage_to_the_progress_display(
    tmp_path: Path,
    database: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setenv("GVCLASS_PROGRESS_EVENTS", "1")

    run_step(tmp_path, database, ["ok_a", "prot_b"])

    events = [
        json.loads(line.removeprefix("GVCLASS_PROGRESS "))
        for line in capsys.readouterr().err.splitlines()
        if line.startswith("GVCLASS_PROGRESS ")
    ]
    assert [(e["query"], e["stage"], e["progress"]) for e in events] == [
        ("ok_a", "contig_screen", 100)
    ]


def test_a_disabled_screen_runs_nothing(tmp_path: Path, database: Path) -> None:
    results = run_step(tmp_path, database, ["ok_a"], enabled=False)

    assert results[0]["summary_data"]["screen_status"] == "disabled"
    assert not (tmp_path / "out").exists()


def test_a_disabled_screen_does_not_report_a_missing_addon(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger=LOGGER.name):
        results = run_step(tmp_path, tmp_path / "db", ["ok_a"], enabled=False)

    assert results[0]["summary_data"]["screen_status"] == "disabled"
    assert caplog.records == []


def test_protein_input_is_not_applicable(tmp_path: Path, database: Path) -> None:
    results = run_step(tmp_path, database, ["prot_a", "ok_b"])

    assert results[0]["summary_data"]["screen_status"] == "not_applicable"
    assert results[1]["summary_data"]["screen_status"] == "ok"
    assert not (tmp_path / "out" / "contig_screen" / "prot_a").exists()


def test_a_protein_fasta_is_not_applicable(tmp_path: Path, database: Path) -> None:
    results, files = make_queries(tmp_path, ["ok_a"])
    protein = tmp_path / "queries" / "ok_p.fasta"
    protein.write_text(
        ">p\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD\n", encoding="utf-8"
    )
    results.append({"query": "ok_p", "status": "complete", "summary_data": {}})
    files["ok_p"] = protein

    step.add_screen_columns(
        results, files, tmp_path / "out", database, options(), LOGGER
    )

    assert results[1]["summary_data"]["screen_status"] == "not_applicable"
    assert results[0]["summary_data"]["screen_status"] == "ok"


def test_contigs_mode_is_not_applicable(tmp_path: Path, database: Path) -> None:
    results = run_step(tmp_path, database, ["ok_a"], contigs_mode=True)

    assert results[0]["summary_data"]["screen_status"] == "not_applicable"
    assert not (tmp_path / "out").exists()


def test_a_missing_addon_is_reported_once_with_the_install_command(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger=LOGGER.name):
        results = run_step(tmp_path, tmp_path / "db", ["ok_a", "ok_b", "ok_c"])

    assert [r["summary_data"]["screen_status"] for r in results] == [
        "not_installed"
    ] * 3
    lines = [r.getMessage() for r in caplog.records if "screen --install" in r.message]
    assert len(lines) == 1
    assert not (tmp_path / "out").exists()


def test_not_installed_leaves_the_percentages_missing(tmp_path: Path) -> None:
    results = run_step(tmp_path, tmp_path / "db", ["ok_a"])

    columns = screen_columns(results[0])
    assert math.isnan(columns["screen_cellular_pct"])
    assert math.isnan(columns["screen_giant_virus_pct"])
    assert math.isnan(columns["screen_flagged_pct"])


def test_a_failed_query_has_no_row_to_screen(tmp_path: Path, database: Path) -> None:
    results, files = make_queries(tmp_path, ["ok_a"])
    results.append({"query": "ok_bad", "status": "failed", "error": "x"})
    files["ok_bad"] = files["ok_a"]

    step.add_screen_columns(
        results, files, tmp_path / "out", database, options(), LOGGER
    )

    assert "summary_data" not in results[1]
    assert not (tmp_path / "out" / "contig_screen" / "ok_bad").exists()


def test_a_crashed_screen_fails_that_bin_and_logs_a_warning(
    tmp_path: Path, database: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING, logger=LOGGER.name):
        results = run_step(tmp_path, database, ["crash_a", "ok_b"])

    assert results[0]["summary_data"]["screen_status"] == "failed"
    assert math.isnan(results[0]["summary_data"]["screen_flagged_pct"])
    assert results[1]["summary_data"]["screen_status"] == "ok"
    warning = " ".join(r.getMessage() for r in caplog.records)
    assert "crash_a" in warning
    assert "out of memory" in warning


def test_a_screen_that_exits_with_an_error_but_wrote_its_summary_counts(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["half_a"])

    assert screen_columns(results[0])["screen_status"] == "partial"
    assert screen_columns(results[0])["screen_cellular_pct"] == 10.0
    assert math.isnan(screen_columns(results[0])["screen_flagged_pct"])


def test_two_failed_screens_make_the_bin_failed(tmp_path: Path, database: Path) -> None:
    results = run_step(tmp_path, database, ["dead_a"])

    assert results[0]["summary_data"]["screen_status"] == "failed"


def test_a_bin_that_one_screen_cannot_score_is_partial(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["single_a", "short_b"])

    first = screen_columns(results[0])
    second = screen_columns(results[1])
    assert (first["screen_status"], first["screen_cellular_pct"]) == ("partial", 10.0)
    assert math.isnan(first["screen_giant_virus_pct"])
    assert (second["screen_status"], second["screen_giant_virus_pct"]) == (
        "partial",
        20.0,
    )
    assert math.isnan(second["screen_cellular_pct"])


@pytest.mark.parametrize(
    ("cellular", "giant_virus", "expected"),
    [
        ("ok", "ok", "ok"),
        ("ok", "not_scored", "partial"),
        ("not_scored", "ok", "partial"),
        ("not_scored", "not_scored", "partial"),
        ("ok", "failed", "partial"),
        ("failed", "ok", "partial"),
        ("failed", "failed", "failed"),
        ("failed", "not_scored", "failed"),
        ("not_scored", "failed", "failed"),
    ],
)
def test_the_status_of_a_bin_follows_the_status_of_its_two_screens(
    cellular: str, giant_virus: str, expected: str
) -> None:
    assert step.screen_status(cellular, giant_virus) == expected


def test_resume_reuses_a_finished_summary_without_running_the_screen(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    seed_summary(
        folder, "ok_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n", MANIFEST_DIGEST
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert screen_columns(results[0]) == {
        "screen_status": "ok",
        "screen_cellular_pct": 7.0,
        "screen_giant_virus_pct": 8.0,
        "screen_flagged_pct": 9.0,
    }
    assert not (folder / "args.json").exists()


def test_resume_runs_the_screen_for_a_bin_without_a_finished_summary(
    tmp_path: Path, database: Path
) -> None:
    seed_summary(
        tmp_path / "out" / "contig_screen" / "ok_a",
        "ok_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, ["ok_a", "ok_b"], resume=True)

    assert results[1]["summary_data"]["screen_flagged_pct"] == 25.0
    assert (tmp_path / "out" / "contig_screen" / "ok_b" / "args.json").exists()


def test_resume_runs_the_screen_again_when_the_earlier_screen_failed(
    tmp_path: Path, database: Path
) -> None:
    seed_summary(
        tmp_path / "out" / "contig_screen" / "ok_a",
        "ok_a\tfailed\tfailed\t1000\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_status"] == "ok"
    assert (tmp_path / "out" / "contig_screen" / "ok_a" / "args.json").exists()


def test_resume_runs_the_screen_again_when_one_of_two_screens_failed(
    tmp_path: Path, database: Path
) -> None:
    seed_summary(
        tmp_path / "out" / "contig_screen" / "ok_a",
        "ok_a\tok\tfailed\t1000\t100\t10.00\tNaN\tNaN\tNaN\tNaN\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_status"] == "ok"


def test_resume_reuses_a_summary_in_which_a_screen_could_not_score_the_bin(
    tmp_path: Path, database: Path
) -> None:
    seed_summary(
        tmp_path / "out" / "contig_screen" / "ok_a",
        "ok_a\tok\tnot_scored\t1000\t100\t10.00\tNaN\tNaN\tNaN\tNaN\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_status"] == "partial"
    assert not (tmp_path / "out" / "contig_screen" / "ok_a" / "args.json").exists()


def test_resume_runs_the_screen_again_for_a_summary_without_the_totals(
    tmp_path: Path, database: Path
) -> None:
    seed_summary(
        tmp_path / "out" / "contig_screen" / "ok_a",
        "ok_a\tok\tok\n",
        None,
        header="bin\tcellular_status\tgiant_virus_status\n",
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_flagged_pct"] == 25.0


def test_resume_runs_the_screen_again_when_the_addon_manifest_changed(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    seed_summary(folder, "ok_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n", "0" * 64)

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_flagged_pct"] == 25.0
    assert (folder / "args.json").exists()
    assert MANIFEST_DIGEST in (folder / "contig_screen_summary.tsv").read_text()


def test_a_summary_written_by_a_run_is_reused_by_the_next_resume(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    run_step(tmp_path, database, ["ok_a"])
    (folder / "args.json").unlink()

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_status"] == "ok"
    assert not (folder / "args.json").exists()


def test_a_bin_that_was_not_screened_for_a_missing_addon_is_screened_by_a_later_resume(
    tmp_path: Path,
) -> None:
    database = tmp_path / "db"
    database.mkdir()
    first = run_step(tmp_path, database, ["ok_a"])
    (database / "contig_screen").mkdir()
    (database / "contig_screen" / "manifest.json").write_text("{}", encoding="utf-8")

    second = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert first[0]["summary_data"]["screen_status"] == "not_installed"
    assert second[0]["summary_data"]["screen_status"] == "ok"
    assert second[0]["summary_data"]["screen_flagged_pct"] == 25.0


def test_a_disabled_resume_neither_reuses_nor_deletes_saved_results(
    tmp_path: Path, database: Path
) -> None:
    run_step(tmp_path, database, ["ok_a"])
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    (folder / "ok_a.cellular.tsv").write_text("kept\n", encoding="utf-8")
    before = {path.name: path.read_bytes() for path in folder.iterdir()}

    results = run_step(tmp_path, database, ["ok_a"], resume=True, enabled=False)

    assert screen_columns(results[0])["screen_status"] == "disabled"
    assert math.isnan(screen_columns(results[0])["screen_flagged_pct"])
    assert {path.name: path.read_bytes() for path in folder.iterdir()} == before


def test_without_resume_the_screen_runs_again_and_a_stale_summary_is_not_used(
    tmp_path: Path, database: Path
) -> None:
    seed_summary(
        tmp_path / "out" / "contig_screen" / "crash_a",
        "crash_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, ["crash_a"])

    assert results[0]["summary_data"]["screen_status"] == "failed"
    assert not (
        tmp_path / "out" / "contig_screen" / "crash_a" / "contig_screen_summary.tsv"
    ).exists()


def test_the_pool_has_at_most_as_many_threads_as_the_worker_plan(
    tmp_path: Path, database: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[int] = []

    class Spy(ThreadPoolExecutor):
        def __init__(self, max_workers: int) -> None:
            seen.append(max_workers)
            super().__init__(max_workers=max_workers)

    monkeypatch.setattr(step, "ThreadPoolExecutor", Spy)

    run_step(
        tmp_path, database, ["ok_a", "ok_b", "ok_c"], plan_workers=lambda n: (2, 3)
    )
    run_step(tmp_path, database, ["ok_a"], plan_workers=lambda n: (4, 3))

    assert seen == [2, 1]


def test_the_worker_plan_is_asked_for_the_bins_that_need_a_run(
    tmp_path: Path, database: Path
) -> None:
    asked: list[int] = []
    folder = tmp_path / "out" / "contig_screen" / "ok_c"
    seed_summary(
        folder, "ok_c\tok\tok\t1\t1\t1.00\t1\t1.00\t1\t1.00\n", MANIFEST_DIGEST
    )

    run_step(
        tmp_path,
        database,
        ["ok_a", "prot_b", "ok_c", "ok_d"],
        resume=True,
        plan_workers=lambda n: asked.append(n) or (3, 5),
    )

    assert asked == [2]


def test_the_threads_of_the_plan_go_to_every_screen(
    tmp_path: Path, database: Path
) -> None:
    run_step(tmp_path, database, ["ok_a"], plan_workers=lambda n: (1, 7))

    args = json.loads(
        (tmp_path / "out" / "contig_screen" / "ok_a" / "args.json").read_text()
    )
    assert args[args.index("-t") + 1] == "7"


def test_no_plan_is_asked_for_when_no_bin_needs_a_run(
    tmp_path: Path, database: Path
) -> None:
    def refuse(n: int) -> tuple[int, int]:
        raise AssertionError(n)

    results = run_step(tmp_path, database, ["prot_a"], plan_workers=refuse)

    assert results[0]["summary_data"]["screen_status"] == "not_applicable"


def test_an_unexpected_error_fails_that_bin_and_never_the_run(
    tmp_path: Path,
    database: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def explode(*args: object, **kwargs: object) -> None:
        msg = "disk full"
        raise OSError(msg)

    monkeypatch.setattr(step.subprocess, "Popen", explode)

    with caplog.at_level(logging.WARNING, logger=LOGGER.name):
        results = run_step(tmp_path, database, ["ok_a", "prot_b"])

    assert results[0]["summary_data"]["screen_status"] == "failed"
    assert results[1]["summary_data"]["screen_status"] == "not_applicable"
    assert "disk full" in " ".join(r.getMessage() for r in caplog.records)


def test_the_step_does_not_import_torch() -> None:
    code = (
        "import sys; import src.pipeline.contig_screen_step; "
        "print('torch' in sys.modules, 'src.contig_screen.cellular' in sys.modules)"
    )

    shown = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
        cwd=Path(__file__).resolve().parents[1],
    )

    assert shown.stdout.strip() == "False False"


def refuse_manifest(database: Path) -> str:
    raise PermissionError(13, "Permission denied")


def test_an_unreadable_manifest_fails_the_bins_that_need_it_and_not_the_run(
    tmp_path: Path,
    database: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    monkeypatch.setattr(step.layout, "manifest_digest", refuse_manifest)

    with caplog.at_level(logging.WARNING, logger=LOGGER.name):
        results = run_step(tmp_path, database, ["ok_a", "prot_b"])

    assert results[0]["summary_data"]["screen_status"] == "failed"
    assert math.isnan(results[0]["summary_data"]["screen_flagged_pct"])
    assert results[1]["summary_data"]["screen_status"] == "not_applicable"
    assert "Permission denied" in " ".join(r.getMessage() for r in caplog.records)
    assert not (tmp_path / "out").exists()


def test_an_unreadable_manifest_does_not_touch_a_disabled_screen(
    tmp_path: Path,
    database: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    monkeypatch.setattr(step.layout, "manifest_digest", refuse_manifest)

    with caplog.at_level(logging.INFO, logger=LOGGER.name):
        results = run_step(tmp_path, database, ["ok_a"], enabled=False)

    assert results[0]["summary_data"]["screen_status"] == "disabled"
    assert caplog.records == []


def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def wait_until_gone(pid: int, seconds: float = 10.0) -> bool:
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if not process_exists(pid):
            return True
        time.sleep(0.1)
    return False


def test_the_time_limit_of_a_screen_is_four_hours() -> None:
    assert step.SCREEN_TIMEOUT_S == 4 * 60 * 60


def test_a_screen_that_hangs_is_killed_with_its_children_and_fails_that_bin(
    tmp_path: Path,
    database: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    monkeypatch.setattr(step, "SCREEN_TIMEOUT_S", 1.0)
    started = time.monotonic()

    with caplog.at_level(logging.WARNING, logger=LOGGER.name):
        results = run_step(
            tmp_path, database, ["hang_a", "ok_b"], plan_workers=lambda n: (2, 3)
        )

    child = int(
        (tmp_path / "out" / "contig_screen" / "hang_a" / "child.pid").read_text()
    )
    assert time.monotonic() - started < 60
    assert results[0]["summary_data"]["screen_status"] == "failed"
    assert math.isnan(results[0]["summary_data"]["screen_flagged_pct"])
    assert results[1]["summary_data"]["screen_status"] == "ok"
    assert wait_until_gone(child)
    assert "timed out" in " ".join(r.getMessage() for r in caplog.records)


def test_the_output_of_a_screen_goes_to_a_log_file_and_only_its_tail_is_logged(
    tmp_path: Path, database: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING, logger=LOGGER.name):
        run_step(tmp_path, database, ["noisy_a"])

    log = tmp_path / "out" / "contig_screen" / "noisy_a" / "screen.log"
    warning = " ".join(r.getMessage() for r in caplog.records)
    assert log.stat().st_size > 5_000_000
    assert "last line of noise" in warning
    assert len(warning) < 2_000


def test_resume_runs_the_screen_again_when_the_bin_input_changed(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    seed_summary(
        folder,
        "ok_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n",
        MANIFEST_DIGEST,
        input_digest=hashlib.sha256(b"an earlier input").hexdigest(),
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_flagged_pct"] == 25.0
    assert (folder / "args.json").exists()
    assert INPUT_DIGEST in (folder / "contig_screen_summary.tsv").read_text()


def test_resume_runs_the_screen_again_for_a_summary_that_names_another_bin(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    seed_summary(
        folder,
        "other\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_flagged_pct"] == 25.0
    assert (folder / "args.json").exists()


def test_resume_runs_the_screen_again_for_a_summary_without_the_input_digest(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    seed_summary(
        folder,
        "ok_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\t" + MANIFEST_DIGEST + "\n",
        None,
        header=(
            "bin\tcellular_status\tgiant_virus_status\tinput_bp\tcellular_bp\t"
            "cellular_pct\tgiant_virus_bp\tgiant_virus_pct\tflagged_bp\t"
            "flagged_pct\tmanifest_sha256\n"
        ),
    )

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert results[0]["summary_data"]["screen_flagged_pct"] == 25.0
    assert (folder / "args.json").exists()


@pytest.mark.parametrize("name", ["NA", "007"])
def test_resume_reads_the_bin_name_of_a_summary_as_text(
    tmp_path: Path, database: Path, name: str
) -> None:
    folder = tmp_path / "out" / "contig_screen" / name
    seed_summary(
        folder,
        f"{name}\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n",
        MANIFEST_DIGEST,
    )

    results = run_step(tmp_path, database, [name], resume=True)

    assert screen_columns(results[0])["screen_flagged_pct"] == 9.0
    assert not (folder / "args.json").exists()


def test_an_interrupt_stops_the_running_screens_and_their_children(
    tmp_path: Path, database: Path
) -> None:
    def interrupt(signum: int, frame: object) -> None:
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGALRM, interrupt)
    signal.alarm(2)
    started = time.monotonic()
    try:
        with pytest.raises(KeyboardInterrupt):
            run_step(tmp_path, database, ["hang_a"])
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)

    child = int(
        (tmp_path / "out" / "contig_screen" / "hang_a" / "child.pid").read_text()
    )
    assert time.monotonic() - started < 60
    assert wait_until_gone(child)


def child_of(tmp_path: Path, name: str) -> int:
    return int((tmp_path / "out" / "contig_screen" / name / "child.pid").read_text())


def test_a_child_that_outlives_a_crashed_screen_is_killed(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["orphan_a"])

    assert results[0]["summary_data"]["screen_status"] == "failed"
    assert wait_until_gone(child_of(tmp_path, "orphan_a"))


def spawn_in_own_session() -> subprocess.Popen[bytes]:
    return subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(600)"], start_new_session=True
    )


def test_a_screen_registered_after_the_stop_is_killed_at_once() -> None:
    running = step._Running()
    running.stop()
    process = spawn_in_own_session()

    running.add(process)

    assert process.wait(timeout=10) != 0


def test_a_screen_registered_before_the_stop_is_killed_by_it() -> None:
    running = step._Running()
    process = spawn_in_own_session()
    running.add(process)

    running.stop()

    assert process.wait(timeout=10) != 0


def interrupt_after(seconds: int) -> object:
    """Raise KeyboardInterrupt in the main thread after ``seconds``."""

    def interrupt(number: int, frame: object) -> None:
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGALRM, interrupt)
    signal.alarm(seconds)
    return previous


def test_an_interrupt_cancels_the_bins_that_have_not_started(
    tmp_path: Path, database: Path
) -> None:
    previous = interrupt_after(2)
    try:
        with pytest.raises(KeyboardInterrupt):
            run_step(
                tmp_path,
                database,
                ["hang_a", "hang_b", "hang_c"],
                plan_workers=lambda n: (1, 3),
            )
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)

    folders = tmp_path / "out" / "contig_screen"
    assert (folders / "hang_a" / "args.json").exists()
    assert not (folders / "hang_b" / "args.json").exists()
    assert not (folders / "hang_c" / "args.json").exists()
    assert wait_until_gone(child_of(tmp_path, "hang_a"))


def test_sigterm_stops_the_running_screens_and_hands_the_signal_on(
    tmp_path: Path, database: Path
) -> None:
    seen: list[int] = []

    def earlier_handler(signum: int, frame: object) -> None:
        seen.append(signum)
        raise SystemExit(128 + signum)

    previous = signal.signal(signal.SIGTERM, earlier_handler)
    timer = threading.Timer(2.0, os.kill, (os.getpid(), signal.SIGTERM))
    timer.start()
    try:
        with pytest.raises(SystemExit) as exited:
            run_step(tmp_path, database, ["hang_a"])
        restored = signal.getsignal(signal.SIGTERM)
    finally:
        timer.cancel()
        signal.signal(signal.SIGTERM, previous)

    assert exited.value.code == 128 + signal.SIGTERM
    assert seen == [signal.SIGTERM]
    assert restored is earlier_handler
    assert wait_until_gone(child_of(tmp_path, "hang_a"))


def test_the_sigterm_handler_is_put_back_after_a_normal_run(
    tmp_path: Path, database: Path
) -> None:
    before = signal.getsignal(signal.SIGTERM)

    run_step(tmp_path, database, ["ok_a"])

    assert signal.getsignal(signal.SIGTERM) is before


def test_a_step_in_a_worker_thread_leaves_sigterm_alone(
    tmp_path: Path, database: Path
) -> None:
    before = signal.getsignal(signal.SIGTERM)
    errors: list[BaseException] = []

    def work() -> None:
        try:
            run_step(tmp_path, database, ["ok_a"])
        except BaseException as err:  # noqa: BLE001
            errors.append(err)

    thread = threading.Thread(target=work)
    thread.start()
    thread.join()

    assert errors == []
    assert signal.getsignal(signal.SIGTERM) is before


@pytest.mark.parametrize("mode", ["torn", "wronginput", "garbled"])
def test_a_summary_that_is_incomplete_or_for_other_bytes_is_not_accepted(
    tmp_path: Path, database: Path, caplog: pytest.LogCaptureFixture, mode: str
) -> None:
    with caplog.at_level(logging.WARNING, logger=LOGGER.name):
        results = run_step(tmp_path, database, [f"{mode}_a"])

    columns = screen_columns(results[0])
    assert columns["screen_status"] == "failed"
    assert math.isnan(columns["screen_cellular_pct"])
    assert math.isnan(columns["screen_flagged_pct"])
    assert caplog.records


DRIVER = """
import logging
import os
import signal
import sys
import threading
import time
from pathlib import Path

from src.pipeline import contig_screen_step as step

step.SCREEN_COMMAND = (sys.executable, sys.argv[1])
mode = os.environ.get("DRIVER_MODE", "plain")
work = Path(sys.argv[2])
query = work / "queries" / "hang_a.fna"
query.parent.mkdir(parents=True)
query.write_text(">c1\\nACGT\\n", encoding="utf-8")
manifest = work / "db" / "contig_screen" / "manifest.json"
manifest.parent.mkdir(parents=True)
manifest.write_text("{}", encoding="utf-8")
if mode == "slow_add":
    # A window of three seconds between the start of a screen and its registration.
    original_add = step._Running.add

    def slow_add(self, process):
        time.sleep(3)
        original_add(self, process)

    step._Running.add = slow_add
if mode == "stop_signal":
    # SIGTERM arrives while stop() holds the lock, on the path of a KeyboardInterrupt.
    original_init = step._Running.__init__

    class Trigger(set):
        def __iter__(self):
            os.kill(os.getpid(), signal.SIGTERM)
            return super().__iter__()

    def init(self):
        original_init(self)
        self._processes = Trigger()

    step._Running.__init__ = init

    def interrupt(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGALRM, interrupt)
    signal.alarm(2)
results = [{"query": "hang_a", "status": "complete", "summary_data": {}}]
options = step.ScreenOptions(True, False, False, lambda n: (1, 1))
try:
    step.add_screen_columns(
        results, {"hang_a": query}, work / "out", work / "db", options,
        logging.getLogger("driver"),
    )
except KeyboardInterrupt:
    print("KeyboardInterrupt")
"""


def start_driver(tmp_path: Path, mode: str) -> tuple[subprocess.Popen[bytes], Path]:
    """Run the step in a child process; return it and the pid file of its screen."""
    (tmp_path / "driver.py").write_text(DRIVER, encoding="utf-8")
    work = tmp_path / "work"
    root = Path(__file__).resolve().parents[1]
    runner = subprocess.Popen(
        [
            sys.executable,
            str(tmp_path / "driver.py"),
            str(tmp_path / "fake_screen.py"),
            str(work),
        ],
        cwd=root,
        env={**os.environ, "PYTHONPATH": str(root), "DRIVER_MODE": mode},
        stdout=subprocess.PIPE,
    )
    return runner, work / "out" / "contig_screen" / "hang_a" / "child.pid"


def wait_for_file(path: Path) -> int:
    deadline = time.monotonic() + 60
    while not path.exists() and time.monotonic() < deadline:
        time.sleep(0.2)
    return int(path.read_text())


def test_sigterm_with_the_default_handler_still_stops_the_screens_and_ends_the_run(
    tmp_path: Path,
) -> None:
    runner, pid_file = start_driver(tmp_path, "plain")
    child = wait_for_file(pid_file)

    runner.send_signal(signal.SIGTERM)

    assert runner.wait(timeout=30) == -signal.SIGTERM
    assert wait_until_gone(child)


def test_sigterm_between_the_start_of_a_screen_and_its_registration_orphans_nothing(
    tmp_path: Path,
) -> None:
    runner, pid_file = start_driver(tmp_path, "slow_add")
    child = wait_for_file(pid_file)

    runner.send_signal(signal.SIGTERM)

    assert runner.wait(timeout=60) == -signal.SIGTERM
    assert wait_until_gone(child)


def test_sigterm_during_the_stop_after_an_interrupt_does_not_deadlock(
    tmp_path: Path,
) -> None:
    runner, pid_file = start_driver(tmp_path, "stop_signal")
    child = wait_for_file(pid_file)

    try:
        code = runner.wait(timeout=30)
    except subprocess.TimeoutExpired:
        runner.kill()
        pytest.fail("the step hangs when SIGTERM arrives during stop()")

    assert code == 0
    assert b"KeyboardInterrupt" in runner.stdout.read()
    assert wait_until_gone(child)


COUNT_COLUMNS = ["input_bp", "cellular_bp", "giant_virus_bp", "flagged_bp"]
SAVED_ROW = {
    "bin": "ok_a",
    "cellular_status": "ok",
    "giant_virus_status": "ok",
    "input_bp": "1000",
    "cellular_bp": "70",
    "cellular_pct": "7.00",
    "giant_virus_bp": "80",
    "giant_virus_pct": "8.00",
    "flagged_bp": "90",
    "flagged_pct": "9.00",
    "manifest_sha256": MANIFEST_DIGEST,
    "input_sha256": INPUT_DIGEST,
}


def save_row(folder: Path, **changes: str | None) -> None:
    """Save a finished row; a value of None leaves that column out of the file."""
    row = {**SAVED_ROW, **changes}
    kept = {column: value for column, value in row.items() if value is not None}
    folder.mkdir(parents=True)
    (folder / "contig_screen_summary.tsv").write_text(
        "\t".join(kept) + "\n" + "\t".join(kept.values()) + "\n", encoding="utf-8"
    )


def test_a_saved_row_with_every_column_is_reused(
    tmp_path: Path, database: Path
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    save_row(folder)

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert screen_columns(results[0])["screen_flagged_pct"] == 9.0
    assert not (folder / "args.json").exists()


@pytest.mark.parametrize("column", COUNT_COLUMNS)
def test_resume_screens_again_when_a_saved_row_lacks_a_count_column(
    tmp_path: Path, database: Path, column: str
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    save_row(folder, **{column: None})

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert screen_columns(results[0])["screen_status"] == "ok"
    assert screen_columns(results[0])["screen_flagged_pct"] == 25.0
    assert (folder / "args.json").exists()


@pytest.mark.parametrize("column", COUNT_COLUMNS)
def test_resume_screens_again_when_a_saved_count_does_not_parse(
    tmp_path: Path, database: Path, column: str
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    save_row(folder, **{column: "oops"})

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert screen_columns(results[0])["screen_flagged_pct"] == 25.0
    assert (folder / "args.json").exists()


@pytest.mark.parametrize(
    ("cellular", "giant_virus", "estimate", "contamination_type"),
    [
        (10.0, 0.0, 10.0, "cellular"),
        (55.5, 80.0, 55.5, "cellular"),
        (100.0, 100.0, 100.0, "cellular"),
        (9.99, 10.0, 9.99, "mixed_viral"),
        (0.0, 42.0, 0.0, "mixed_viral"),
        (9.99, 9.99, 9.99, "clean"),
        (0.0, 0.0, 0.0, "clean"),
        (4.0, NAN, 4.0, "clean"),
        (12.5, NAN, 12.5, "cellular"),
        (12.3456, 0.0, 12.35, "cellular"),
    ],
)
def test_a_cellular_share_becomes_the_contamination_estimate_of_the_bin(
    cellular: float, giant_virus: float, estimate: float, contamination_type: str
) -> None:
    result = step.ScreenResult("ok", cellular, giant_virus, NAN)

    assert result.contamination_columns() == {
        "estimated_contamination": estimate,
        "estimated_contamination_strategy": "contig_screen_cellular_v1",
        "contamination_query_support": "supported",
        "contamination_model_reliability": "simulated_benchmark",
        "contamination_type": contamination_type,
    }


@pytest.mark.parametrize(
    "result",
    [
        step.ScreenResult("disabled"),
        step.ScreenResult("not_installed"),
        step.ScreenResult("not_applicable"),
        step.ScreenResult("failed"),
        step.ScreenResult("partial", giant_virus_pct=20.0),
        step.ScreenResult("partial", NAN, 0.0, NAN),
    ],
)
def test_a_bin_without_a_cellular_share_has_no_screen_contamination_columns(
    result: step.ScreenResult,
) -> None:
    assert result.contamination_columns() == {}


def test_a_scored_bin_replaces_all_five_marker_model_fields(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["ok_a"])

    assert contamination_columns(results[0]) == {
        "estimated_contamination": 10.0,
        "estimated_contamination_strategy": "contig_screen_cellular_v1",
        "contamination_query_support": "supported",
        "contamination_model_reliability": "simulated_benchmark",
        "contamination_type": "cellular",
    }


@pytest.mark.parametrize(
    ("row", "estimate", "contamination_type"),
    [
        ("ok_a\tok\tok\t1000\t70\t7.00\t150\t15.00\t200\t20.00\n", 7.0, "mixed_viral"),
        ("ok_a\tok\tok\t1000\t70\t7.00\t80\t8.00\t90\t9.00\n", 7.0, "clean"),
        ("ok_a\tok\tnot_scored\t1000\t70\t7.00\tNaN\tNaN\tNaN\tNaN\n", 7.0, "clean"),
        (
            "ok_a\tok\tnot_scored\t1000\t300\t30.00\tNaN\tNaN\tNaN\tNaN\n",
            30.0,
            "cellular",
        ),
    ],
)
def test_a_reused_screen_result_sets_the_contamination_of_the_bin_too(
    tmp_path: Path, database: Path, row: str, estimate: float, contamination_type: str
) -> None:
    folder = tmp_path / "out" / "contig_screen" / "ok_a"
    seed_summary(folder, row, MANIFEST_DIGEST)

    results = run_step(tmp_path, database, ["ok_a"], resume=True)

    assert not (folder / "args.json").exists()
    columns = contamination_columns(results[0])
    assert columns["estimated_contamination"] == estimate
    assert columns["contamination_type"] == contamination_type
    assert columns["estimated_contamination_strategy"] == "contig_screen_cellular_v1"


@pytest.mark.parametrize("name", ["crash_a", "dead_a", "short_a", "prot_a"])
def test_a_bin_without_a_cellular_share_keeps_the_marker_model_values(
    tmp_path: Path, database: Path, name: str
) -> None:
    results = run_step(tmp_path, database, [name])

    assert contamination_columns(results[0]) == MARKER_MODEL


def test_a_disabled_screen_keeps_the_marker_model_values(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["ok_a"], enabled=False)

    assert contamination_columns(results[0]) == MARKER_MODEL


def test_a_missing_addon_keeps_the_marker_model_values(tmp_path: Path) -> None:
    results = run_step(tmp_path, tmp_path / "db", ["ok_a"])

    assert contamination_columns(results[0]) == MARKER_MODEL


def test_the_bins_of_one_run_take_their_own_contamination(
    tmp_path: Path, database: Path
) -> None:
    results = run_step(tmp_path, database, ["ok_a", "crash_b", "single_c"])

    assert [r["summary_data"]["estimated_contamination"] for r in results] == [
        10.0,
        41.5,
        10.0,
    ]
