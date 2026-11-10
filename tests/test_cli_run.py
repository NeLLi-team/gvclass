"""A command-line run: exit status, the runner process and the temporary files.

The database setup, the runner and the output validation are replaced, so the tests
see what the command decides, not what the pipeline computes.
"""

from __future__ import annotations

import os
import re
import signal
import sys
import tempfile
import threading
from pathlib import Path
from unittest.mock import Mock

import pytest

from src.bin import gvclass_cli
from src.bin.gvclass_cli import (
    MONITOR_THREAD_NAME,
    THREAD_ENV_VARS,
    RunnerOutcome,
    run_runner,
)
from src.bin.progress_monitor import ResourceMonitor
from src.pipeline.run_config import plan_workers
from tests.cli_helpers import SUCCESS, run_cli

CONTIG_FASTA = ">c1\n" + "ACGT" * 3000 + "\n"
FAILED = RunnerOutcome(return_code=1, tail=(), interrupted=False, peak_memory_mb=0.0)
INTERRUPTED = RunnerOutcome(
    return_code=-15, tail=(), interrupted=True, peak_memory_mb=0.0
)


def test_the_banner_and_the_runner_cap_workers_at_the_number_of_queries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """``-j 8 -t 16`` with two queries: two workers with eight threads each."""
    run = run_cli(tmp_path, monkeypatch, "-j", "8", "-t", "16")

    config = run.config
    assert run.exit_code == 0
    assert re.search(r"Workers: 2 \S 8 threads", capsys.readouterr().out)
    assert plan_workers(
        2, config.threads, config.max_workers, config.threads_per_worker
    ) == (2, 8)


def test_the_runner_gets_the_planned_per_marker_threads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    for name in THREAD_ENV_VARS:
        monkeypatch.delenv(name, raising=False)

    run = run_cli(tmp_path, monkeypatch, "-t", "16", "-j", "8")

    env = run.runner_calls[0]["env"]
    assert {env[name] for name in THREAD_ENV_VARS} == {"4"}


def test_the_command_starts_the_runner_on_the_written_run_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(tmp_path, monkeypatch, "-v")

    [call] = run.runner_calls
    assert call["command"] == [
        sys.executable,
        "-m",
        "src.bin.gvclass_runner",
        str(run.output_dir / "run_config.json"),
        "--verbose",
    ]
    assert call["verbose"] is True
    assert run.config.query_path == tmp_path / "queries"


def test_the_species_tree_flags_reach_the_run_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    on = run_cli(
        tmp_path,
        monkeypatch,
        "--species-tree-combined",
        "--species-tree-trim",
        "pytrimal",
    ).config

    assert (on.species_tree, on.species_tree_combined) == (True, True)
    assert on.species_tree_trim == "pytrimal"

    off = run_cli(tmp_path / "again", monkeypatch).config

    assert (off.species_tree, off.species_tree_combined) == (False, False)
    assert off.species_tree_trim == "witchi"


def test_the_exit_status_is_0_when_the_runner_and_the_validation_succeed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert run_cli(tmp_path, monkeypatch).exit_code == 0


def test_the_exit_status_is_1_and_the_issues_are_printed_when_validation_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    validation = {"success": False, "issues": ["q0: archive missing"], "warnings": []}

    run = run_cli(tmp_path, monkeypatch, validation=validation)

    assert run.exit_code == 1
    assert "q0: archive missing" in capsys.readouterr().out


def test_a_failed_runner_gives_exit_status_1_and_its_last_lines(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    outcome = RunnerOutcome(
        return_code=3, tail=("last line\n",), interrupted=False, peak_memory_mb=1.0
    )

    run = run_cli(tmp_path, monkeypatch, outcome=outcome)

    printed = capsys.readouterr().out
    assert run.exit_code == 1
    assert "exit code: 3" in printed
    assert "last line" in printed


def test_a_single_file_without_contigs_is_an_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    source = tmp_path / "assembly.fna"
    source.write_text(">c\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(
        gvclass_cli, "check_and_setup_database", Mock(return_value=True)
    )
    argv = [str(source), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]

    code = gvclass_cli.main(argv)

    assert code == 1
    assert "--contigs" in capsys.readouterr().err


def test_the_run_config_carries_the_effective_min_length(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(
        tmp_path,
        monkeypatch,
        "--contigs",
        "--contigs-min-length",
        "5000",
        "--min-length",
        "30000",
        suffix=".fna",
        body=CONTIG_FASTA,
    )

    assert run.config.min_sequence_length == 5000
    assert run.config.contigs_min_length == 5000
    assert run.config.contigs_mode is True


@pytest.mark.parametrize(
    "outcome", [SUCCESS, FAILED, INTERRUPTED], ids=["success", "failure", "interrupt"]
)
def test_the_temporary_contig_files_are_removed_when_the_runner_ends(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, outcome: RunnerOutcome
) -> None:
    run = run_cli(
        tmp_path,
        monkeypatch,
        "--contigs",
        suffix=".fna",
        body=CONTIG_FASTA,
        outcome=outcome,
    )

    assert run.config.contigs_mode is True
    assert not run.config.query_path.exists()


def test_the_temporary_contig_files_are_removed_when_the_runner_cannot_start(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(
        tmp_path,
        monkeypatch,
        "--contigs",
        suffix=".fna",
        body=CONTIG_FASTA,
        runner_error=OSError("no such interpreter"),
    )

    assert run.exit_code == 1
    assert not run.config.query_path.exists()


def test_the_temporary_contig_files_are_removed_when_sigterm_comes_before_the_runner(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", None)
    monkeypatch.setenv("TMPDIR", str(scratch))
    handler_before = signal.getsignal(signal.SIGTERM)

    def sigterm_during_the_database_check(*_args: object, **_kwargs: object) -> bool:
        assert any(scratch.iterdir())
        signal.raise_signal(signal.SIGTERM)
        return True

    run = run_cli(
        tmp_path,
        monkeypatch,
        "--contigs",
        suffix=".fna",
        body=CONTIG_FASTA,
        database_check=sigterm_during_the_database_check,
    )

    assert run.exit_code == 1
    assert run.runner_calls == []
    assert list(scratch.iterdir()) == []
    assert signal.getsignal(signal.SIGTERM) == handler_before


def _is_running(pid: int) -> bool:
    """Tell whether process ``pid`` is alive, a zombie not counted."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
    except FileNotFoundError:
        return False
    return stat.rsplit(")", 1)[1].split()[0] != "Z"


CHILD_WITH_GRANDCHILD = """
import os, subprocess, sys, time
grandchild = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
print(f"PIDS {os.getpid()} {grandchild.pid}", flush=True)
time.sleep(60)
"""


@pytest.mark.parametrize("stop", ["interrupt", "sigterm"])
def test_run_runner_stops_the_process_group_and_the_monitor_on_interrupt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, stop: str
) -> None:
    seen: dict = {}
    handler_before = signal.getsignal(signal.SIGTERM)

    def interrupt_on_pids(_monitor: ResourceMonitor, line: str) -> bool:
        if not line.startswith("PIDS"):
            return False
        seen["pids"] = [int(pid) for pid in line.split()[1:]]
        seen["monitor"] = [
            thread
            for thread in threading.enumerate()
            if thread.name == MONITOR_THREAD_NAME
        ]
        if stop == "sigterm":
            signal.raise_signal(signal.SIGTERM)
        raise KeyboardInterrupt

    monkeypatch.setattr(ResourceMonitor, "handle_progress_line", interrupt_on_pids)

    outcome = run_runner(
        [sys.executable, "-c", CHILD_WITH_GRANDCHILD],
        env=os.environ.copy(),
        cwd=tmp_path,
        output_dir=tmp_path,
        queries_to_process=1,
        verbose=False,
    )

    assert outcome.interrupted is True
    assert not any(_is_running(pid) for pid in seen["pids"])
    assert [thread.daemon for thread in seen["monitor"]] == [True]
    assert MONITOR_THREAD_NAME not in {t.name for t in threading.enumerate()}
    assert signal.getsignal(signal.SIGTERM) == handler_before


def test_run_runner_keeps_the_last_200_lines_that_are_not_progress(
    tmp_path: Path,
) -> None:
    code = (
        "import sys\n"
        "for i in range(300): print(f'line {i}', file=sys.stderr)\n"
        "print('GVCLASS_PROGRESS {}', file=sys.stderr)\n"
        "sys.exit(3)\n"
    )

    outcome = run_runner(
        [sys.executable, "-c", code],
        env=os.environ.copy(),
        cwd=tmp_path,
        output_dir=tmp_path,
        queries_to_process=1,
        verbose=False,
    )

    assert outcome.return_code == 3
    assert outcome.interrupted is False
    assert len(outcome.tail) == 200
    assert outcome.tail[-1] == "line 299\n"
    assert not any("GVCLASS_PROGRESS" in line for line in outcome.tail)
