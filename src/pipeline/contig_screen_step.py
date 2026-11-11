"""Screen the bins of a run and add the totals to their summary rows.

The step runs after the queries and before the final summary, so it also covers the
queries that ``--resume`` skipped. Every bin that needs a screen gets one ``python -m
src.bin.gvclass_screen`` subprocess, run by a small thread pool whose size and threads
per bin come from the worker plan of the run for that many bins. The subprocess keeps
torch and its process-wide thread limit out of the process that runs the queries in a
thread pool. A screen failure sets ``screen_status`` and logs a warning; it never
fails a query or the run. A screen gets its own session and a time limit
(``SCREEN_TIMEOUT_S``); on a timeout, an interrupt or SIGTERM the whole process group is
killed, so a DIAMOND child dies too. Its output goes to ``screen.log`` in the bin
folder, not to memory. The files of a bin go to ``<output>/contig_screen/<query>/``,
outside the query folder, because query folders are archived into a tar. Under
``--resume`` a saved result is reused only if no screen failed, it names the same
bin, and it carries the SHA-256 of the bin file and of the manifest of the add-on
that is installed now.

A bin that the cellular screen scored also takes the cellular share as its
``estimated_contamination``, with the strategy ``contig_screen_cellular_v1`` (see
``ScreenResult.contamination_columns``). Every other bin keeps the values of the
marker model, because missing evidence is not zero.
"""

from __future__ import annotations

import contextlib
import logging
import math
import os
import signal
import subprocess
import sys
import threading
from collections.abc import Callable, Iterator, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any

import pandas as pd

from src.contig_screen import layout
from src.contig_screen.tables import read_table
from src.pipeline.progress_events import emit_progress_event
from src.pipeline.query_processing_engine import is_nucleotide_input
from src.utils.common import sha256_file

SCREEN_COMMAND = (sys.executable, "-m", "src.bin.gvclass_screen")
REPO_DIR = Path(__file__).resolve().parents[2]
SUMMARY_FILE = "contig_screen_summary.tsv"
INSTALL_COMMAND = "gvclass screen --install"
TEXT_COLUMNS = ["bin", "manifest_sha256", "input_sha256"]
NUMBER_COLUMNS = (
    "input_bp",
    "cellular_bp",
    "cellular_pct",
    "giant_virus_bp",
    "giant_virus_pct",
    "flagged_bp",
    "flagged_pct",
)
#: The columns of a saved row. A row that lacks one is not read: it has no usable
#: numbers, so a ``--resume`` screens the bin again.
SUMMARY_COLUMNS = frozenset(
    {"cellular_status", "giant_virus_status", *TEXT_COLUMNS, *NUMBER_COLUMNS}
)
SCREEN_STATUSES = frozenset({"ok", "failed", "not_scored"})
STDERR_TAIL_LINES = 3
LOG_TAIL_BYTES = 4096
LOG_LINE_CHARS = 300
LOG_FILE = "screen.log"
#: A screen that runs longer than this is stopped and the bin is ``failed``. A bin of
#: 176 kb takes about 2 minutes, so a run of 4 hours is stuck.
SCREEN_TIMEOUT_S = 4 * 60 * 60

OK = "ok"
PARTIAL = "partial"
FAILED = "failed"
DISABLED = "disabled"
NOT_INSTALLED = "not_installed"
NOT_APPLICABLE = "not_applicable"

#: A cellular or giant-virus share at or above this percent is a contamination call.
CONTAMINATION_CALL_PCT = 10.0
CONTAMINATION_STRATEGY = "contig_screen_cellular_v1"
CONTAMINATION_SUPPORT = "supported"
CONTAMINATION_RELIABILITY = "simulated_benchmark"
CONTAMINATION_CELLULAR = "cellular"
CONTAMINATION_MIXED_VIRAL = "mixed_viral"
CONTAMINATION_CLEAN = "clean"


#: Maps the number of bins to screen to (bins at once, threads for each).
WorkerPlan = Callable[[int], tuple[int, int]]


@dataclass(frozen=True, slots=True)
class ScreenOptions:
    """What the run asks of the step.

    ``plan_workers`` sizes the pool from the number of bins that need a screen run,
    so bins that are skipped or reused under ``--resume`` take no threads.
    """

    enabled: bool
    contigs_mode: bool
    resume: bool
    plan_workers: WorkerPlan


@dataclass(frozen=True, slots=True)
class ScreenResult:
    """The status of one bin and the shares of its base pairs that the screens call."""

    status: str
    cellular_pct: float = math.nan
    giant_virus_pct: float = math.nan
    flagged_pct: float = math.nan

    def columns(self) -> dict[str, str | float]:
        """Return the ``screen_*`` columns of the summary row."""
        return {
            "screen_status": self.status,
            "screen_cellular_pct": self.cellular_pct,
            "screen_giant_virus_pct": self.giant_virus_pct,
            "screen_flagged_pct": self.flagged_pct,
        }

    def contamination_columns(self) -> dict[str, str | float]:
        """Return the contamination columns that the cellular share sets.

        The cellular share is the contamination estimate; it counts no second giant
        virus, so a share of the giant-virus screen at the call threshold gives the type
        ``mixed_viral`` when the cellular share is below it. The result is empty when
        the cellular screen did not score the bin: missing evidence is not zero, so the
        caller keeps the values of the marker model.
        """
        if not math.isfinite(self.cellular_pct):
            return {}
        estimate = round(self.cellular_pct, 2)
        # A NaN share of the giant-virus screen compares false, so it is not mixed.
        if estimate >= CONTAMINATION_CALL_PCT:
            contamination_type = CONTAMINATION_CELLULAR
        elif self.giant_virus_pct >= CONTAMINATION_CALL_PCT:
            contamination_type = CONTAMINATION_MIXED_VIRAL
        else:
            contamination_type = CONTAMINATION_CLEAN
        return {
            "estimated_contamination": estimate,
            "estimated_contamination_strategy": CONTAMINATION_STRATEGY,
            "contamination_query_support": CONTAMINATION_SUPPORT,
            "contamination_model_reliability": CONTAMINATION_RELIABILITY,
            "contamination_type": contamination_type,
        }


class _Running:
    """The screen processes that run now, and whether the step is stopping.

    ``stop`` sets the flag first; a screen that registers after it is killed at
    once, so none escapes between the snapshot and the shutdown of the pool.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._processes: set[subprocess.Popen[bytes]] = set()
        self.stopping = False

    def add(self, process: subprocess.Popen[bytes]) -> None:
        with self._lock:
            self._processes.add(process)
            late = self.stopping
        if late:
            _kill_group(process)

    def discard(self, process: subprocess.Popen[bytes]) -> None:
        with self._lock:
            self._processes.discard(process)

    def stop(self) -> None:
        with self._lock:
            self.stopping = True
            processes = list(self._processes)
        for process in processes:
            _kill_group(process)


@dataclass(frozen=True, slots=True)
class _Context:
    """What every bin of one step shares."""

    query_files: Mapping[str, Path]
    screen_dir: Path
    database: Path
    options: ScreenOptions
    manifest_sha256: str | None
    manifest_error: str | None
    running: _Running
    logger: logging.Logger

    @property
    def installed(self) -> bool:
        """True if the add-on is installed."""
        return self.manifest_sha256 is not None


def add_screen_columns(
    results: list[dict[str, Any]],
    query_files: Mapping[str, Path],
    output_dir: Path,
    database: Path,
    options: ScreenOptions,
    logger: logging.Logger,
) -> None:
    """Set the ``screen_*`` columns in the summary data of every completed query.

    A bin that the cellular screen scored also gets its contamination columns from the
    screen. A failed query has no summary row and is not screened.
    """
    rows = [r for r in results if r.get("status") == "complete" and "summary_data" in r]
    if not rows:
        return
    digest, manifest_error = _read_manifest_digest(database)
    context = _Context(
        query_files,
        output_dir / layout.SCREEN_DIR,
        database,
        options,
        digest,
        manifest_error,
        _Running(),
        logger,
    )
    if options.enabled and manifest_error is not None:
        logger.warning(
            "Contig screen add-on manifest in %s is not readable: %s",
            database,
            manifest_error,
        )
    elif options.enabled and not context.installed:
        logger.info(
            "Contig screen add-on not installed in %s; run `%s` to add the "
            "screen_* columns",
            database,
            INSTALL_COMMAND,
        )
    outcomes: dict[str, ScreenResult] = {}
    pending: list[dict[str, Any]] = []
    for row in rows:
        decided = _guarded(row["query"], logger, partial(_decide, row, context))
        if decided is None:
            pending.append(row)
        else:
            outcomes[row["query"]] = decided
    if pending:
        outcomes.update(_run_pending(pending, context))
    for row in rows:
        outcome = outcomes[row["query"]]
        row["summary_data"].update(outcome.columns())
        row["summary_data"].update(outcome.contamination_columns())


def _read_manifest_digest(database: Path) -> tuple[str | None, str | None]:
    """Return the digest of the add-on manifest and the error that blocked it."""
    try:
        return layout.manifest_digest(database), None
    except OSError as err:
        return None, str(err)


class _Terminated(BaseException):
    """SIGTERM arrived. Raised by the handler, which does nothing else."""


def _run_pending(
    pending: list[dict[str, Any]], context: _Context
) -> dict[str, ScreenResult]:
    """Run the screens of the bins in a pool; an interrupt or SIGTERM stops them."""
    workers, threads = context.options.plan_workers(len(pending))
    run = partial(_run_row, context=context, threads=threads)
    pool_size = max(1, min(workers, len(pending)))
    with _ended_by_sigterm(), ThreadPoolExecutor(pool_size) as pool:
        try:
            ran = list(pool.map(run, pending))
        except BaseException:
            try:
                _ignore_sigterm()  # no signal re-enters the cleanup
            finally:
                context.running.stop()
                pool.shutdown(wait=False, cancel_futures=True)
            raise
    return dict(zip((row["query"] for row in pending), ran, strict=True))


def _ignore_sigterm() -> None:
    if threading.current_thread() is threading.main_thread():
        signal.signal(signal.SIGTERM, signal.SIG_IGN)


@contextlib.contextmanager
def _ended_by_sigterm() -> Iterator[None]:
    """In the main thread, turn SIGTERM into ``_Terminated`` and end the run after it.

    The screens lead their own sessions, so the signal does not reach them. The
    handler only ignores further SIGTERM and raises; the caller stops the screens and
    leaves its pool, so every worker has finished. Then the old handler returns and
    the signal is sent again, so the run ends as it would have. Other threads cannot
    set handlers.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return

    def on_sigterm(_signum: int, _frame: object) -> None:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        raise _Terminated

    handler = signal.signal(signal.SIGTERM, on_sigterm)
    earlier = signal.SIG_DFL if handler is None else handler
    try:
        yield
    except _Terminated:
        signal.signal(signal.SIGTERM, earlier)
        os.kill(os.getpid(), signal.SIGTERM)
        raise SystemExit(128 + signal.SIGTERM) from None
    finally:
        signal.signal(signal.SIGTERM, earlier)


def screen_status(cellular: str, giant_virus: str) -> str:
    """Return the status of a bin from the statuses of its two screens.

    ``ok`` needs both screens to score the bin; ``failed`` means no screen scored it
    and one failed; anything else is ``partial``.
    """
    statuses = (cellular, giant_virus)
    if statuses == (OK, OK):
        return OK
    if OK not in statuses and FAILED in statuses:
        return FAILED
    return PARTIAL


def _guarded(
    name: str,
    logger: logging.Logger,
    call: Callable[[], ScreenResult | None],
) -> ScreenResult | None:
    """Return ``call()``; an error in the step fails only this bin."""
    try:
        return call()
    except Exception as err:  # noqa: BLE001  # the screen is optional: it never fails the run
        logger.warning("Contig screen failed for %s: %s", name, err)
        return ScreenResult(FAILED)


def _decide(row: dict[str, Any], context: _Context) -> ScreenResult | None:
    """Return the result of a bin that needs no run, or None when it must run."""
    query_file = context.query_files[row["query"]]
    skipped = _skip_status(query_file, context)
    if skipped is not None:
        return ScreenResult(skipped)
    if context.options.resume:
        name = row["query"]
        reused = _read_summary_row(context.screen_dir / name)
        if reused is not None and _is_reusable(reused, name, query_file, context):
            return _result_from_row(reused)
    return None


def _skip_status(query_file: Path, context: _Context) -> str | None:
    """Return why the bin is not screened, or None when it is."""
    if not context.options.enabled:
        return DISABLED
    if context.options.contigs_mode or not is_nucleotide_input(query_file):
        return NOT_APPLICABLE
    if context.manifest_error is not None:
        return FAILED
    if not context.installed:
        return NOT_INSTALLED
    return None


def _run_row(row: dict[str, Any], context: _Context, threads: int) -> ScreenResult:
    """Run the screen on the bin of one row; a failure is a ``failed`` result."""
    name = row["query"]
    result = _guarded(
        name,
        context.logger,
        partial(
            _run_screen,
            name,
            context.query_files[name],
            context.screen_dir / name,
            context,
            threads,
        ),
    )
    return result or ScreenResult(FAILED)


def _run_screen(
    name: str, query_file: Path, folder: Path, context: _Context, threads: int
) -> ScreenResult:
    """Run the screen on one bin and return its result."""
    # The queries are at 100% by now; the stage keeps the display from looking stuck.
    emit_progress_event(
        {
            "event": "contig_screen",
            "query": name,
            "progress": 100,
            "stage": "contig_screen",
        }
    )
    if context.running.stopping:
        return ScreenResult(FAILED)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    command = [
        *SCREEN_COMMAND,
        str(query_file),
        "-o",
        str(folder),
        "-d",
        str(context.database),
        "-t",
        str(threads),
    ]
    log_path = folder / LOG_FILE
    with log_path.open("wb") as log:
        process = subprocess.Popen(  # a fixed command with checked paths
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=REPO_DIR,
            start_new_session=True,
        )
        context.running.add(process)
        try:
            returncode = process.wait(timeout=SCREEN_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            _kill_group(process)
            (folder / SUMMARY_FILE).unlink(missing_ok=True)
            context.logger.warning(
                "Contig screen of %s timed out after %d s and was stopped: %s",
                name,
                SCREEN_TIMEOUT_S,
                _log_tail(log_path),
            )
            return ScreenResult(FAILED)
        finally:
            context.running.discard(process)
    if returncode != 0:
        _kill_group(process)  # what the crashed screen left running, DIAMOND too
    return _result_after_exit(name, query_file, folder, returncode, context)


def _result_after_exit(
    name: str, query_file: Path, folder: Path, returncode: int, context: _Context
) -> ScreenResult:
    """Return the result of a screen that ended; log what its log file ends with."""
    row = _read_summary_row(folder)
    if row is None or not _is_valid(row, name, query_file, context):
        context.logger.warning(
            "Contig screen of %s failed (exit code %d, summary missing, incomplete "
            "or not of this bin file): %s",
            name,
            returncode,
            _log_tail(folder / LOG_FILE),
        )
        return ScreenResult(FAILED)
    if returncode != 0:
        context.logger.warning(
            "Contig screen of %s ended with exit code %d: %s",
            name,
            returncode,
            _log_tail(folder / LOG_FILE),
        )
    return _result_from_row(row)


def _kill_group(process: subprocess.Popen[bytes]) -> None:
    """Kill a screen's process group: the screen and everything it started.

    The screen leads its own session, so its group id is its pid and a DIAMOND child
    is in the group. The group is signalled also after the leader has gone, because
    a child can outlive it. A group that is already empty is not an error.
    """
    with contextlib.suppress(ProcessLookupError, PermissionError):
        os.killpg(process.pid, signal.SIGKILL)
    process.wait()


def _read_summary_row(folder: Path) -> pd.Series | None:
    """Return the row of ``contig_screen_summary.tsv`` if it has the totals."""
    path = folder / SUMMARY_FILE
    if not path.is_file():
        return None
    try:
        table = read_table(path, TEXT_COLUMNS)
    except (OSError, ValueError):
        return None
    if table.empty or not set(table.columns) >= SUMMARY_COLUMNS:
        return None
    return table.iloc[0]


def _is_valid(row: pd.Series, name: str, query_file: Path, context: _Context) -> bool:
    """Return True if a summary row is whole and made for this bin and this add-on.

    The row names the bin; it carries the SHA-256 of the add-on manifest and of the
    bin file as they are now (an empty or missing digest never matches); every number
    parses; and a screen with status ``ok`` has its percentage.
    """
    if row["bin"] != name or row["manifest_sha256"] != context.manifest_sha256:
        return False
    if row["input_sha256"] != sha256_file(query_file):
        return False
    try:
        numbers = {column: float(row[column]) for column in NUMBER_COLUMNS}
    except (TypeError, ValueError):
        return False
    statuses = (str(row["cellular_status"]), str(row["giant_virus_status"]))
    scored = (
        (statuses[0] == OK, numbers["cellular_pct"]),
        (statuses[1] == OK, numbers["giant_virus_pct"]),
        (statuses == (OK, OK), numbers["flagged_pct"]),
    )
    return set(statuses) <= SCREEN_STATUSES and all(
        math.isfinite(value) for is_ok, value in scored if is_ok
    )


def _is_reusable(
    row: pd.Series, name: str, query_file: Path, context: _Context
) -> bool:
    """Return True if a saved result may stand under ``--resume``: valid, no failure."""
    failed = FAILED in (row["cellular_status"], row["giant_virus_status"])
    return not failed and _is_valid(row, name, query_file, context)


def _result_from_row(row: pd.Series) -> ScreenResult:
    return ScreenResult(
        screen_status(str(row["cellular_status"]), str(row["giant_virus_status"])),
        float(row["cellular_pct"]),
        float(row["giant_virus_pct"]),
        float(row["flagged_pct"]),
    )


def _log_tail(path: Path) -> str:
    """Return the last lines of a log file, read from its end, each cut short."""
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            handle.seek(max(0, handle.tell() - LOG_TAIL_BYTES))
            text = handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""
    lines = [line.strip()[-LOG_LINE_CHARS:] for line in text.splitlines()]
    return " | ".join([line for line in lines if line][-STDERR_TAIL_LINES:])
