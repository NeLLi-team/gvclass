"""The completion record and the status files of a run.

``results.jsonl`` is the one record of finished queries: an append-only file with one
JSON line per finished query, flushed and fsynced, so a crash leaves every earlier
line readable. The last line of a query wins. ``--resume`` reads it once and reuses
a query only if nothing that shaped its result has changed.

``run_status.json`` describes the run and is written twice, when the run starts and
when it ends. ``run.log`` keeps one event line per query event.
"""

from __future__ import annotations

import json
import logging
import os
import threading
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from src.__version__ import __version__
from src.pipeline.run_config import RunConfig
from src.pipeline.summary_writer import final_summary_record
from src.utils.common import sha256_file

logger = logging.getLogger("gvclass_runner")

RESULTS_FILENAME = "results.jsonl"
RUN_STATUS_FILENAME = "run_status.json"
RUN_LOG_FILENAME = "run.log"
RUN_STATUS_VERSION = 2
SOFTWARE_VERSION = f"v{__version__}"

COMPLETE = "complete"
FAILED = "failed"
INVALID_INPUT = "invalid_input"


@dataclass(frozen=True, slots=True)
class QueryRecord:
    """One line of ``results.jsonl``: how one query finished.

    ``settings`` are the result-changing settings of the run. ``archive`` (file name,
    size, SHA-256) and ``summary`` (the final-schema row) exist only for a
    ``complete`` query.
    """

    query: str
    status: str
    error: str | None
    input_sha256: str | None
    settings: dict[str, Any]
    database_version: str
    software_version: str
    finished_at: str
    archive: dict[str, Any] | None
    summary: dict[str, str] | None


def utc_now() -> str:
    """Return the current UTC time as an ISO 8601 string."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_results(output_dir: Path) -> dict[str, QueryRecord]:
    """Return the last record of each query in ``results.jsonl``.

    A line that cannot be read, such as the torn last line of a killed run, is
    skipped with a warning.
    """
    path = output_dir / RESULTS_FILENAME
    if not path.is_file():
        return {}
    records: dict[str, QueryRecord] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            record = QueryRecord(**json.loads(line))
        except (TypeError, ValueError):
            logger.warning("Ignoring unreadable line %d of %s", number, path)
            continue
        records[record.query] = record
    return records


def archive_is_intact(output_dir: Path, record: QueryRecord) -> bool:
    """Return True if the archive of ``record`` exists with the recorded size."""
    if record.archive is None:
        return False
    path = output_dir / record.archive["filename"]
    return path.is_file() and path.stat().st_size == record.archive["bytes"]


def resumable_queries(
    config: RunConfig, query_files: Iterable[Path], database_version: str
) -> dict[str, QueryRecord]:
    """Return the records of the queries that ``--resume`` may skip.

    A query is skipped only if its last record is ``complete``; the SHA-256 of its
    input, the result settings, the database version and the software version equal
    the current ones; and its archive exists with the recorded size. The others are
    reprocessed, and the log names every field that differed.
    """
    records = read_results(config.output_dir)
    settings = config.result_settings()
    reusable: dict[str, QueryRecord] = {}
    for query_file in query_files:
        name = query_file.stem
        record = records.get(name)
        differences = _differences(
            record, query_file, config.output_dir, settings, database_version
        )
        if differences:
            logger.info("Resume: reprocessing %s (%s)", name, "; ".join(differences))
        elif record is not None:
            reusable[name] = record
    return reusable


def _differences(
    record: QueryRecord | None,
    query_file: Path,
    output_dir: Path,
    settings: Mapping[str, Any],
    database_version: str,
) -> list[str]:
    """Return what differs between the saved record and the run now; empty if none."""
    if record is None:
        return ["no earlier record"]
    if record.status != COMPLETE:
        return [f"last record is {record.status}"]
    differences = []
    if record.input_sha256 != sha256_file(query_file):
        differences.append("input_sha256 changed")
    differences.extend(
        f"{key} changed"
        for key in sorted(settings.keys() | record.settings.keys())
        if settings.get(key) != record.settings.get(key)
    )
    if record.database_version != database_version:
        differences.append(
            f"database_version changed from {record.database_version} "
            f"to {database_version}"
        )
    if record.software_version != SOFTWARE_VERSION:
        differences.append(
            f"software_version changed from {record.software_version} "
            f"to {SOFTWARE_VERSION}"
        )
    if not archive_is_intact(output_dir, record):
        differences.append("archive missing or changed in size")
    return differences


class RunRecorder:
    """Thread-safe writer for ``results.jsonl``, ``run_status.json`` and ``run.log``."""

    def __init__(
        self,
        config: RunConfig,
        database_version: str,
        query_names: Iterable[str],
        skipped: Mapping[str, QueryRecord],
    ) -> None:
        self.config = config
        self.database_version = database_version
        self.skipped = dict(skipped)
        self._output_dir = config.output_dir
        self._lock = threading.Lock()
        self._started_at = utc_now()
        self._finished = False
        self._queries: dict[str, dict[str, Any]] = {
            name: {"status": "pending", "resume_skipped": False} for name in query_names
        }
        for name, record in self.skipped.items():
            self._queries[name] = {
                "status": COMPLETE,
                "resume_skipped": True,
                "finished_at": record.finished_at,
            }

    def start(self) -> None:
        """Write the first ``run_status.json`` and the start lines of ``run.log``."""
        lines = [
            (
                f"{self._started_at}\tRUN_START\tsoftware={SOFTWARE_VERSION}\t"
                f"database={self.database_version}\tqueries={len(self._queries)}\n"
            )
        ]
        lines.extend(
            f"{self._started_at}\tQUERY_RESUME_SKIPPED\t{name}\n"
            for name in sorted(self.skipped)
        )
        with self._lock:
            self._output_dir.mkdir(parents=True, exist_ok=True)
            self._end_torn_results_line()
            mode = "a" if self.config.resume else "w"
            with (self._output_dir / RUN_LOG_FILENAME).open(mode) as log:
                log.write("".join(lines))
            self._write_status("running")

    def record_query_started(self, query_name: str) -> None:
        """Log that a query started; the status file is not rewritten."""
        now = utc_now()
        with self._lock:
            self._queries[query_name]["status"] = "running"
            self._log(f"{now}\tQUERY_START\t{query_name}\n")

    def record_complete(self, query_file: Path, result: Mapping[str, Any]) -> None:
        """Append the ``complete`` record of a query whose archive is published."""
        name = query_file.stem
        archive = self._output_dir / f"{name}.tar.gz"
        record = self._record(
            name,
            COMPLETE,
            input_sha256=sha256_file(query_file),
            archive={
                "filename": archive.name,
                "bytes": archive.stat().st_size,
                "sha256": sha256_file(archive),
            },
            summary=final_summary_record(result["summary_data"]),
        )
        self._append(record, f"QUERY_COMPLETE\t{name}\tarchive={archive.name}")

    def record_failed(self, query_name: str, error: str, status: str = FAILED) -> None:
        """Append the record of a query that failed or was rejected as input."""
        record = self._record(query_name, status, error=error)
        event = "QUERY_FAILED" if status == FAILED else "QUERY_REJECTED"
        self._append(record, f"{event}\t{query_name}\t{error}")

    def finish(self, status: str) -> None:
        """Write the final ``run_status.json``; later calls do nothing."""
        with self._lock:
            if self._finished:
                return
            self._finished = True
            self._write_status(status, completed_at=utc_now())
            self._log(f"{utc_now()}\tRUN_{status.upper()}\n")

    def _record(
        self,
        query_name: str,
        status: str,
        *,
        error: str | None = None,
        input_sha256: str | None = None,
        archive: dict[str, Any] | None = None,
        summary: dict[str, str] | None = None,
    ) -> QueryRecord:
        return QueryRecord(
            query=query_name,
            status=status,
            error=error,
            input_sha256=input_sha256,
            settings=self.config.result_settings(),
            database_version=self.database_version,
            software_version=SOFTWARE_VERSION,
            finished_at=utc_now(),
            archive=archive,
            summary=summary,
        )

    def _append(self, record: QueryRecord, event: str) -> None:
        """Append ``record`` to ``results.jsonl``, fsynced, and log ``event``."""
        line = json.dumps(asdict(record), sort_keys=True) + "\n"
        with self._lock:
            with (self._output_dir / RESULTS_FILENAME).open(
                "a", encoding="utf-8"
            ) as handle:
                handle.write(line)
                handle.flush()
                os.fsync(handle.fileno())
            self._queries[record.query] = {
                "status": record.status,
                "resume_skipped": False,
                "finished_at": record.finished_at,
                **({"error": record.error} if record.error else {}),
            }
            self._log(f"{record.finished_at}\t{event}\n")

    def _end_torn_results_line(self) -> None:
        """Start a new line if a killed run left a partial last line."""
        path = self._output_dir / RESULTS_FILENAME
        if path.is_file() and path.stat().st_size > 0:
            with path.open("rb+") as handle:
                handle.seek(-1, os.SEEK_END)
                if handle.read(1) != b"\n":
                    handle.write(b"\n")

    def _log(self, line: str) -> None:
        with (self._output_dir / RUN_LOG_FILENAME).open("a") as log:
            log.write(line)

    def _write_status(self, status: str, completed_at: str | None = None) -> None:
        """Replace ``run_status.json`` atomically. The caller holds the lock."""
        queries = self._queries
        payload = {
            "schema_version": RUN_STATUS_VERSION,
            "run_id": self._started_at.replace(":", "").replace("+00:00", "Z"),
            "status": status,
            "started_at": self._started_at,
            "completed_at": completed_at,
            "software_version": SOFTWARE_VERSION,
            "database": {
                "path": str(self.config.database_path),
                "version": self.database_version,
            },
            "settings": self.config.to_dict(),
            "totals": _totals(queries),
            "queries": queries,
        }
        path = self._output_dir / RUN_STATUS_FILENAME
        temporary = path.with_suffix(".json.part")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)


def _totals(queries: Mapping[str, Mapping[str, Any]]) -> dict[str, int]:
    states = [query["status"] for query in queries.values()]
    return {
        "queries": len(queries),
        "complete": states.count(COMPLETE),
        "failed": states.count(FAILED) + states.count(INVALID_INPUT),
        "resume_skipped": sum(query["resume_skipped"] for query in queries.values()),
    }
