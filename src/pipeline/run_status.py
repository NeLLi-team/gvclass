"""Single-file run status and resume metadata for GVClass outputs."""

from __future__ import annotations

import json
import os
import hashlib
import tarfile
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from src.__version__ import __version__ as _GVCLASS_VERSION
from src.utils.common import sha256_file
from src.utils.database_manager import DatabaseManager

RUN_STATUS_FILENAME = "run_status.json"
RUN_LOG_FILENAME = "run.log"
RUN_STATUS_VERSION = 1
_SOFTWARE_VERSION = f"v{_GVCLASS_VERSION}"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def run_status_path(output_dir: Path) -> Path:
    return output_dir / RUN_STATUS_FILENAME


def run_log_path(output_dir: Path) -> Path:
    return output_dir / RUN_LOG_FILENAME


def load_run_status(output_dir: Path) -> Optional[Dict[str, Any]]:
    path = run_status_path(output_dir)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _atomic_write_json(path: Path, payload: Dict[str, Any]) -> None:
    tmp_path = path.with_suffix(path.suffix + ".part")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(tmp_path, path)


def _artifact_payload(path: Path, checksum: bool = True) -> Dict[str, Any]:
    payload = {
        "filename": path.name,
        "bytes": path.stat().st_size,
    }
    if checksum:
        payload["sha256"] = sha256_file(path)
    return payload


def _tar_member_payload(
    archive_path: Path,
    member_name: str,
    filename: str,
    checksum: bool = True,
) -> Dict[str, Any]:
    with tarfile.open(archive_path, "r:gz") as tar_handle:
        member = tar_handle.getmember(member_name)
        payload = {
            "filename": filename,
            "archive": archive_path.name,
            "member": member_name,
            "bytes": member.size,
        }
        if checksum:
            member_file = tar_handle.extractfile(member)
            if member_file is None:
                raise RuntimeError(f"Archive member is not readable: {member_name}")
            payload["sha256"] = hashlib.sha256(member_file.read()).hexdigest()
        return payload


def _artifact_exists(output_dir: Path, artifact: Optional[Dict[str, Any]]) -> bool:
    if not artifact:
        return False
    archive_name = artifact.get("archive")
    member_name = artifact.get("member")
    if archive_name and member_name:
        return _tar_member_exists(output_dir, artifact)

    path = output_dir / artifact.get("filename", "")
    if not path.exists():
        return False
    expected_size = artifact.get("bytes")
    if expected_size is not None and path.stat().st_size != expected_size:
        return False
    expected_checksum = artifact.get("sha256")
    return expected_checksum is None or sha256_file(path) == expected_checksum


def _tar_member_exists(output_dir: Path, artifact: Dict[str, Any]) -> bool:
    archive_path = output_dir / artifact["archive"]
    if not archive_path.exists():
        return False
    try:
        with tarfile.open(archive_path, "r:gz") as tar_handle:
            member = tar_handle.getmember(artifact["member"])
            expected_size = artifact.get("bytes")
            if expected_size is not None and member.size != expected_size:
                return False
            expected_checksum = artifact.get("sha256")
            if expected_checksum is None:
                return True
            member_file = tar_handle.extractfile(member)
            if member_file is None:
                return False
            return hashlib.sha256(member_file.read()).hexdigest() == expected_checksum
    except (KeyError, OSError, tarfile.TarError, RuntimeError):
        return False


def query_completion_state_in_run_status(
    output_dir: Path, query_name: str
) -> Optional[bool]:
    """Return manifest completion state for a query.

    ``None`` means there is no usable manifest entry, so callers may consult
    legacy resume markers. ``False`` means the manifest explicitly knows this
    query and it is not safe to skip.
    """
    status = load_run_status(output_dir)
    if not status:
        return None
    queries = status.get("queries", {})
    if query_name not in queries:
        return None
    query_status = queries.get(query_name, {})
    if query_status.get("status") != "complete":
        return False
    artifacts = query_status.get("artifacts", {})
    return _artifact_exists(output_dir, artifacts.get("summary")) and _artifact_exists(
        output_dir, artifacts.get("archive")
    )


def query_is_complete_in_run_status(output_dir: Path, query_name: str) -> bool:
    return query_completion_state_in_run_status(output_dir, query_name) is True


class RunStatusRecorder:
    """Thread-safe writer for ``run_status.json`` and ``run.log``."""

    def __init__(
        self,
        output_dir: Path,
        database_path: Path,
        query_names: Iterable[str],
        settings: Dict[str, Any],
        resume_skipped_query_names: Optional[Iterable[str]] = None,
    ):
        self.output_dir = output_dir
        self.database_path = database_path
        self.query_names = list(query_names)
        self.settings = settings
        self.resume_skipped_query_names = set(resume_skipped_query_names or [])
        self.path = run_status_path(output_dir)
        self.log_path = run_log_path(output_dir)
        self.lock = threading.Lock()
        self.previous_status = load_run_status(output_dir) or {}
        self.status: Dict[str, Any] = {}

    def start(self) -> None:
        started_at = utc_now()
        previous_queries = self.previous_status.get("queries", {})
        queries = {}
        for query_name in self.query_names:
            if query_name in self.resume_skipped_query_names:
                previous = dict(previous_queries.get(query_name, {}))
                if previous.get("status") == "complete":
                    previous["status"] = "complete"
                    previous["resume_skipped"] = True
                    previous["resume_skipped_at"] = started_at
                    queries[query_name] = previous
                    continue
                queries[query_name] = self._legacy_resume_entry(query_name, started_at)
                continue
            queries[query_name] = {"status": "pending"}

        self.status = {
            "schema_version": RUN_STATUS_VERSION,
            "run_id": started_at.replace(":", "").replace("+00:00", "Z"),
            "status": "running",
            "started_at": started_at,
            "updated_at": started_at,
            "completed_at": None,
            "software_version": _SOFTWARE_VERSION,
            "database": {
                "path": str(self.database_path),
                "version": DatabaseManager.get_database_version(self.database_path),
            },
            "settings": self.settings,
            "totals": {},
            "queries": queries,
        }
        self._refresh_totals()
        with self.lock:
            log_lines = [
                f"{started_at}\tRUN_START\tsoftware={_SOFTWARE_VERSION}\t"
                f"database={self.status['database']['version']}\t"
                f"queries={len(self.query_names)}\n"
            ]
            for query_name, entry in sorted(queries.items()):
                if entry.get("resume_skipped"):
                    log_lines.append(
                        f"{started_at}\tQUERY_RESUME_SKIPPED\t{query_name}\n"
                    )
            self.log_path.write_text("".join(log_lines))
            _atomic_write_json(self.path, self.status)

    def record_query_started(self, query_name: str) -> None:
        with self.lock:
            entry = self._entry(query_name)
            if entry.get("resume_skipped"):
                return
            now = utc_now()
            entry["status"] = "running"
            entry.setdefault("started_at", now)
            entry["updated_at"] = now
            self._write_locked(f"{now}\tQUERY_START\t{query_name}\n")

    def record_query_complete(self, query_name: str) -> None:
        summary_file = self.output_dir / f"{query_name}.summary.tab"
        final_summary_file = self.output_dir / f"{query_name}.final_summary.tsv"
        archive_file = self.output_dir / f"{query_name}.tar.gz"
        missing = []
        if not archive_file.exists():
            missing.append(archive_file.name)
        summary_artifact = self._query_summary_artifact(
            query_name, summary_file, archive_file
        )
        if summary_artifact is None:
            missing.append(summary_file.name)
        if missing:
            raise RuntimeError(
                f"Cannot mark {query_name} complete; missing artifacts: {missing}"
            )

        artifacts = {
            "archive": _artifact_payload(archive_file),
            "summary": summary_artifact,
        }
        final_summary_artifact = self._query_summary_artifact(
            query_name,
            final_summary_file,
            archive_file,
            suffix="final_summary.tsv",
        )
        if final_summary_artifact is not None:
            artifacts["final_summary"] = final_summary_artifact

        with self.lock:
            entry = self._entry(query_name)
            now = utc_now()
            entry.update(
                {
                    "status": "complete",
                    "updated_at": now,
                    "completed_at": now,
                    "software_version": _SOFTWARE_VERSION,
                    "database_version": self.status["database"]["version"],
                    "artifacts": artifacts,
                }
            )
            entry.setdefault("started_at", now)
            self._write_locked(
                f"{now}\tQUERY_COMPLETE\t{query_name}\t"
                f"summary={summary_artifact['filename']}\t"
                f"archive={archive_file.name}\n"
            )

    def record_query_failed(self, query_name: str, error: str) -> None:
        with self.lock:
            entry = self._entry(query_name)
            now = utc_now()
            entry.update(
                {
                    "status": "failed",
                    "updated_at": now,
                    "completed_at": now,
                    "error": error,
                }
            )
            entry.setdefault("started_at", now)
            self._write_locked(f"{now}\tQUERY_FAILED\t{query_name}\t{error}\n")

    def finish(self, status: str) -> None:
        with self.lock:
            now = utc_now()
            self.status["status"] = status
            self.status["updated_at"] = now
            self.status["completed_at"] = now
            self._write_locked(f"{now}\tRUN_{status.upper()}\n")

    def _entry(self, query_name: str) -> Dict[str, Any]:
        return self.status.setdefault("queries", {}).setdefault(query_name, {})

    def _legacy_resume_entry(self, query_name: str, timestamp: str) -> Dict[str, Any]:
        artifacts = {}
        summary_file = self.output_dir / f"{query_name}.summary.tab"
        final_summary_file = self.output_dir / f"{query_name}.final_summary.tsv"
        archive_file = self.output_dir / f"{query_name}.tar.gz"
        if summary_file.exists():
            artifacts["summary"] = _artifact_payload(summary_file, checksum=False)
        if final_summary_file.exists():
            artifacts["final_summary"] = _artifact_payload(
                final_summary_file, checksum=False
            )
        if archive_file.exists():
            artifacts["archive"] = _artifact_payload(archive_file, checksum=False)
        return {
            "status": "complete",
            "resume_skipped": True,
            "resume_skipped_at": timestamp,
            "software_version": _SOFTWARE_VERSION,
            "database_version": DatabaseManager.get_database_version(
                self.database_path
            ),
            "artifacts": artifacts,
        }

    def _query_summary_artifact(
        self,
        query_name: str,
        summary_file: Path,
        archive_file: Path,
        suffix: str = "summary.tab",
    ) -> Optional[Dict[str, Any]]:
        if summary_file.exists():
            return _artifact_payload(summary_file)
        if not archive_file.exists():
            return None
        member_name = f"{query_name}/{query_name}.{suffix}"
        try:
            return _tar_member_payload(archive_file, member_name, summary_file.name)
        except (KeyError, OSError, tarfile.TarError, RuntimeError):
            return None

    def _write_locked(self, log_line: str) -> None:
        self.status["updated_at"] = utc_now()
        self._refresh_totals()
        with self.log_path.open("a") as log_handle:
            log_handle.write(log_line)
        _atomic_write_json(self.path, self.status)

    def _refresh_totals(self) -> None:
        counts = {
            "queries": len(self.status.get("queries", {})),
            "pending": 0,
            "running": 0,
            "complete": 0,
            "resume_skipped": 0,
            "failed": 0,
        }
        for entry in self.status.get("queries", {}).values():
            state = entry.get("status", "pending")
            if state in counts:
                counts[state] += 1
            if entry.get("resume_skipped"):
                counts["resume_skipped"] += 1
        counts["done"] = counts["complete"] + counts["failed"]
        self.status["totals"] = counts
