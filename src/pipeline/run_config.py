"""The settings of one run, and the plan that turns threads into workers.

``RunConfig`` holds every resolved setting. The command line builds it once; the
runner and everything below it read it. It has no defaults: the caller names each
setting, so a default lives in one place only.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

#: Without ``-j`` or ``--threads-per-worker``, aim for this many threads per worker.
TARGET_THREADS_PER_WORKER = 4
#: A batch of more than ``LARGE_BATCH`` queries gets at most this many workers.
LARGE_BATCH = 100
MAX_LARGE_BATCH_WORKERS = 20

_PATH_FIELDS = ("query_path", "output_dir", "database_path")


@dataclass(frozen=True, slots=True)
class RunConfig:
    """Every resolved setting of a run."""

    query_path: Path
    output_dir: Path
    database_path: Path
    threads: int
    max_workers: int | None
    threads_per_worker: int | None
    tree_method: str
    iqtree_mode: str
    mode_fast: bool
    completeness_mode: str
    sensitive_mode: bool
    resume: bool
    allow_short: bool
    min_sequence_length: int | None
    species_tree: bool
    species_tree_combined: bool
    species_tree_trim: str
    contig_screen: bool
    contigs_mode: bool
    contigs_min_length: int

    def __post_init__(self) -> None:
        """Make ``--species-tree-combined`` imply ``--species-tree``."""
        if self.species_tree_combined and not self.species_tree:
            object.__setattr__(self, "species_tree", True)

    def result_settings(self) -> dict[str, Any]:
        """Return the settings that change the scientific result of a query.

        ``--resume`` reuses a saved query only when these are equal.
        """
        return {
            "tree_method": self.tree_method,
            "iqtree_mode": self.iqtree_mode,
            "mode_fast": self.mode_fast,
            "completeness_mode": self.completeness_mode,
            "sensitive_mode": self.sensitive_mode,
            "species_tree": self.species_tree,
            "species_tree_trim": self.species_tree_trim,
        }

    def to_dict(self) -> dict[str, Any]:
        """Return all settings, with paths as text, ready for JSON."""
        values = asdict(self)
        for name in _PATH_FIELDS:
            values[name] = str(values[name])
        return values

    def to_json(self) -> str:
        """Return the configuration as JSON text."""
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_json(cls, text: str) -> RunConfig:
        """Build a configuration from the text that ``to_json`` wrote."""
        values = json.loads(text)
        for name in _PATH_FIELDS:
            values[name] = Path(values[name])
        return cls(**values)


def plan_workers(
    n_jobs: int,
    total_threads: int,
    max_workers: int | None = None,
    threads_per_worker: int | None = None,
) -> tuple[int, int]:
    """Return ``(workers, threads_per_worker)`` for ``n_jobs`` jobs.

    The worker count is capped at ``n_jobs`` before the threads are divided, so
    ``-j 8`` with two jobs and 16 threads gives two workers with eight threads each.
    The plan uses no more than ``total_threads`` unless the caller gave both
    ``max_workers`` and ``threads_per_worker``; then the caller owns the total.
    Without either, a batch in which every job can have a worker with at least two
    threads gets one worker per job; a larger batch gets workers of about four threads.

    Raises:
        ValueError: If an argument is below 1.
    """
    for name, value in (
        ("n_jobs", n_jobs),
        ("total_threads", total_threads),
        ("max_workers", max_workers),
        ("threads_per_worker", threads_per_worker),
    ):
        if value is not None and value < 1:
            raise ValueError(f"{name} must be at least 1, got {value}")

    if max_workers is not None and threads_per_worker is not None:
        return min(max_workers, n_jobs), threads_per_worker
    if threads_per_worker is not None:
        per_worker = min(threads_per_worker, total_threads)
        return min(n_jobs, total_threads // per_worker), per_worker
    if max_workers is not None:
        workers = min(max_workers, n_jobs, total_threads)
        return workers, total_threads // workers
    return _automatic_plan(n_jobs, total_threads)


def _automatic_plan(n_jobs: int, total_threads: int) -> tuple[int, int]:
    """Plan workers when the caller set neither ``-j`` nor ``--threads-per-worker``."""
    if total_threads >= 2 * n_jobs:
        return n_jobs, total_threads // n_jobs
    workers = max(1, min(n_jobs, total_threads // TARGET_THREADS_PER_WORKER))
    if n_jobs > LARGE_BATCH:
        workers = min(workers, MAX_LARGE_BATCH_WORKERS)
    return workers, total_threads // workers
