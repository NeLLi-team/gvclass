"""Builders and fakes for the tests of the run flow.

``fake_flow`` replaces the three steps of the flow that need a database or a
long compute: the input check, the marker-panel check and the per-query work. The fake
query step does what the real one promises: it publishes ``<query>.tar.gz`` with
the two summary members and returns the result of the query.
"""

from __future__ import annotations

import io
import tarfile
import tempfile
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pytest

from src.pipeline import parallel_runner
from src.pipeline.parallel_runner import RunInputs
from src.pipeline.run_config import RunConfig
from src.pipeline.summary_writer import write_individual_final_summary_file

CURRENT_QUALITY_VALUES = {
    "completeness_query_support": "supported",
    "estimated_completeness_strategy": "novelty_aware_v1",
    "completeness_model_group": "Algavirales",
    "completeness_model_reliability": "unvalidated_legacy",
    "contamination_query_support": "supported",
    "estimated_contamination_strategy": "extra_trees_v1",
    "contamination_model_reliability": "unvalidated_legacy",
}


def make_config(tmp_path: Path, **overrides: Any) -> RunConfig:
    """Return a configuration for one small run under ``tmp_path``."""
    values: dict[str, Any] = {
        "query_path": tmp_path / "queries",
        "output_dir": tmp_path / "out",
        "database_path": tmp_path / "db",
        "threads": 2,
        "max_workers": 1,
        "threads_per_worker": 1,
        "tree_method": "veryfasttree",
        "iqtree_mode": "fast",
        "mode_fast": True,
        "completeness_mode": "novelty-aware",
        "sensitive_mode": True,
        "resume": False,
        "allow_short": True,
        "min_sequence_length": None,
        "species_tree": False,
        "species_tree_combined": False,
        "species_tree_trim": "witchi",
        "contig_screen": True,
        "contigs_mode": False,
        "contigs_min_length": 0,
    }
    return RunConfig(**{**values, **overrides})


def write_queries(
    config: RunConfig, names: list[str], suffix: str = ".faa", text: str | None = None
) -> list[Path]:
    """Write one small query file per name; return the paths."""
    config.query_path.mkdir(parents=True, exist_ok=True)
    body = text or (">c\nACGT\n" if suffix == ".fna" else ">p\nM\n")
    paths = [config.query_path / f"{name}{suffix}" for name in names]
    for path in paths:
        path.write_text(body, encoding="utf-8")
    return paths


def fake_setup(config: RunConfig) -> RunInputs:
    """Stand in for the input check: every ``.faa`` and ``.fna`` is valid."""
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.database_path.mkdir(parents=True, exist_ok=True)
    files = sorted([*config.query_path.glob("*.faa"), *config.query_path.glob("*.fna")])
    return RunInputs(tuple(files), (), config.database_path)


def write_query_archive(
    output_dir: Path, name: str, summary: Mapping[str, Any]
) -> Path:
    """Publish ``<name>.tar.gz`` with the two per-query summary files."""
    archive = output_dir / f"{name}.tar.gz"
    with tempfile.TemporaryDirectory() as scratch:
        work = Path(scratch) / name
        work.mkdir()
        (work / f"{name}.summary.tab").write_text(f"query\n{name}\n", encoding="utf-8")
        write_individual_final_summary_file(work, name, dict(summary))
        with tarfile.open(archive, "w:gz") as handle:
            handle.add(work, arcname=name)
    return archive


def make_fake_process(
    calls: list[str] | None = None,
    summaries: Mapping[str, Mapping[str, Any]] | None = None,
    failing: Mapping[str, BaseException] | None = None,
    threads_seen: list[int] | None = None,
) -> Callable[[Path, RunConfig, int], dict[str, Any]]:
    """Return a stand-in for ``process_single_query`` that logs its queries."""

    def process(query_file: Path, config: RunConfig, threads: int) -> dict[str, Any]:
        name = query_file.stem
        if calls is not None:
            calls.append(name)
        if threads_seen is not None:
            threads_seen.append(threads)
        if failing and name in failing:
            raise failing[name]
        summary = {
            "query": name,
            "taxonomy_majority": f"processed_{name}",
            "estimated_completeness": 91.25,
            **CURRENT_QUALITY_VALUES,
            **(summaries or {}).get(name, {}),
        }
        write_query_archive(config.output_dir, name, summary)
        return {
            "query": name,
            "status": "complete",
            "summary_data": summary,
            "genetic_code": "N/A",
            "output_dir": str(config.output_dir),
        }

    return process


def fake_flow(
    monkeypatch: pytest.MonkeyPatch,
    calls: list[str] | None = None,
    summaries: Mapping[str, Mapping[str, Any]] | None = None,
    failing: Mapping[str, BaseException] | None = None,
    threads_seen: list[int] | None = None,
) -> None:
    """Replace the input check, the panel check and the query step of the flow."""
    monkeypatch.setattr(parallel_runner, "validate_and_setup_task", fake_setup)
    monkeypatch.setattr(parallel_runner, "check_marker_panels", lambda _database: None)
    monkeypatch.setattr(
        parallel_runner,
        "process_single_query",
        make_fake_process(calls, summaries, failing, threads_seen),
    )


def read_jsonl_lines(path: Path) -> list[str]:
    """Return the non-empty lines of a text file."""
    return [line for line in path.read_text(encoding="utf-8").splitlines() if line]


def tar_member_text(archive: Path, member: str) -> str:
    """Return the text of one archive member."""
    with tarfile.open(archive, "r:gz") as handle:
        extracted = handle.extractfile(member)
        assert extracted is not None
        return io.TextIOWrapper(extracted, encoding="utf-8").read()
