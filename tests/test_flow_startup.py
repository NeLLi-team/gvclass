"""The flow checks the database and plans the workers before it processes a query."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.config.marker_sets import PANEL_MODELS, PANEL_PREFIXES
from src.pipeline import parallel_runner
from src.utils.error_handling import ProcessingError
from tests.pipeline_helpers import (
    fake_flow,
    fake_setup,
    make_config,
    make_fake_process,
    write_queries,
)


def panel_models() -> set[str]:
    """Return one model name for every panel entry and every prefix panel."""
    named = {model for models in PANEL_MODELS.values() for model in models}
    return named | {f"{prefix}1" for prefix in PANEL_PREFIXES.values()}


def write_hmm(database: Path, names: set[str]) -> None:
    """Write a ``hmm/combined.hmm`` that holds only the ``NAME`` lines of ``names``."""
    (database / "hmm").mkdir(parents=True, exist_ok=True)
    text = "".join(f"HMMER3/f\nNAME  {name}\nLENG  10\n//\n" for name in sorted(names))
    (database / "hmm" / "combined.hmm").write_text(text, encoding="utf-8")


def test_the_query_pool_plans_workers_from_the_number_of_queries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    threads_seen: list[int] = []
    fake_flow(monkeypatch, threads_seen=threads_seen)
    config = make_config(
        tmp_path,
        threads=16,
        max_workers=8,
        threads_per_worker=None,
        contig_screen=False,
    )
    write_queries(config, ["q1", "q2"])

    parallel_runner.gvclass_flow(config)

    assert threads_seen == [8, 8]


def test_the_panel_check_accepts_a_database_with_every_panel_model(
    tmp_path: Path,
) -> None:
    write_hmm(tmp_path, panel_models() | {"OG1091"})

    parallel_runner.check_marker_panels(tmp_path)


def test_the_panel_check_names_the_panel_and_the_missing_model(
    tmp_path: Path,
) -> None:
    missing = PANEL_MODELS["phage"][0]
    write_hmm(tmp_path, panel_models() - {missing})

    with pytest.raises(ProcessingError, match=rf"phage: {missing}"):
        parallel_runner.check_marker_panels(tmp_path)


def test_the_panel_check_names_the_missing_capsid_model(
    tmp_path: Path,
) -> None:
    missing = PANEL_MODELS["capsid"][0]
    write_hmm(tmp_path, panel_models() - {missing})

    with pytest.raises(ProcessingError, match=rf"capsid: {missing}"):
        parallel_runner.check_marker_panels(tmp_path)


def test_the_panel_check_names_a_prefix_panel_that_matches_no_model(
    tmp_path: Path,
) -> None:
    prefix = PANEL_PREFIXES["plv"]
    write_hmm(
        tmp_path, {name for name in panel_models() if not name.startswith(prefix)}
    )

    with pytest.raises(ProcessingError, match=rf"plv: {prefix}\*"):
        parallel_runner.check_marker_panels(tmp_path)


def test_the_flow_stops_before_any_query_when_a_panel_model_is_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[str] = []
    monkeypatch.setattr(parallel_runner, "validate_and_setup_task", fake_setup)
    monkeypatch.setattr(
        parallel_runner, "process_single_query", make_fake_process(calls)
    )
    config = make_config(tmp_path)
    write_queries(config, ["q1"])
    write_hmm(config.database_path, panel_models() - {PANEL_MODELS["mrya"][0]})

    with pytest.raises(ProcessingError, match="mrya"):
        parallel_runner.gvclass_flow(config)

    assert calls == []
    assert not (config.output_dir / "results.jsonl").exists()
