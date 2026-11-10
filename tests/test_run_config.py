"""The run configuration and the worker plan."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from src.pipeline.run_config import RunConfig, plan_workers
from tests.pipeline_helpers import make_config


def test_run_config_json_round_trip_keeps_every_setting(tmp_path: Path) -> None:
    config = make_config(
        tmp_path, max_workers=3, species_tree=True, min_sequence_length=0
    )

    assert RunConfig.from_json(config.to_json()) == config


def test_run_config_json_gives_paths_as_text(tmp_path: Path) -> None:
    config = make_config(tmp_path)

    assert f'"output_dir": "{config.output_dir}"' in config.to_json()


def test_run_config_is_frozen(tmp_path: Path) -> None:
    config = make_config(tmp_path)

    with pytest.raises(FrozenInstanceError):
        config.threads = 4  # type: ignore[misc]


def test_run_config_species_tree_combined_implies_species_tree(
    tmp_path: Path,
) -> None:
    config = make_config(tmp_path, species_tree=False, species_tree_combined=True)

    assert config.species_tree is True


def test_run_config_result_settings_lists_only_what_changes_results(
    tmp_path: Path,
) -> None:
    config = make_config(tmp_path, threads=64, max_workers=8, resume=True)

    assert config.result_settings() == {
        "tree_method": config.tree_method,
        "iqtree_mode": config.iqtree_mode,
        "mode_fast": config.mode_fast,
        "completeness_mode": config.completeness_mode,
        "sensitive_mode": config.sensitive_mode,
        "species_tree": config.species_tree,
        "species_tree_trim": config.species_tree_trim,
    }


def test_run_config_result_settings_ignore_thread_and_resume_settings(
    tmp_path: Path,
) -> None:
    first = make_config(tmp_path, threads=2, resume=False, contig_screen=True)
    second = make_config(tmp_path, threads=32, resume=True, contig_screen=False)

    assert first.result_settings() == second.result_settings()


@pytest.mark.parametrize(
    ("n_jobs", "threads", "max_workers", "per_worker", "expected"),
    [
        (2, 16, 8, None, (2, 8)),  # -j 8 with 2 jobs: cap before dividing
        (2, 16, None, None, (2, 8)),
        (3, 8, None, None, (3, 2)),  # a small batch: a worker for each job
        (50, 64, None, None, (16, 4)),  # a large batch: 4 threads each
        (1, 8, None, None, (1, 8)),
        (3, 4, None, None, (1, 4)),
        (200, 128, None, None, (20, 6)),  # no more than 20 workers above 100 jobs
        (50, 8, 100, None, (8, 1)),  # -j above the threads: one thread each
        (50, 10, None, 3, (3, 3)),
        (4, 8, None, 16, (1, 8)),  # a worker never gets more than the run has
        (2, 16, None, 4, (2, 4)),
        (50, 64, 8, 3, (8, 3)),
        (2, 8, 4, 4, (2, 4)),
        (50, 8, 4, 4, (4, 4)),  # both given: the user owns the total
    ],
)
def test_plan_workers_gives_workers_and_threads(
    n_jobs: int,
    threads: int,
    max_workers: int | None,
    per_worker: int | None,
    expected: tuple[int, int],
) -> None:
    assert plan_workers(n_jobs, threads, max_workers, per_worker) == expected


@pytest.mark.parametrize("n_jobs", [1, 3, 7, 50, 300])
@pytest.mark.parametrize("threads", [1, 2, 8, 16, 63, 64])
@pytest.mark.parametrize("max_workers", [None, 1, 5, 100])
@pytest.mark.parametrize("per_worker", [None, 1, 3, 100])
def test_plan_workers_never_plans_more_than_the_threads_unless_both_given(
    n_jobs: int, threads: int, max_workers: int | None, per_worker: int | None
) -> None:
    workers, per = plan_workers(n_jobs, threads, max_workers, per_worker)

    assert 1 <= workers <= n_jobs
    assert per >= 1
    if max_workers is None or per_worker is None:
        assert workers * per <= threads


@pytest.mark.parametrize(
    ("n_jobs", "threads", "max_workers", "per_worker"),
    [(0, 8, None, None), (2, 0, None, None), (2, 8, 0, None), (2, 8, None, -1)],
)
def test_plan_workers_rejects_values_below_one(
    n_jobs: int, threads: int, max_workers: int | None, per_worker: int | None
) -> None:
    with pytest.raises(ValueError, match="at least 1"):
        plan_workers(n_jobs, threads, max_workers, per_worker)
