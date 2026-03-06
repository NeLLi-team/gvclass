"""Tests for sensitive-mode CLI plumbing."""

from pathlib import Path
from types import SimpleNamespace

from src.bin.gvclass_cli import (
    ContigInput,
    PipelineContext,
    WorkerPlan,
    build_pipeline_command,
    resolve_sensitive_mode,
)


def _make_args(**overrides):
    base = {
        "cluster_type": "local",
        "cluster_queue": None,
        "cluster_project": None,
        "cluster_walltime": "04:00:00",
        "verbose": False,
        "resume": False,
        "sensitive": False,
        "extended": False,
        "mode_fast": False,
    }
    base.update(overrides)
    return SimpleNamespace(**base)


def _make_context(sensitive_mode: bool) -> PipelineContext:
    return PipelineContext(
        query_dir=Path("query"),
        active_query_dir=Path("query"),
        output_dir=Path("out"),
        database=Path("resources"),
        threads=24,
        tree_method="fasttree",
        mode_fast=True,
        completeness_mode="legacy",
        sensitive_mode=sensitive_mode,
        contigs=ContigInput(enabled=False),
        contigs_min_length=10000,
        temp_contigs_dir=None,
        n_queries=2,
        n_input_files=1,
        n_skipped=0,
        workers=WorkerPlan(workers=2, threads_per_worker=12),
    )


def test_resolve_sensitive_mode_from_args() -> None:
    args = _make_args(sensitive=True)
    config = {"pipeline": {"sensitive_mode": False}}
    assert resolve_sensitive_mode(args, config) is True


def test_resolve_sensitive_mode_from_config() -> None:
    args = _make_args(sensitive=False)
    config = {"pipeline": {"sensitive_mode": True}}
    assert resolve_sensitive_mode(args, config) is True


def test_build_pipeline_command_includes_sensitive_flag() -> None:
    args = _make_args()
    context = _make_context(sensitive_mode=True)
    cmd = build_pipeline_command(
        args=args,
        context=context,
        python_exe="python",
        pixi_cmd=None,
        use_pixi_run=False,
    )
    assert "--sensitive" in cmd


def test_build_pipeline_command_omits_sensitive_flag() -> None:
    args = _make_args()
    context = _make_context(sensitive_mode=False)
    cmd = build_pipeline_command(
        args=args,
        context=context,
        python_exe="python",
        pixi_cmd=None,
        use_pixi_run=False,
    )
    assert "--sensitive" not in cmd
