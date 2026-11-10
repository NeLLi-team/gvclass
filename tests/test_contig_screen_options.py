"""The options that switch the contig screen: CLI flag, config key, runner flags."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
import yaml
from click.testing import CliRunner

from src.bin import gvclass_runner, setup_database
from src.bin.gvclass_cli import (
    BUILTIN_CONFIG_DEFAULTS,
    ContigInput,
    PipelineContext,
    WorkerPlan,
    build_parser,
    build_pipeline_command,
    describe_contig_screen,
    resolve_contig_screen,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def make_context(*, contig_screen: bool = True, contigs: bool = False) -> Any:
    return PipelineContext(
        query_dir=Path("query"),
        active_query_dir=Path("query"),
        output_dir=Path("out"),
        database=Path("resources"),
        threads=8,
        tree_method="veryfasttree",
        iqtree_mode="fast",
        mode_fast=True,
        completeness_mode="novelty-aware",
        sensitive_mode=True,
        contigs=ContigInput(enabled=contigs),
        input_min_length=20000,
        contigs_min_length=10000,
        validation_min_length=20000,
        temp_contigs_dir=None,
        n_queries=2,
        n_input_files=1,
        n_skipped=0,
        workers=WorkerPlan(workers=2, threads_per_worker=4),
        contig_screen=contig_screen,
    )


def command_for(context: PipelineContext) -> list[str]:
    args = SimpleNamespace(
        cluster_type="local",
        cluster_queue=None,
        cluster_project=None,
        cluster_walltime=None,
        verbose=False,
        resume=False,
    )
    return build_pipeline_command(
        args=args,
        context=context,
        python_exe="python",
        pixi_cmd=None,
        use_pixi_run=False,
    )


def test_the_parser_takes_no_screen() -> None:
    parser = build_parser()

    assert parser.parse_args(["bins", "--no-screen"]).no_screen is True
    assert parser.parse_args(["bins"]).no_screen is False


def test_the_screen_is_on_by_default_and_the_flag_turns_it_off() -> None:
    config = {"pipeline": {"contig_screen": True}}

    assert resolve_contig_screen(SimpleNamespace(no_screen=False), config) is True
    assert resolve_contig_screen(SimpleNamespace(no_screen=True), config) is False


def test_the_config_key_turns_the_screen_off() -> None:
    config = {"pipeline": {"contig_screen": False}}

    assert resolve_contig_screen(SimpleNamespace(no_screen=False), config) is False


def test_a_config_without_the_key_keeps_the_screen_on() -> None:
    assert (
        resolve_contig_screen(SimpleNamespace(no_screen=False), {"pipeline": {}})
        is True
    )


def test_the_builtin_defaults_and_both_yaml_files_enable_the_screen() -> None:
    shipped = yaml.safe_load((REPO_ROOT / "config" / "gvclass_config.yaml").read_text())
    test_config = yaml.safe_load(
        (REPO_ROOT / "config" / "gvclass_config_test.yaml").read_text()
    )

    assert BUILTIN_CONFIG_DEFAULTS["pipeline"]["contig_screen"] is True
    assert shipped["pipeline"]["contig_screen"] is True
    assert test_config["pipeline"]["contig_screen"] is True
    assert shipped == BUILTIN_CONFIG_DEFAULTS


def test_the_runner_command_omits_the_screen_flags_by_default() -> None:
    command = command_for(make_context())

    assert "--no-screen" not in command
    assert "--contigs" not in command


def test_the_runner_command_forwards_no_screen() -> None:
    command = command_for(make_context(contig_screen=False))

    assert "--no-screen" in command


def test_the_runner_command_tells_the_runner_about_contigs_mode() -> None:
    command = command_for(make_context(contigs=True))

    assert "--contigs" in command


def run_runner(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, *flags: str
) -> dict[str, Any]:
    seen: dict[str, Any] = {}
    monkeypatch.setattr(
        gvclass_runner,
        "gvclass_flow",
        lambda config: (
            seen.update(
                contig_screen=config.contig_screen, contigs_mode=config.contigs_mode
            )
            or "s.tsv"
        ),
    )
    queries = tmp_path / "queries"
    queries.mkdir(exist_ok=True)
    outcome = CliRunner().invoke(
        gvclass_runner.main, [str(queries), "-o", str(tmp_path / "out"), *flags]
    )
    assert outcome.exit_code == 0, outcome.output
    return seen


def test_the_runner_screens_by_default(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    seen = run_runner(monkeypatch, tmp_path)

    assert seen["contig_screen"] is True
    assert seen["contigs_mode"] is False


def test_the_runner_takes_no_screen_and_contigs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    seen = run_runner(monkeypatch, tmp_path, "--no-screen", "--contigs")

    assert seen["contig_screen"] is False
    assert seen["contigs_mode"] is True


def test_setup_db_installs_the_addon_only_when_asked(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    installed: list[Path] = []
    monkeypatch.setattr(
        setup_database, "check_and_setup_database", lambda *a, **k: True
    )
    monkeypatch.setattr(
        "src.bin.gvclass_screen.install_bundle",
        lambda database: installed.append(database) or 0,
    )
    monkeypatch.setattr(sys, "argv", ["setup_database", "-d", str(tmp_path / "db")])
    assert setup_database.main() == 0
    assert installed == []

    monkeypatch.setattr(
        sys, "argv", ["setup_database", "-d", str(tmp_path / "db"), "--with-screen"]
    )

    assert setup_database.main() == 0
    assert installed == [(tmp_path / "db").resolve()]


def test_setup_db_reports_a_failed_addon_install(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(
        setup_database, "check_and_setup_database", lambda *a, **k: True
    )
    monkeypatch.setattr("src.bin.gvclass_screen.install_bundle", lambda database: 1)
    monkeypatch.setattr(
        sys, "argv", ["setup_database", "-d", str(tmp_path / "db"), "--with-screen"]
    )

    assert setup_database.main() == 1


def test_setup_db_never_asks_a_question() -> None:
    source = (REPO_ROOT / "src" / "bin" / "setup_database.py").read_text()

    assert "input(" not in source
    assert "prompt_yes_no" not in source


def test_the_startup_lines_say_what_the_screen_does(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    context = make_context()
    context.database = tmp_path
    assert "gvclass screen --install" in describe_contig_screen(context)

    (tmp_path / "contig_screen").mkdir()
    (tmp_path / "contig_screen" / "manifest.json").write_text("{}")
    assert describe_contig_screen(context) == "enabled"

    context.contig_screen = False
    assert describe_contig_screen(context) == "disabled"

    context.contig_screen = True
    context.contigs = ContigInput(enabled=True)
    assert "--contigs" in describe_contig_screen(context)


def test_the_startup_lines_survive_an_unreadable_manifest(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def refuse(database: Path) -> bool:
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr("src.contig_screen.layout.bundle_present", refuse)
    context = make_context()
    context.database = tmp_path

    assert "not readable" in describe_contig_screen(context)

    context.contig_screen = False
    assert describe_contig_screen(context) == "disabled"
