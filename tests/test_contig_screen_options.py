"""The options that switch the contig screen: flag, config key, run configuration."""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
import yaml

from src.bin import setup_database
from src.bin.gvclass_cli import (
    BUILTIN_CONFIG_DEFAULTS,
    build_parser,
    describe_contig_screen,
    resolve_contig_screen,
)
from tests.cli_helpers import run_cli
from tests.pipeline_helpers import make_config

REPO_ROOT = Path(__file__).resolve().parents[1]


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


def test_the_builtin_defaults_and_both_yaml_files_enable_the_screen() -> None:
    shipped = yaml.safe_load((REPO_ROOT / "config" / "gvclass_config.yaml").read_text())
    test_config = yaml.safe_load(
        (REPO_ROOT / "config" / "gvclass_config_test.yaml").read_text()
    )

    assert BUILTIN_CONFIG_DEFAULTS["pipeline"]["contig_screen"] is True
    assert shipped["pipeline"]["contig_screen"] is True
    assert test_config["pipeline"]["contig_screen"] is True
    assert shipped == BUILTIN_CONFIG_DEFAULTS


def test_the_run_config_screens_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(tmp_path, monkeypatch)

    assert run.config.contig_screen is True
    assert run.config.contigs_mode is False


def test_the_run_config_carries_no_screen(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(tmp_path, monkeypatch, "--no-screen")

    assert run.config.contig_screen is False


def test_the_run_config_tells_the_runner_about_contigs_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(
        tmp_path,
        monkeypatch,
        "--contigs",
        suffix=".fna",
        body=">c1\n" + "ACGT" * 3000 + "\n",
    )

    assert run.config.contigs_mode is True


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


@pytest.mark.parametrize(("flags", "update"), [([], False), (["--update"], True)])
def test_setup_db_passes_the_update_consent_to_the_installer(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    flags: list[str],
    update: bool,
) -> None:
    check = Mock(return_value=True)
    monkeypatch.setattr(setup_database, "check_and_setup_database", check)
    monkeypatch.setattr(
        sys, "argv", ["setup_database", "-d", str(tmp_path / "db"), *flags]
    )

    assert setup_database.main() == 0
    assert check.call_args.kwargs["update"] is update


def test_the_startup_lines_say_what_the_screen_does(tmp_path: Path) -> None:
    config = make_config(tmp_path, database_path=tmp_path, contig_screen=True)
    assert "gvclass screen --install" in describe_contig_screen(config)

    (tmp_path / "contig_screen").mkdir()
    (tmp_path / "contig_screen" / "manifest.json").write_text("{}")
    assert describe_contig_screen(config) == "enabled"

    assert describe_contig_screen(replace(config, contig_screen=False)) == "disabled"
    assert "--contigs" in describe_contig_screen(replace(config, contigs_mode=True))


def test_the_startup_lines_survive_an_unreadable_manifest(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def refuse(database: Path) -> bool:
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr("src.contig_screen.layout.bundle_present", refuse)
    config = make_config(tmp_path, database_path=tmp_path, contig_screen=True)

    assert "not readable" in describe_contig_screen(config)
    assert describe_contig_screen(replace(config, contig_screen=False)) == "disabled"
