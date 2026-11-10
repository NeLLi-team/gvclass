"""A complete database is replaced only with consent: a prompt or ``--update``."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pytest

from src.bin import gvclass_cli
from src.utils.database_manager import DatabaseManager

PINNED = {"version": "v2.1.0", "url": "https://example.test/db.tar.gz"}


def make_database(path: Path, version: str | None) -> Path:
    for relative in DatabaseManager.REQUIRED_FILES:
        target = path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("x\n")
    if version is not None:
        (path / "DB_VERSION").write_text(f"{version}\n")
    return path


@pytest.fixture
def installer(monkeypatch: pytest.MonkeyPatch) -> Mock:
    install = Mock()
    monkeypatch.setattr(DatabaseManager, "setup_database", install)
    monkeypatch.setattr(
        DatabaseManager, "get_effective_database_source", Mock(return_value=PINNED)
    )
    monkeypatch.setattr(
        DatabaseManager, "get_latest_database_source", Mock(return_value=None)
    )
    return install


def set_terminal(monkeypatch: pytest.MonkeyPatch, *, interactive: bool) -> None:
    monkeypatch.setattr(gvclass_cli, "is_interactive_terminal", lambda: interactive)


def run_check(db: Path, *, update: bool = False) -> bool:
    output = gvclass_cli.CliOutput(plain_output=True)
    return gvclass_cli.check_and_setup_database(db, {}, output, update=update)


def test_a_run_without_a_terminal_keeps_an_older_database_and_warns_once(
    tmp_path, monkeypatch, installer, capsys
):
    db = make_database(tmp_path / "resources", "v2.0.0")
    set_terminal(monkeypatch, interactive=False)

    assert run_check(db) is True

    installer.assert_not_called()
    warnings = [
        line for line in capsys.readouterr().err.splitlines() if "v2.0.0" in line
    ]
    assert len(warnings) == 1
    assert "v2.1.0" in warnings[0]
    assert "setup-db --update" in warnings[0]


def test_a_terminal_run_asks_and_the_default_answer_keeps_the_database(
    tmp_path, monkeypatch, installer
):
    db = make_database(tmp_path / "resources", "v2.0.0")
    set_terminal(monkeypatch, interactive=True)
    monkeypatch.setattr("builtins.input", lambda prompt: "")

    assert run_check(db) is True

    installer.assert_not_called()


def test_a_terminal_run_updates_after_a_yes(tmp_path, monkeypatch, installer):
    db = make_database(tmp_path / "resources", "v2.0.0")
    set_terminal(monkeypatch, interactive=True)
    monkeypatch.setattr("builtins.input", lambda prompt: "y")

    assert run_check(db) is True

    assert installer.call_args.kwargs["force"] is True


@pytest.mark.parametrize("version", [None, "dev-build"])
def test_a_database_without_a_release_version_is_never_updated(
    tmp_path, monkeypatch, installer, capsys, version
):
    db = make_database(tmp_path / "resources", version)
    set_terminal(monkeypatch, interactive=True)
    monkeypatch.setattr("builtins.input", lambda prompt: pytest.fail("asked"))

    assert run_check(db) is True

    installer.assert_not_called()
    custom = [line for line in capsys.readouterr().out.splitlines() if "custom" in line]
    assert len(custom) == 1


@pytest.mark.parametrize("version", ["v2.1.0", "v2.2.0"])
def test_a_current_database_is_left_alone(tmp_path, monkeypatch, installer, version):
    db = make_database(tmp_path / "resources", version)
    set_terminal(monkeypatch, interactive=False)

    assert run_check(db) is True

    installer.assert_not_called()


@pytest.mark.parametrize("version", [None, "v2.0.0", "v2.1.0"])
def test_the_explicit_update_replaces_a_complete_database(
    tmp_path, monkeypatch, installer, version
):
    db = make_database(tmp_path / "resources", version)
    set_terminal(monkeypatch, interactive=False)

    assert run_check(db, update=True) is True

    assert installer.call_args.kwargs == {"preferred_source": PINNED, "force": True}


def test_a_missing_database_is_installed(tmp_path, monkeypatch, installer):
    set_terminal(monkeypatch, interactive=False)

    assert run_check(tmp_path / "resources") is True

    assert installer.call_args.kwargs == {"preferred_source": PINNED, "force": False}


def test_an_incomplete_database_is_installed(tmp_path, monkeypatch, installer):
    db = make_database(tmp_path / "resources", "v2.0.0")
    (db / "labels.tsv").unlink()
    set_terminal(monkeypatch, interactive=False)

    assert run_check(db) is True

    assert installer.call_args.kwargs == {"preferred_source": PINNED, "force": False}


def test_a_failed_update_reports_failure(tmp_path, monkeypatch, installer):
    db = make_database(tmp_path / "resources", "v2.0.0")
    set_terminal(monkeypatch, interactive=False)
    installer.side_effect = RuntimeError("down")

    assert run_check(db, update=True) is False


@pytest.mark.parametrize("version", ["unknown", "not installed", "dev-build"])
def test_no_update_is_offered_for_a_version_that_cannot_be_compared(version):
    assert gvclass_cli.should_offer_database_update(version, PINNED) is False
