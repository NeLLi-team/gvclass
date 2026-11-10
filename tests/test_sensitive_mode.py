"""Tests for sensitive-mode CLI plumbing."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from src.bin.gvclass_cli import (
    BUILTIN_CONFIG_DEFAULTS,
    build_parser,
    resolve_sensitive_mode,
)
from tests.cli_helpers import run_cli


def _make_args(**overrides):
    base = {"sensitive": False, "no_sensitive": False}
    base.update(overrides)
    return SimpleNamespace(**base)


def test_resolve_sensitive_mode_from_args() -> None:
    args = _make_args(sensitive=True)
    config = {"pipeline": {"sensitive_mode": False}}
    assert resolve_sensitive_mode(args, config) is True


def test_resolve_sensitive_mode_from_config() -> None:
    args = _make_args(sensitive=False)
    config = {"pipeline": {"sensitive_mode": True}}
    assert resolve_sensitive_mode(args, config) is True


def test_the_run_config_is_sensitive_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert run_cli(tmp_path, monkeypatch).config.sensitive_mode is True


def test_the_run_config_turns_sensitive_mode_off_with_no_sensitive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = run_cli(tmp_path, monkeypatch, "--no-sensitive")

    assert run.config.sensitive_mode is False


def test_cli_no_sensitive_flag() -> None:
    """--no-sensitive should resolve sensitive_mode to False."""
    parser = build_parser()
    default_args = parser.parse_args(["example", "-o", "out"])
    assert resolve_sensitive_mode(default_args, BUILTIN_CONFIG_DEFAULTS) is True

    no_sens_args = parser.parse_args(["example", "-o", "out", "--no-sensitive"])
    assert resolve_sensitive_mode(no_sens_args, BUILTIN_CONFIG_DEFAULTS) is False
