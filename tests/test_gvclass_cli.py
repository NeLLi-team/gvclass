from pathlib import Path

from src.bin.gvclass_cli import (
    CliOutput,
    inspect_local_database_state,
    load_config,
    prompt_yes_no,
    should_offer_database_update,
)
from src.utils.database_manager import DatabaseManager


def test_inspect_local_database_state_requires_database_directory(tmp_path):
    db_path = tmp_path / "resources"
    db_path.mkdir()
    for relative_path in DatabaseManager.REQUIRED_FILES:
        target = db_path / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("test\n")

    state = inspect_local_database_state(db_path)

    assert state["exists"] is True
    assert state["complete"] is False
    assert "database" in state["missing_files"]


def test_prompt_yes_no_accepts_enter_as_default_yes(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda _: "")

    assert prompt_yes_no("Update database?", default=True) is True


def test_should_offer_database_update_uses_version_comparison():
    assert should_offer_database_update("v1.2.2", {"version": "v1.4.0"}) is True
    assert should_offer_database_update("v1.4.0", {"version": "v1.4.0"}) is False


def _load(tmp_path: Path, text: str):
    cfg = tmp_path / "gvclass_config.yaml"
    cfg.write_text(text)
    return load_config(str(cfg), tmp_path, CliOutput(plain_output=True))


def test_load_config_partial_pipeline_uses_defaults(tmp_path):
    # Only one pipeline key set; the rest must fall back to defaults.
    config = _load(tmp_path, "pipeline:\n  threads: 4\n")
    assert config["pipeline"]["threads"] == 4
    assert config["pipeline"]["tree_method"] == "iqtree"
    assert config["pipeline"]["mode_fast"] is True
    assert config["pipeline"]["completeness_mode"] == "novelty-aware"
    # The whole database section must still be present.
    assert config["database"]["download_version"] == "v1.6.0"


def test_load_config_missing_database_uses_defaults(tmp_path):
    config = _load(tmp_path, "pipeline:\n  tree_method: iqtree\n")
    assert config["pipeline"]["tree_method"] == "iqtree"
    assert "database" in config
    assert config["database"]["download_url"].startswith("https://")


def test_load_config_partial_database_preserves_download_fields(tmp_path):
    config = _load(tmp_path, "database:\n  path: /custom/db\n")
    assert config["database"]["path"] == "/custom/db"
    # Setting only path must not drop the download fields needed for setup-db.
    assert config["database"]["download_url"].startswith("https://")
    assert config["database"]["download_sha256"]
    assert config["pipeline"]["sensitive_mode"] is True


def test_load_config_empty_file_returns_defaults(tmp_path):
    config = _load(tmp_path, "# only a comment\n")
    assert config["pipeline"]["threads"] == 16
    assert config["database"]["download_version"] == "v1.6.0"
