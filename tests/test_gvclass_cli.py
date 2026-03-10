from src.bin.gvclass_cli import (
    inspect_local_database_state,
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
