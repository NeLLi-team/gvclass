"""Security-focused tests for database extraction helpers."""

from pathlib import Path

import pytest

from src.utils.database_manager import DatabaseManager


def test_normalize_member_path_rejects_parent_components() -> None:
    with pytest.raises(RuntimeError, match="path traversal"):
        DatabaseManager._normalize_member_path("../escape/file", strip_components=0)


def test_resolve_safe_target_rejects_escape_path(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="path traversal"):
        DatabaseManager._resolve_safe_target(
            base_path=tmp_path,
            relative_path=Path("../outside.txt"),
            member_name="../outside.txt",
        )


def test_normalize_member_path_strips_top_level_directory() -> None:
    normalized = DatabaseManager._normalize_member_path(
        "resources_v1_2_1/models/combined.hmm", strip_components=1
    )
    assert normalized == Path("models/combined.hmm")
