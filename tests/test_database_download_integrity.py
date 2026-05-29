"""M30: fail-closed database download integrity.

- An env-provided GVCLASS_DB_URL without GVCLASS_DB_SHA256 must raise at source
  resolution (before any download/fallback), unless GVCLASS_DB_ALLOW_UNVERIFIED=1.
- Non-https schemes are rejected.
- No built-in source ships without a pinned sha256.
"""

from __future__ import annotations

import pytest

from src.utils.database_manager import DatabaseManager


def test_env_url_without_sha_raises(monkeypatch):
    monkeypatch.setenv("GVCLASS_DB_URL", "https://example.com/db.tar.gz")
    monkeypatch.delenv("GVCLASS_DB_SHA256", raising=False)
    monkeypatch.delenv("GVCLASS_DB_ALLOW_UNVERIFIED", raising=False)
    with pytest.raises(RuntimeError, match="GVCLASS_DB_SHA256"):
        DatabaseManager._source_from_env()


def test_env_url_allow_unverified_opt_in(monkeypatch):
    monkeypatch.setenv("GVCLASS_DB_URL", "https://example.com/db.tar.gz")
    monkeypatch.delenv("GVCLASS_DB_SHA256", raising=False)
    monkeypatch.setenv("GVCLASS_DB_ALLOW_UNVERIFIED", "1")
    source = DatabaseManager._source_from_env()
    assert source is not None
    assert source["sha256"] is None
    assert source["url"].startswith("https://")


def test_env_url_with_sha_ok(monkeypatch):
    monkeypatch.setenv("GVCLASS_DB_URL", "https://example.com/db.tar.gz")
    monkeypatch.setenv("GVCLASS_DB_SHA256", "a" * 64)
    source = DatabaseManager._source_from_env()
    assert source["sha256"] == "a" * 64


def test_env_url_absent_returns_none(monkeypatch):
    monkeypatch.delenv("GVCLASS_DB_URL", raising=False)
    assert DatabaseManager._source_from_env() is None


@pytest.mark.parametrize("url", ["http://x/y.tar.gz", "file:///tmp/y.tar.gz", "ftp://x/y.tar.gz"])
def test_non_https_scheme_rejected(url):
    with pytest.raises(ValueError, match="https"):
        DatabaseManager._normalize_source({"url": url, "sha256": "a" * 64})


def test_no_unpinned_sources():
    for source in DatabaseManager.DATABASE_SOURCES:
        assert source.get("sha256"), f"unpinned source: {source.get('version')}"
