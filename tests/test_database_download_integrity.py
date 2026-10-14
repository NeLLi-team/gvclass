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
    assert source["allow_unverified"] is True


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


def test_source_resolution_uses_only_environment_override(monkeypatch):
    monkeypatch.setenv("GVCLASS_DB_URL", "https://example.com/env.tar.gz")
    monkeypatch.setenv("GVCLASS_DB_SHA256", "a" * 64)
    preferred = {
        "version": "v8.0.0",
        "url": "https://example.com/config.tar.gz",
        "sha256": "b" * 64,
    }

    sources = DatabaseManager._resolve_database_sources(preferred)

    assert [source["url"] for source in sources] == ["https://example.com/env.tar.gz"]


def test_source_resolution_uses_only_preferred_source(monkeypatch):
    monkeypatch.delenv("GVCLASS_DB_URL", raising=False)
    preferred = {
        "version": "v8.0.0",
        "url": "https://example.com/config.tar.gz",
        "sha256": "b" * 64,
    }

    sources = DatabaseManager._resolve_database_sources(preferred)

    assert sources == [
        {
            "version": "v8.0.0",
            "url": "https://example.com/config.tar.gz",
            "filename": "config.tar.gz",
            "sha256": "b" * 64,
        }
    ]


def test_source_resolution_uses_only_pinned_default(monkeypatch):
    monkeypatch.delenv("GVCLASS_DB_URL", raising=False)

    sources = DatabaseManager._resolve_database_sources()

    assert sources == [DatabaseManager.DEFAULT_DATABASE_SOURCE]


def test_preferred_source_without_checksum_is_rejected(monkeypatch):
    monkeypatch.delenv("GVCLASS_DB_URL", raising=False)
    preferred = {"url": "https://example.com/config.tar.gz"}

    with pytest.raises(ValueError, match="no SHA-256 checksum"):
        DatabaseManager._resolve_database_sources(preferred)


@pytest.mark.parametrize("checksum", ["", "abc123", "g" * 64])
def test_zenodo_record_requires_valid_sha256_sidecar(checksum, monkeypatch):
    record = {
        "id": 123,
        "metadata": {"version": "2.1.0"},
        "files": [
            {
                "key": "resources.tar.gz",
                "links": {"self": "https://zenodo.org/archive"},
            },
            {
                "key": "resources.tar.gz.sha256",
                "links": {"self": "https://zenodo.org/checksum"},
            },
        ],
    }
    monkeypatch.setattr(
        DatabaseManager, "_fetch_text", staticmethod(lambda url: checksum)
    )

    with pytest.raises(ValueError, match="invalid SHA-256 sidecar"):
        DatabaseManager._source_from_zenodo_record(record)


def test_zenodo_record_requires_sha256_sidecar():
    record = {
        "id": 123,
        "metadata": {"version": "2.1.0"},
        "files": [
            {
                "key": "resources.tar.gz",
                "links": {"self": "https://zenodo.org/archive"},
            }
        ],
    }

    with pytest.raises(ValueError, match="no SHA-256 sidecar"):
        DatabaseManager._source_from_zenodo_record(record)
