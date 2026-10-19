"""Shared helpers: a small resource bundle with a manifest."""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Mapping
from pathlib import Path

import pytest

from src.contig_screen import layout

BundleFactory = Callable[[Path, Mapping[str, bytes]], layout.ScreenResources]


def write_bundle(root: Path, files: Mapping[str, bytes]) -> str:
    """Write ``files`` under ``root`` with a manifest; return the manifest digest."""
    for relative, content in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    manifest = layout.build_manifest(root, files)
    (root / layout.MANIFEST_FILE).parent.mkdir(parents=True, exist_ok=True)
    (root / layout.MANIFEST_FILE).write_bytes(manifest)
    return hashlib.sha256(manifest).hexdigest()


@pytest.fixture
def make_resources() -> BundleFactory:
    """Return a factory that writes a bundle and opens it as the screen does."""

    def factory(root: Path, files: Mapping[str, bytes]) -> layout.ScreenResources:
        digest = write_bundle(root, files)
        return layout.ScreenResources.open(root, manifest_sha256=digest)

    return factory
