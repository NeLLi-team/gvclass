"""Opt-in resource audit for labels.tsv and reference FAA headers."""

from __future__ import annotations

import os
import re
from collections import defaultdict
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
RESOURCE_DIR = Path(
    os.environ.get("GVCLASS_RESOURCE_AUDIT_DIR", REPO_ROOT / "resources")
)
LABELS = RESOURCE_DIR / "labels.tsv"
ALIASES = RESOURCE_DIR / "aliases.tsv"
FAA_DIR = Path(
    os.environ.get(
        "GVCLASS_RESOURCE_AUDIT_FAA_DIR",
        RESOURCE_DIR / "database" / "faa",
    )
)


def _load_label_ids() -> set[str]:
    label_ids: set[str] = set()
    duplicates: set[str] = set()
    with LABELS.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            label_id = line.split("\t", 1)[0]
            if label_id in label_ids:
                duplicates.add(label_id)
            label_ids.add(label_id)
    assert duplicates == set()
    return label_ids


def _load_aliases(label_ids: set[str]) -> dict[str, str]:
    aliases: dict[str, str] = {}
    if os.environ.get("GVCLASS_RESOURCE_AUDIT_ALLOW_ALIASES", "1") == "0":
        return aliases
    if not ALIASES.exists():
        return aliases

    with ALIASES.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            alias_idx = header.index("alias_id")
            canonical_idx = header.index("canonical_id")
        except ValueError as exc:
            raise AssertionError(
                "aliases.tsv must include alias_id and canonical_id columns"
            ) from exc
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            alias_id = fields[alias_idx]
            canonical_id = fields[canonical_idx]
            aliases[alias_id] = canonical_id

    bad_targets = {
        alias_id: canonical_id
        for alias_id, canonical_id in aliases.items()
        if canonical_id not in label_ids
    }
    assert bad_targets == {}
    return aliases


def _candidate_label_ids(protein_id: str) -> list[str]:
    candidates = [protein_id]

    stripped = re.sub(r"_\d+$", "", protein_id)
    if stripped != protein_id:
        candidates.append(stripped)

    if "|" in protein_id:
        base, payload = protein_id.split("|", 1)
        candidates.append(base)
        contig = payload.rsplit("_", 1)[0] if "_" in payload else payload
        if contig:
            candidates.append(f"{base}_{contig}")

    return list(dict.fromkeys(candidates))


def _iter_faa_header_ids() -> list[str]:
    header_ids: list[str] = []
    for faa in sorted(FAA_DIR.glob("*.faa")):
        with faa.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    header_ids.append(line[1:].split()[0])
    return header_ids


def _resolve_header_id(
    header_id: str, label_ids: set[str], aliases: dict[str, str]
) -> str | None:
    for candidate in _candidate_label_ids(header_id):
        if candidate in label_ids:
            return candidate
        if candidate in aliases:
            return aliases[candidate]
    return None


@pytest.mark.requires_db
@pytest.mark.skipif(
    os.environ.get("GVCLASS_RUN_RESOURCE_AUDIT") != "1",
    reason="Set GVCLASS_RUN_RESOURCE_AUDIT=1 to audit the full resource bundle",
)
def test_reference_faa_headers_and_labels_are_bijective() -> None:
    """Every active label should have FAA evidence and every FAA header should resolve."""
    if not LABELS.exists() or not FAA_DIR.exists():
        pytest.skip("GVClass resource bundle is not installed")

    label_ids = _load_label_ids()
    aliases = _load_aliases(label_ids)
    headers = _iter_faa_header_ids()
    resolved_labels: dict[str, set[str]] = defaultdict(set)
    unresolved_headers: list[str] = []

    for header_id in headers:
        resolved = _resolve_header_id(header_id, label_ids, aliases)
        if resolved is None:
            unresolved_headers.append(header_id)
        else:
            resolved_labels[resolved].add(header_id)

    labels_without_faa = sorted(label_ids - set(resolved_labels))

    assert not unresolved_headers, (
        f"{len(unresolved_headers)} FAA headers did not resolve to labels; "
        f"examples: {unresolved_headers[:20]}"
    )
    assert not labels_without_faa, (
        f"{len(labels_without_faa)} labels had no FAA evidence; "
        f"examples: {labels_without_faa[:20]}"
    )
