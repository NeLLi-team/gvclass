"""Opt-in resource audit for labels.tsv and reference FAA headers."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Iterator

import pytest

from src.utils.resource_store import ResourceStore


REPO_ROOT = Path(__file__).resolve().parents[1]
RESOURCE_DIR = Path(
    os.environ.get("GVCLASS_RESOURCE_AUDIT_DIR", REPO_ROOT / "resources")
)
STORE = ResourceStore(RESOURCE_DIR)
FAA_DIR = Path(
    os.environ.get(
        "GVCLASS_RESOURCE_AUDIT_FAA_DIR",
        RESOURCE_DIR / "database" / "faa",
    )
)


def _load_label_ids() -> set[str]:
    labels = STORE.label_path("labels.tsv")
    label_ids: set[str] = set()
    duplicates: set[str] = set()
    with labels.open() as handle:
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
    aliases_path = STORE.label_path("aliases.tsv")
    if not aliases_path.exists():
        return aliases

    with aliases_path.open() as handle:
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


def _iter_faa_header_ids() -> Iterator[str]:
    if FAA_DIR.exists():
        for faa in sorted(FAA_DIR.glob("*.faa")):
            with faa.open() as handle:
                for line in handle:
                    if line.startswith(">"):
                        yield line[1:].split()[0]
        return

    if STORE.has_parquet_faa():
        pq, _ds = STORE._pyarrow_modules()
        parquet = pq.ParquetFile(STORE.faa_parquet)
        for batch in parquet.iter_batches(columns=["seq_id"]):
            yield from batch.column(0).to_pylist()


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
    if not STORE.labels_satisfied(RESOURCE_DIR) or not (
        FAA_DIR.exists() or STORE.has_parquet_faa()
    ):
        pytest.skip("GVClass resource bundle is not installed")

    label_ids = _load_label_ids()
    aliases = _load_aliases(label_ids)
    resolved_labels: set[str] = set()
    unresolved_count = 0
    unresolved_examples: list[str] = []

    for header_id in _iter_faa_header_ids():
        resolved = _resolve_header_id(header_id, label_ids, aliases)
        if resolved is None:
            unresolved_count += 1
            if len(unresolved_examples) < 20:
                unresolved_examples.append(header_id)
        else:
            resolved_labels.add(resolved)

    labels_without_faa = sorted(label_ids - resolved_labels)

    assert unresolved_count == 0, (
        f"{unresolved_count} FAA headers did not resolve to labels; "
        f"examples: {unresolved_examples}"
    )
    assert not labels_without_faa, (
        f"{len(labels_without_faa)} labels had no FAA evidence; "
        f"examples: {labels_without_faa[:20]}"
    )
