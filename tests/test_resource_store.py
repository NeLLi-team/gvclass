from __future__ import annotations

import os
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from src.core.marker_processing import MarkerProcessor
from src.utils.resource_store import ResourceStore

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")


def _write_parquet(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


def test_resource_store_prefers_legacy_labels_and_faa(tmp_path: Path) -> None:
    labels = tmp_path / "labels.tsv"
    labels.write_text("NCLDV__ref\tNCLDV|p|c|o|f|g|s\n")
    faa = tmp_path / "database" / "faa" / "GVOGm0001.faa"
    faa.parent.mkdir(parents=True)
    faa.write_text(">NCLDV__ref_1\nMKT\n")

    store = ResourceStore(tmp_path)

    assert store.label_path("labels.tsv") == labels
    assert store.marker_faa_path("GVOGm0001") == faa


def test_resource_store_materializes_parquet_labels_and_marker_faa(
    tmp_path: Path,
) -> None:
    _write_parquet(
        tmp_path / "parquet" / "labels" / "labels.parquet",
        [
            {"label_id": "NCLDV__ref", "taxonomy": "NCLDV|p|c|o|f|g|s"},
            {
                "label_id": "EUK-pEVE__ref",
                "taxonomy": "EUK-pEVE|p-pEVE|c-pEVE|o-pEVE|f-pEVE|g-pEVE|s-pEVE",
            },
        ],
    )
    _write_parquet(
        tmp_path / "parquet" / "faa.parquet",
        [
            {
                "marker": "GVOGm0001.faa",
                "seq_id": "NCLDV__ref_1",
                "label_id": "NCLDV__ref",
                "domain": "NCLDV",
                "sequence": "MKT",
            },
            {
                "marker": "GVOGm0001.faa",
                "seq_id": "EUK-pEVE__ref_1",
                "label_id": "EUK-pEVE__ref",
                "domain": "EUK-pEVE",
                "sequence": "MKK",
            },
            {
                "marker": "GVOGm0002.faa",
                "seq_id": "NCLDV__other_1",
                "label_id": "NCLDV__other",
                "domain": "NCLDV",
                "sequence": "MAA",
            },
        ],
    )

    store = ResourceStore(tmp_path)

    labels = store.label_path("labels.tsv")
    assert labels.read_text().splitlines() == [
        "NCLDV__ref\tNCLDV|p|c|o|f|g|s",
        "EUK-pEVE__ref\tEUK-pEVE|p-pEVE|c-pEVE|o-pEVE|f-pEVE|g-pEVE|s-pEVE",
    ]

    faa = store.marker_faa_path("GVOGm0001")
    assert faa.read_text().splitlines() == [
        ">NCLDV__ref_1",
        "MKT",
        ">EUK-pEVE__ref_1",
        "MKK",
    ]

    materialized_dir = store.materialized_faa_dir(["GVOGm0001"])
    assert (materialized_dir / "GVOGm0001.faa").exists()


def test_resource_store_concurrent_materialization_uses_unique_temp_files(
    tmp_path: Path, monkeypatch
) -> None:
    _write_parquet(
        tmp_path / "parquet" / "labels" / "labels.parquet",
        [{"label_id": "NCLDV__ref", "taxonomy": "NCLDV|p|c|o|f|g|s"}],
    )
    store = ResourceStore(tmp_path)
    target = store._label_cache_dir() / "labels.tsv"
    real_replace = os.replace
    replace_barrier = threading.Barrier(2)

    def synchronized_replace(source, destination):
        if Path(destination) == target:
            replace_barrier.wait(timeout=5)
        return real_replace(source, destination)

    monkeypatch.setattr(os, "replace", synchronized_replace)
    with ThreadPoolExecutor(max_workers=2) as executor:
        paths = list(executor.map(lambda _: store.label_path("labels.tsv"), range(2)))

    assert paths == [target, target]
    assert target.read_text() == "NCLDV__ref\tNCLDV|p|c|o|f|g|s\n"
    assert target.stat().st_mode & 0o777 == 0o644
    assert list(target.parent.glob(f".{target.name}.*.tmp")) == []


def test_resource_store_materialized_dir_can_mix_legacy_and_parquet(
    tmp_path: Path,
) -> None:
    legacy = tmp_path / "database" / "faa" / "LEGACY.faa"
    legacy.parent.mkdir(parents=True)
    legacy.write_text(">LEGACY__ref_1\nMKT\n")
    _write_parquet(
        tmp_path / "parquet" / "faa.parquet",
        [
            {
                "marker": "PARQUET.faa",
                "seq_id": "NCLDV__ref_1",
                "label_id": "NCLDV__ref",
                "domain": "NCLDV",
                "sequence": "MKK",
            },
        ],
    )

    materialized_dir = ResourceStore(tmp_path).materialized_faa_dir(
        ["LEGACY", "PARQUET"]
    )

    assert (materialized_dir / "LEGACY.faa").read_text() == ">LEGACY__ref_1\nMKT\n"
    assert (materialized_dir / "PARQUET.faa").read_text() == ">NCLDV__ref_1\nMKK\n"


def test_resource_store_materialized_dir_can_skip_missing_markers(
    tmp_path: Path,
) -> None:
    _write_parquet(
        tmp_path / "parquet" / "faa.parquet",
        [
            {
                "marker": "PRESENT.faa",
                "seq_id": "NCLDV__ref_1",
                "label_id": "NCLDV__ref",
                "domain": "NCLDV",
                "sequence": "MKK",
            },
        ],
    )

    store = ResourceStore(tmp_path)
    with pytest.raises(FileNotFoundError):
        store.materialized_faa_dir(["PRESENT", "MISSING"])

    materialized_dir = store.materialized_faa_dir(["PRESENT", "MISSING"], strict=False)

    assert (materialized_dir / "PRESENT.faa").exists()
    assert not (materialized_dir / "MISSING.faa").exists()


def _faa_row(seq_id: str, sequence: str = "MKK") -> dict[str, str]:
    label_id = seq_id.split("|", 1)[0]
    return {
        "marker": "PRESENT.faa",
        "seq_id": seq_id,
        "label_id": label_id,
        "domain": "EUK",
        "sequence": sequence,
    }


def test_resource_store_rebuilds_cache_when_parquet_is_replaced_by_an_older_file(
    tmp_path: Path,
) -> None:
    parquet = tmp_path / "parquet" / "faa.parquet"
    _write_parquet(parquet, [_faa_row("EUK__contaminated|1")])
    cached = ResourceStore(tmp_path).marker_faa_path("PRESENT")
    assert cached.read_text() == ">EUK__contaminated|1\nMKK\n"

    # A cleaned archive keeps its build-time mtime, which is older than the cache.
    _write_parquet(
        parquet, [_faa_row("EUK__clean|1"), _faa_row("EUK__clean|2")]
    )
    older = cached.stat().st_mtime - 3600
    os.utime(parquet, (older, older))

    assert ResourceStore(tmp_path).marker_faa_path("PRESENT").read_text() == (
        ">EUK__clean|1\nMKK\n>EUK__clean|2\nMKK\n"
    )


def test_resource_store_keeps_trees_that_share_a_cache_apart(
    tmp_path: Path, monkeypatch
) -> None:
    shared = tmp_path / "shared-cache"
    monkeypatch.setenv("GVCLASS_RESOURCE_CACHE", str(shared))
    first, second = tmp_path / "first", tmp_path / "second"
    _write_parquet(first / "parquet" / "faa.parquet", [_faa_row("EUK__one|1")])
    _write_parquet(second / "parquet" / "faa.parquet", [_faa_row("EUK__two|1")])
    # Same size and the same modification time: only the manifest tells them apart.
    stamp = (first / "parquet" / "faa.parquet").stat().st_mtime
    os.utime(second / "parquet" / "faa.parquet", (stamp, stamp))
    (first / "RESOURCE_MANIFEST.json").write_text('{"members": [1]}')
    (second / "RESOURCE_MANIFEST.json").write_text('{"members": [2]}')

    one = ResourceStore(first).marker_faa_path("PRESENT")
    two = ResourceStore(second).marker_faa_path("PRESENT")

    assert one.read_text() == ">EUK__one|1\nMKK\n"
    assert two.read_text() == ">EUK__two|1\nMKK\n"
    assert one.parent != two.parent
    assert one.parent.parent.parent == shared


def test_marker_processor_preserves_missing_marker_fallback(tmp_path: Path) -> None:
    processor = MarkerProcessor(
        "MISSING",
        database_path=tmp_path / "db",
        output_dir=tmp_path / "out",
    )

    assert processor.get_database_path() == (
        tmp_path / "db" / "database" / "faa" / "MISSING.faa"
    )
