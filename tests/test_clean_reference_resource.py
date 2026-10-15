"""Tests for the reference-resource cleaning script."""

from __future__ import annotations

import importlib.util
import json
import random
from pathlib import Path
from types import ModuleType

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "clean_reference_resource.py"
MEMBERS = (
    "path\tcategory\tpurpose\n"
    "DB_VERSION\tidentity\tversion\n"
    "parquet/faa.parquet\truntime\tproteins\n"
    "parquet/labels/labels.parquet\truntime\tlabels\n"
    "parquet/labels/inactive_labels.parquet\tcuration\tinactive\n"
    "parquet/labels/label_context.parquet\tcuration\tcontext\n"
    "markers/stats.tsv\truntime\tstats\n"
    "RESOURCE_METADATA.json\tidentity\tmetadata\n"
    "RESOURCE_MANIFEST.json\tmanifest\tmanifest\n"
)


def _load_script() -> ModuleType:
    spec = importlib.util.spec_from_file_location("clean_reference_resource", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _compressible_rows(count: int) -> list[dict[str, str]]:
    """Return rows of point-mutated copies of a few proteins.

    Such sequences compress far better at zstd level 9 than at the default level,
    so a rewrite at the default level is visible in the file size.
    """
    rng = random.Random(11)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    bases = ["".join(rng.choices(alphabet, k=350)) for _ in range(300)]

    def mutated(sequence: str) -> str:
        residues = list(sequence)
        for position in rng.sample(range(len(residues)), 12):
            residues[position] = rng.choice(alphabet)
        return "".join(residues)

    return [
        {
            "marker": f"M{index % 5}.faa",
            "seq_id": f"EUK__G{index % 50}|p{index}",
            "label_id": f"EUK__G{index % 50}",
            "domain": "EUK",
            "sequence": mutated(rng.choice(bases)),
        }
        for index in range(count)
    ]


def _row(marker: str, seq_id: str, domain: str = "EUK") -> dict[str, str]:
    return {
        "marker": marker,
        "seq_id": seq_id,
        "label_id": seq_id.split("|", 1)[0],
        "domain": domain,
        "sequence": "MKT",
    }


@pytest.fixture
def source_tree(tmp_path: Path) -> Path:
    tree = tmp_path / "source"
    (tree / "parquet" / "labels").mkdir(parents=True)
    (tree / "markers").mkdir()
    rows = [
        _row("M1.faa", "EUK__A|a1"),
        _row("M2.faa", "EUK__A|a1"),
        _row("M1.faa", "EUK__A|a2"),
        _row("M1.faa", "EUK__B|b1"),
        _row("M1.faa", "BAC__C|c1", "BAC"),
    ]
    pq.write_table(pa.Table.from_pylist(rows), tree / "parquet" / "faa.parquet")
    labels = [{"label_id": x, "taxonomy": "t"} for x in ("EUK__A", "EUK__B", "BAC__C")]
    pq.write_table(
        pa.Table.from_pylist(labels), tree / "parquet" / "labels" / "labels.parquet"
    )
    inactive = [{"label_id": "OLD", "taxonomy": "t", "reason": "earlier"}]
    pq.write_table(
        pa.Table.from_pylist(inactive),
        tree / "parquet" / "labels" / "inactive_labels.parquet",
    )
    context = [
        {"label_id": "EUK__B", "context_type": "c", "context_value": "v", "source": "s"}
    ]
    pq.write_table(
        pa.Table.from_pylist(context),
        tree / "parquet" / "labels" / "label_context.parquet",
    )
    (tree / "markers" / "stats.tsv").write_text("marker\tn\nM1\t3\n", encoding="utf-8")
    (tree / "DB_VERSION").write_text("v2.1.0\n", encoding="utf-8")
    (tree / "RESOURCE_METADATA.json").write_text('{"version": "v2.1.0"}', "utf-8")
    (tree / "RESOURCE_MANIFEST.json").write_text('{"schema": "x"}', encoding="utf-8")
    return tree


def _run(tmp_path: Path, tree: Path) -> Path:
    members = tmp_path / "members.tsv"
    members.write_text(MEMBERS, encoding="utf-8")
    removal = tmp_path / "remove.tsv"
    removal.write_text(
        "seq_id\tsource\tevidence\n"
        "EUK__A|a1\tlisted\te1\n"
        "EUK__B|b1\tscreen\te2\n"
        "EUK__ABSENT|x\tscreen\te3\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    code = _load_script().main(
        [
            "--source",
            str(tree),
            "--out",
            str(out),
            "--remove",
            str(removal),
            "--members",
            str(members),
        ]
    )
    assert code == 0
    return out


def test_cleaning_drops_listed_proteins_under_every_marker(
    tmp_path: Path, source_tree: Path
) -> None:
    out = _run(tmp_path, source_tree)

    remaining = pq.read_table(out / "parquet" / "faa.parquet").column("seq_id")
    assert sorted(remaining.to_pylist()) == ["BAC__C|c1", "EUK__A|a2"]


def test_cleaning_moves_an_emptied_genome_to_the_inactive_table(
    tmp_path: Path, source_tree: Path
) -> None:
    out = _run(tmp_path, source_tree)

    labels = pq.read_table(out / "parquet" / "labels" / "labels.parquet")
    inactive = pq.read_table(out / "parquet" / "labels" / "inactive_labels.parquet")
    context = pq.read_table(out / "parquet" / "labels" / "label_context.parquet")
    assert sorted(labels.column("label_id").to_pylist()) == ["BAC__C", "EUK__A"]
    assert inactive.column("label_id").to_pylist() == ["OLD", "EUK__B"]
    assert context.num_rows == 0


def test_cleaning_records_provenance_metadata_and_manifest(
    tmp_path: Path, source_tree: Path
) -> None:
    out = _run(tmp_path, source_tree)

    metadata = json.loads((out / "RESOURCE_METADATA.json").read_text("utf-8"))
    cleaning = metadata["v2_1_0_reference_cleaning"]
    assert cleaning["removed_unique_proteins"] == 2
    assert cleaning["removed_reference_rows"] == 3
    assert cleaning["genomes_emptied"] == ["EUK__B"]
    assert cleaning["genome_marker_pairs_without_protein"] == 2
    assert metadata["v2_1_0_staging_status"] == "Release archive."
    assert "zenodo" not in metadata["download_archive_status"].lower()
    manifest = json.loads((out / "RESOURCE_MANIFEST.json").read_text("utf-8"))
    listed = {member["path"] for member in manifest["members"]}
    assert "provenance/v2.1.0/removed_reference_sequences.tsv" in listed
    assert manifest["member_count"] == len(listed)
    removed = (out / "provenance/v2.1.0/removed_reference_sequences.tsv").read_text(
        "utf-8"
    )
    assert "EUK__A|a1\tEUK__A\tEUK\tM1.faa,M2.faa\tlisted\te1" in removed


def test_cleaning_links_unchanged_files_and_refuses_an_existing_output(
    tmp_path: Path, source_tree: Path
) -> None:
    out = _run(tmp_path, source_tree)

    assert (out / "markers" / "stats.tsv").samefile(source_tree / "markers/stats.tsv")
    with pytest.raises(SystemExit, match="output tree exists"):
        _run(tmp_path, source_tree)


def test_cleaning_keeps_the_compression_level_of_the_source(
    tmp_path: Path, source_tree: Path
) -> None:
    source = source_tree / "parquet" / "faa.parquet"
    pq.write_table(
        pa.Table.from_pylist(_compressible_rows(20_000)),
        source,
        compression="zstd",
        compression_level=9,
    )

    out = _run(tmp_path, source_tree)

    assert (
        out / "parquet" / "faa.parquet"
    ).stat().st_size <= source.stat().st_size * 1.01
