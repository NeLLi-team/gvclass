"""Build a GVClass resource tree without listed prokaryotic reference proteins.

Reference proteins named in the removal tables leave ``parquet/faa.parquet``. A
reference genome that keeps no protein moves from the active to the inactive
label table. The script writes the removal evidence into ``provenance/v2.1.0/``
and rebuilds ``RESOURCE_MANIFEST.json`` and the resource metadata. Unchanged
files are hard-linked from the source tree. The 2.1.0 resource was built with:

    python scripts/clean_reference_resource.py --source <staged v2.1.0 tree> \
        --out <cleaned tree> --remove <removal tables> --members <members.tsv>
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Final

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

logger = logging.getLogger(__name__)

FAA: Final = Path("parquet/faa.parquet")
LABELS: Final = Path("parquet/labels/labels.parquet")
INACTIVE: Final = Path("parquet/labels/inactive_labels.parquet")
CONTEXT: Final = Path("parquet/labels/label_context.parquet")
PROVENANCE: Final = Path("provenance/v2.1.0")
REMOVED_TABLE: Final = PROVENANCE / "removed_reference_sequences.tsv"
SUMMARY_TABLE: Final = PROVENANCE / "reference_cleaning_summary.tsv"
INACTIVE_REASON: Final = (
    "all_reference_sequences_removed_prokaryotic_in_eukaryote_v2_1_0"
)
REWRITTEN: Final = frozenset(
    {
        FAA,
        LABELS,
        INACTIVE,
        CONTEXT,
        Path("RESOURCE_METADATA.json"),
        Path("RESOURCE_MANIFEST.json"),
    }
)
BATCH_ROWS: Final = 100_000
# The v2.0.0 build wrote its Parquet files at zstd level 9. The pyarrow default
# (level 1) made faa.parquet 7% larger for the same rows.
ZSTD_LEVEL: Final = 9
NEW_MEMBERS: Final = (
    (
        REMOVED_TABLE,
        "provenance",
        "Reference proteins removed from the eukaryotic references, with evidence.",
    ),
    (
        SUMMARY_TABLE,
        "provenance",
        "Counts of reference proteins removed per reference genome and source.",
    ),
)


def main(argv: Sequence[str] | None = None) -> int:
    """Build the cleaned resource tree and print its identity."""
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    removal = read_removal_tables(args.remove)
    logger.info("%d proteins to remove", len(removal))
    if args.out.exists():
        message = f"output tree exists: {args.out}"
        raise SystemExit(message)
    link_unchanged(args.source, args.out)
    removed_rows = rewrite_faa(args.source, args.out, removal)
    orphaned, lost_pairs = update_labels(args.source, args.out, removed_rows)
    write_provenance(args.out, removal, removed_rows, orphaned)
    update_metadata(args.source, args.out, removal, removed_rows, orphaned, lost_pairs)
    members = write_manifest(args.source, args.out, args.members)
    logger.info("manifest members: %d", members)
    return 0


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--remove",
        type=Path,
        nargs="+",
        required=True,
        help="TSV with seq_id, source and evidence columns",
    )
    parser.add_argument(
        "--members",
        type=Path,
        required=True,
        help="members TSV of the source tree (path, category, purpose)",
    )
    return parser.parse_args(argv)


def read_removal_tables(paths: Iterable[Path]) -> dict[str, dict[str, str]]:
    """Return one record per protein; the first table that lists it wins."""
    removal: dict[str, dict[str, str]] = {}
    for path in paths:
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                removal.setdefault(
                    row["seq_id"],
                    {"source": row["source"], "evidence": row["evidence"]},
                )
    return removal


def link_unchanged(source: Path, out: Path) -> None:
    """Hard-link every file the cleaning does not rewrite."""
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        target = out / relative
        if path.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        elif relative not in REWRITTEN:
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(path, target)
            except OSError:
                shutil.copy2(path, target)


def rewrite_faa(
    source: Path, out: Path, removal: dict[str, dict[str, str]]
) -> list[dict[str, str]]:
    """Write faa.parquet without the listed proteins; return the removed rows."""
    wanted = pa.array(sorted(removal), type=pa.string())
    reader = pq.ParquetFile(source / FAA)
    target = out / FAA
    target.parent.mkdir(parents=True, exist_ok=True)
    removed: list[dict[str, str]] = []
    with pq.ParquetWriter(
        target,
        reader.schema_arrow,
        compression="zstd",
        compression_level=ZSTD_LEVEL,
        write_statistics=True,
    ) as writer:
        for batch in reader.iter_batches(batch_size=BATCH_ROWS):
            table = pa.Table.from_batches([batch])
            hit = pc.is_in(table["seq_id"], value_set=wanted)
            if pc.any(hit).as_py():
                dropped = table.filter(hit).select(
                    ["seq_id", "label_id", "domain", "marker"]
                )
                removed.extend(dropped.to_pylist())
                table = table.filter(pc.invert(hit))
            if table.num_rows:
                writer.write_table(table)
    logger.info("removed %d reference rows", len(removed))
    return removed


def update_labels(
    source: Path, out: Path, removed_rows: list[dict[str, str]]
) -> tuple[list[str], int]:
    """Move genomes without any protein left to the inactive label table.

    Also returns how many (genome, marker) pairs lost their only protein, which
    the marker conservation tables do not reflect.
    """
    touched = {row["label_id"] for row in removed_rows}
    table = pq.read_table(out / FAA, columns=["label_id", "marker"])
    remaining = set(pc.unique(table["label_id"]).to_pylist())
    orphaned = sorted(touched - remaining)
    kept = table.filter(
        pc.is_in(table["label_id"], value_set=pa.array(sorted(touched), pa.string()))
    )
    kept_pairs = set(
        zip(kept["label_id"].to_pylist(), kept["marker"].to_pylist(), strict=True)
    )
    lost_pairs = {(row["label_id"], row["marker"]) for row in removed_rows} - kept_pairs
    labels = pq.read_table(source / LABELS)
    inactive = pq.read_table(source / INACTIVE)
    context = pq.read_table(source / CONTEXT)
    if orphaned:
        gone = pa.array(orphaned, type=pa.string())
        moved = labels.filter(pc.is_in(labels["label_id"], value_set=gone))
        rows = [
            {
                "label_id": row["label_id"],
                "taxonomy": row["taxonomy"],
                "reason": INACTIVE_REASON,
            }
            for row in moved.to_pylist()
        ]
        inactive = pa.concat_tables(
            [inactive, pa.Table.from_pylist(rows, schema=inactive.schema)]
        )
        labels = labels.filter(pc.invert(pc.is_in(labels["label_id"], value_set=gone)))
        context = context.filter(
            pc.invert(pc.is_in(context["label_id"], value_set=gone))
        )
    (out / LABELS).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        labels,
        out / LABELS,
        compression="zstd",
        compression_level=ZSTD_LEVEL,
    )
    pq.write_table(
        inactive,
        out / INACTIVE,
        compression="zstd",
        compression_level=ZSTD_LEVEL,
    )
    pq.write_table(
        context,
        out / CONTEXT,
        compression="zstd",
        compression_level=ZSTD_LEVEL,
    )
    return orphaned, len(lost_pairs)


def write_provenance(
    out: Path,
    removal: dict[str, dict[str, str]],
    removed_rows: list[dict[str, str]],
    orphaned: list[str],
) -> None:
    """Write the removed-protein table and the per-genome summary."""
    (out / PROVENANCE).mkdir(parents=True, exist_ok=True)
    by_id: dict[str, dict[str, set[str]]] = {}
    for row in removed_rows:
        entry = by_id.setdefault(
            row["seq_id"],
            {"label_id": set(), "domain": set(), "markers": set()},
        )
        entry["label_id"].add(row["label_id"])
        entry["domain"].add(row["domain"])
        entry["markers"].add(row["marker"])
    with (out / REMOVED_TABLE).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["seq_id", "label_id", "domain", "markers", "source", "evidence"]
        )
        for seq_id in sorted(by_id):
            entry = by_id[seq_id]
            writer.writerow(
                [
                    seq_id,
                    ",".join(sorted(entry["label_id"])),
                    ",".join(sorted(entry["domain"])),
                    ",".join(sorted(entry["markers"])),
                    removal[seq_id]["source"],
                    removal[seq_id]["evidence"],
                ]
            )
    counts: Counter[tuple[str, str]] = Counter()
    for seq_id, entry in by_id.items():
        for label_id in entry["label_id"]:
            counts[(label_id, removal[seq_id]["source"])] += 1
    with (out / SUMMARY_TABLE).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["label_id", "source", "removed_proteins", "genome_emptied"])
        for (label_id, source), count in sorted(counts.items()):
            writer.writerow([label_id, source, count, label_id in orphaned])


def update_metadata(
    source: Path,
    out: Path,
    removal: dict[str, dict[str, str]],
    removed_rows: list[dict[str, str]],
    orphaned: list[str],
    lost_pairs: int,
) -> None:
    """Record the cleaning in RESOURCE_METADATA.json."""
    metadata = json.loads((source / "RESOURCE_METADATA.json").read_text())
    by_source = Counter(
        removal[row["seq_id"]]["source"]
        for row in removed_rows
        if row["seq_id"] in removal
    )
    removed_table = out / REMOVED_TABLE
    metadata["v2_1_0_reference_cleaning"] = {
        "reason": (
            "Prokaryotic proteins inside eukaryotic reference genomes gave "
            "prokaryotic genomes a eukaryotic domain call."
        ),
        "removed_unique_proteins": len({row["seq_id"] for row in removed_rows}),
        "removed_reference_rows": len(removed_rows),
        "removed_by_source": dict(sorted(by_source.items())),
        "genomes_emptied": orphaned,
        "genome_marker_pairs_without_protein": lost_pairs,
        "removed_table": REMOVED_TABLE.as_posix(),
        "removed_table_sha256": hashlib.sha256(removed_table.read_bytes()).hexdigest(),
        "unchanged_by_cleaning": (
            "markers/*, completeness/* and the model files were not recomputed."
        ),
    }
    metadata["download_archive_status"] = (
        "Release archive of GVClass reference resources v2.1.0 on Zenodo "
        "(concept DOI 10.5281/zenodo.18662445); the archive checksum is on the record."
    )
    metadata["v2_1_0_staging_status"] = "Release archive."
    (out / "RESOURCE_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n")


def write_manifest(source: Path, out: Path, members_table: Path) -> int:
    """Rebuild RESOURCE_MANIFEST.json for the cleaned tree."""
    old = json.loads((source / "RESOURCE_MANIFEST.json").read_text())
    rows = members_table.read_text().splitlines()
    if rows[0] != "path\tcategory\tpurpose":
        message = f"unexpected members header in {members_table}"
        raise SystemExit(message)
    entries = [tuple(row.split("\t", 2)) for row in rows[1:]]
    known = {path for path, _, _ in entries}
    for path, category, purpose in NEW_MEMBERS:
        if path.as_posix() not in known:
            entries.append((path.as_posix(), category, purpose))
    members = []
    for relative, category, purpose in entries:
        if relative == "RESOURCE_MANIFEST.json":
            continue
        path = out / relative
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        members.append(
            {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": digest,
                "category": category,
                "purpose": purpose,
            }
        )
    manifest = {**old, "member_count": len(members), "members": members}
    (out / "RESOURCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    observed = sorted(
        p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()
    )
    expected = sorted({e[0] for e in entries} | {"RESOURCE_MANIFEST.json"})
    if observed != expected:
        message = f"tree differs from the member list: {set(observed) ^ set(expected)}"
        raise SystemExit(message)
    out_members = out.parent / f"{out.name}.members.tsv"
    lines = ["path\tcategory\tpurpose"]
    lines.extend("\t".join(entry) for entry in entries)
    out_members.write_text("\n".join(lines) + "\n")
    return len(members)


if __name__ == "__main__":
    raise SystemExit(main())
