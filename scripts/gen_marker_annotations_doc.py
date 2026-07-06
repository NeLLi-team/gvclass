#!/usr/bin/env python3
"""Generate the Marker annotations reference page from the resources database.

Reads ``resources/markers/annotations.tsv`` (the majority functional annotation
per HMM marker model) and writes two committed files under ``docs/reference/``:

  - ``marker-annotations.tsv`` -- verbatim snapshot (all columns) for download
  - ``marker-annotations.md``  -- rendered two-column lookup table for the site

``resources/`` is gitignored, so these committed snapshots are what the
published docs serve. Re-run after a database bundle update::

    python scripts/gen_marker_annotations_doc.py

The page is stamped with the bundle name in ``resources/DB_VERSION``.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC_TSV = REPO / "resources" / "markers" / "annotations.tsv"
DB_VERSION_FILE = REPO / "resources" / "DB_VERSION"
OUT_TSV = REPO / "docs" / "reference" / "marker-annotations.tsv"
OUT_MD = REPO / "docs" / "reference" / "marker-annotations.md"


def _esc(text: str) -> str:
    """Make a value safe to drop into a Markdown table cell."""
    return text.replace("|", r"\|").strip()


def main() -> None:
    if not SRC_TSV.exists():
        raise SystemExit(
            f"Source not found: {SRC_TSV}\n"
            "Populate resources/ first (e.g. `pixi run setup-db`)."
        )
    db_version = (
        DB_VERSION_FILE.read_text().strip() if DB_VERSION_FILE.exists() else "unknown"
    )

    with SRC_TSV.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # Verbatim snapshot (all columns) for download.
    shutil.copyfile(SRC_TSV, OUT_TSV)

    lines = [
        "# Marker annotations",
        "",
        (
            f"Majority functional annotation for each of the {len(rows)} HMM marker "
            "models in GVClass, derived from the consensus annotation of the member "
            f"proteins in each model. This is a snapshot of database bundle "
            f"**{db_version}**."
        ),
        "",
        (
            "Download the full table, including per-annotation support counts "
            "(`majority_count`, `total_proteins`), here: "
            "[marker-annotations.tsv](marker-annotations.tsv). Regenerate it with "
            "`python scripts/gen_marker_annotations_doc.py` after a database update. "
            "For the marker panels that feed the summary table and the genetic codes, "
            "see [Marker panels and genetic codes](markers.md)."
        ),
        "",
        "| Marker model | Majority annotation |",
        "| --- | --- |",
    ]
    lines.extend(
        f"| `{_esc(r['marker_id'])}` | {_esc(r['majority_annotation'])} |" for r in rows
    )
    lines.append("")
    OUT_MD.write_text("\n".join(lines))
    print(f"Wrote {OUT_MD} and {OUT_TSV} ({len(rows)} models, bundle {db_version})")


if __name__ == "__main__":
    main()
