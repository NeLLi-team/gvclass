#!/usr/bin/env python3
"""Generate the Marker annotations reference page from the resources database.

Reads ``resources/markers/annotations.tsv`` (the majority functional annotation
per marker) and ``resources/hmm/combined.hmm`` (to count the markers that are
models of the search library) and writes two committed files under ``docs/reference/``:

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
HMM_FILE = REPO / "resources" / "hmm" / "combined.hmm"
OUT_TSV = REPO / "docs" / "reference" / "marker-annotations.tsv"
OUT_MD = REPO / "docs" / "reference" / "marker-annotations.md"


def _esc(text: str) -> str:
    """Make a value safe to drop into a Markdown table cell."""
    return text.replace("|", r"\|").strip()


def _hmm_model_names(path: Path) -> set[str]:
    """Return the model names of an HMMER profile library."""
    with path.open(encoding="utf-8") as handle:
        return {line.split()[1] for line in handle if line.startswith("NAME ")}


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

    hmm_names = _hmm_model_names(HMM_FILE)
    n_models = sum(row["marker_id"] in hmm_names for row in rows)

    lines = [
        "# Marker annotations",
        "",
        (
            f"The table lists the {len(rows)} markers of `markers/annotations.tsv` in "
            f"database {db_version}, each with the majority functional annotation of "
            "its member proteins. A marker is an HMM model name or an orthogroup ID "
            f"(`OG<number>`). Of the {len(rows)} markers, {n_models} are models of "
            f"the search library `hmm/combined.hmm`, which holds {len(hmm_names)} "
            "models."
        ),
        "",
        (
            "The [TSV download](marker-annotations.tsv) includes support counts "
            "(`majority_count`, `total_proteins`). See "
            "[marker panels and genetic codes](markers.md) for the summary metrics."
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
