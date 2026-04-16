"""Golden-file regression test for ``pixi run example``.

This test runs the bundled ``example/`` directory through GVClass end-to-end
and diffs the stable columns of the resulting ``gvclass_summary.tsv``
against ``tests/fixtures/expected_example_summary.tsv``. Any change to
classification, completeness, or the columns the fixture was built on
trips the test — catching regressions that the unit tests miss.

Skipped unless the reference database is installed at ``resources/`` and
the ``GVCLASS_RUN_GOLDEN`` environment variable is set, so CI jobs
opt-in explicitly (the DB download is large and slow; see the
corresponding CI step in ``.github/workflows/ci.yml``).
"""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import pytest


STABLE_COLUMNS = [
    "query",
    "taxonomy_majority",
    "taxonomy_confidence",
    "estimated_completeness",
    "estimated_completeness_quality",
    "estimated_contamination",
    "contamination_type",
    "gvog4_unique",
    "gvog8_unique",
    "ncldv_mcp_total",
    "vp_completeness",
    "mirus_completeness",
    "ttable",
]


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "expected_example_summary.tsv"
DB_PATH = REPO_ROOT / "resources"


def _project_stable_columns(path: Path) -> list[dict[str, str]]:
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = sorted(reader, key=lambda r: r["query"])
    projected = [
        {col: row.get(col, "") for col in STABLE_COLUMNS}
        for row in rows
    ]
    return projected


@pytest.mark.slow
@pytest.mark.requires_db
@pytest.mark.skipif(
    os.environ.get("GVCLASS_RUN_GOLDEN") != "1",
    reason="Set GVCLASS_RUN_GOLDEN=1 to run the full example-pipeline regression test",
)
def test_example_summary_matches_golden_fixture(tmp_path: Path) -> None:
    assert FIXTURE_PATH.exists(), (
        f"Missing fixture at {FIXTURE_PATH}. Regenerate with:\n"
        "  pixi run gvclass example -o /tmp/gv -t 8 --mode-fast"
    )
    if not (DB_PATH / "models" / "combined.hmm").exists():
        pytest.skip(
            f"GVClass database not installed at {DB_PATH}; "
            "set it up with `pixi run setup-db` before running this test."
        )

    output_dir = tmp_path / "results"
    subprocess.run(
        [
            "pixi",
            "run",
            "gvclass",
            "example",
            "-o",
            str(output_dir),
            "-t",
            "4",
            "--mode-fast",
            "--plain-output",
        ],
        cwd=str(REPO_ROOT),
        check=True,
    )

    produced_rows = _project_stable_columns(output_dir / "gvclass_summary.tsv")
    expected_rows = _project_stable_columns(FIXTURE_PATH)

    assert produced_rows == expected_rows, (
        "Example summary drift detected. Inspect diff and regenerate the "
        f"fixture if the change is intentional:\n{FIXTURE_PATH}"
    )
