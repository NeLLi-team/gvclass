"""Compare a SIF run with the golden summary and check the plain-output log.

Runs inside the SIF with the golden test file mounted at ``/golden_test.py``. The
SIF has no pytest, so the column list is read from the test file without
importing it; the test file stays the one place that names the stable columns.
"""

import ast
import csv
from pathlib import Path

from src.bin.gvclass_cli import CliOutput


def read_stable_columns(test_file):
    """Return the ``STABLE_COLUMNS`` list of the golden-file test."""
    for node in ast.parse(Path(test_file).read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "STABLE_COLUMNS"
            for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise SystemExit(f"STABLE_COLUMNS not found in {test_file}")


stable_columns = read_stable_columns("/golden_test.py")


def stable_summary(path):
    with Path(path).open() as handle:
        return sorted(
            tuple(row[column] for column in stable_columns)
            for row in csv.DictReader(handle, delimiter="\t")
        )


expected = stable_summary("/expected.tsv")
observed = stable_summary("/output/gvclass_summary.tsv")
assert expected, "The golden fixture must contain classifications"
assert observed == expected, (observed, expected)

classification_log = Path("/validation/classification.log").read_text(encoding="utf-8")
assert "\x1b" not in classification_log, "ANSI escape in plain output"
decorations = (*CliOutput.ICONS.values(), "⚠", "✓")
assert not any(symbol in classification_log for symbol in decorations), (
    "Decoration in plain output"
)
assert "Running pyswrd search:" in classification_log, "Missing search diagnostics"
assert " - DEBUG - Processing marker " in classification_log, (
    "Missing verbose marker diagnostics"
)
