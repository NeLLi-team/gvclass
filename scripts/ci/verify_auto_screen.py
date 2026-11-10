"""Run one installed-add-on classification, then reuse its screen twice."""

import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

SCREEN_FIELDS = {
    "screen_status",
    "screen_cellular_pct",
    "screen_giant_virus_pct",
    "screen_flagged_pct",
}
QUERY = "GVMAG-S-1096109-37"


def verify_plain_output(path):
    log = path.read_text(encoding="utf-8")
    assert "\x1b" not in log, "ANSI escape in plain output"
    decorations = (
        0x2705,
        0x274C,
        0x26A0,
        0x2713,
        0x1F4E6,
        0x1F4E5,
        0x1F4C2,
        0x1F50D,
        0x1F4CA,
        0x23F1,
        0x1F4BE,
    )
    assert not any(chr(code) in log for code in decorations), (
        "Decoration in plain output"
    )


def read_rows(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_pre_screen_summary(output):
    # The query writer saves this final-schema row before screening,
    # then archives it under the query directory. Read only that member.
    with tarfile.open(output / f"{QUERY}.tar.gz", "r:gz") as archive:
        member = archive.getmember(f"{QUERY}/{QUERY}.final_summary.tsv")
        assert member.isfile(), "Per-query summary must be a regular file"
        with archive.extractfile(member) as handle:
            lines = handle.read().decode("utf-8").splitlines()
            rows = list(csv.DictReader(lines, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["query"] == QUERY, rows
    return rows[0]


def verify_screen_totals(row, source, manifest):
    assert set(row) == {
        "bin",
        "cellular_status",
        "giant_virus_status",
        "input_sha256",
        "manifest_sha256",
        "input_bp",
        "cellular_bp",
        "giant_virus_bp",
        "flagged_bp",
        "cellular_pct",
        "giant_virus_pct",
        "flagged_pct",
    }, "Unexpected standalone summary columns"
    assert row["bin"] == QUERY
    assert row["cellular_status"] == row["giant_virus_status"] == "ok"
    assert row["input_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert row["manifest_sha256"] == hashlib.sha256(manifest.read_bytes()).hexdigest()
    assert int(row["input_bp"]) == 175691, row
    # This pinned fixture is all-negative, not a sensitivity test.
    for field in ("cellular", "giant_virus", "flagged"):
        assert int(row[f"{field}_bp"]) == 0, row
        assert float(row[f"{field}_pct"]) == 0, row


def verify_summary(path, baseline, standalone, expected_status):
    rows = read_rows(path)
    assert len(rows) == 1, rows
    row = rows[0]
    assert row["query"] == QUERY, row
    assert {key for key in row if key.startswith("screen_")} == SCREEN_FIELDS
    assert {key: value for key, value in row.items() if key not in SCREEN_FIELDS} == {
        key: value for key, value in baseline.items() if key not in SCREEN_FIELDS
    }, "Automatic screening changed classification fields"
    assert row["screen_status"] == expected_status, row
    for field in sorted(SCREEN_FIELDS - {"screen_status"}):
        observed = float(row[field])
        if expected_status == "disabled":
            assert math.isnan(observed), (field, row[field])
        else:
            expected = float(standalone[field.removeprefix("screen_")])
            assert math.isfinite(observed) and 0 <= observed <= 100
            assert observed == expected, (field, observed, expected)


def screen_snapshot(folder):
    snapshot = {}
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            stat = path.stat()
            snapshot[str(path.relative_to(folder))] = {
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "mtime_ns": stat.st_mtime_ns,
                "inode": stat.st_ino,
            }
    assert snapshot, "No automatic screen artifacts"
    assert f"{QUERY}/contig_screen_summary.tsv" in snapshot
    assert f"{QUERY}/screen.log" in snapshot
    return snapshot


def query_snapshot(output):
    # Include the archive and the loose sequence/summary sidecars,
    # but not the global summary or the separate screen directory.
    snapshot = {}
    for path in sorted(output.rglob(f"{QUERY}.*")):
        if path.is_file() and "contig_screen" not in path.relative_to(output).parts:
            stat = path.stat()
            with path.open("rb") as handle:
                digest = hashlib.file_digest(handle, "sha256").hexdigest()
            snapshot[str(path.relative_to(output))] = [
                digest,
                stat.st_mtime_ns,
                stat.st_ino,
            ]
    assert f"{QUERY}.tar.gz" in snapshot, "Missing query archive"
    return snapshot


temporary, repository = map(Path, sys.argv[1:])
evidence = temporary / "gvclass-validation"
inputs = temporary / "gvclass-auto-screen-input"
inputs.mkdir()
shutil.copy2(repository / "example" / f"{QUERY}.fna", inputs)
output = temporary / "gvclass-auto-screen-results"
earlier = next(
    row
    for row in read_rows(temporary / "gvclass-results/gvclass_summary.tsv")
    if row["query"] == QUERY
)
classification_columns = set(earlier) - SCREEN_FIELDS
baseline = None
standalone_rows = read_rows(
    temporary / "gvclass-screen-results/contig_screen_summary.tsv"
)
assert len(standalone_rows) == 1
standalone = standalone_rows[0]
verify_screen_totals(
    standalone,
    inputs / f"{QUERY}.fna",
    temporary / "gvclass-screen-cache/contig_screen/manifest.json",
)
command = [
    str(repository / "gvclass"),
    "--image",
    str(repository / "gvclass.sif"),
    str(inputs),
    "-o",
    str(output),
    "-t",
    "4",
    "--mode-fast",
    "--plain-output",
    "--verbose",
    "--resource-cache-dir",
    str(temporary / "gvclass-screen-cache"),
]
environment = dict(os.environ)
environment.pop("GVCLASS_DB", None)
before = None
query_before = None
for name, flags, expected_status in (
    ("automatic", [], "ok"),
    ("automatic-resume", ["--resume"], "ok"),
    ("automatic-disabled", ["--resume", "--no-screen"], "disabled"),
):
    # The host DB selects the add-on; the launcher must reset the
    # container DB to its embedded core resources. Resumes cover
    # discovery from the cache with no explicit host DB.
    if not flags:
        environment["GVCLASS_DB"] = str(temporary / "gvclass-screen-cache")
    else:
        environment.pop("GVCLASS_DB", None)
    log_path = evidence / f"{name}.log"
    with log_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [*command, *flags],
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    verify_plain_output(log_path)
    if baseline is None:
        # Compare exact values within this run; thread allocations can
        # change tree results between the one- and three-query runs.
        individual = read_pre_screen_summary(output)
        baseline = {key: individual[key] for key in classification_columns}
        (evidence / "automatic-pre-screen-summary.json").write_text(
            json.dumps(baseline, indent=2), encoding="utf-8"
        )
    verify_summary(
        output / "gvclass_summary.tsv", baseline, standalone, expected_status
    )
    shutil.copy2(output / "gvclass_summary.tsv", evidence / f"{name}-summary.tsv")
    status = json.loads((output / "run_status.json").read_text(encoding="utf-8"))
    assert status["status"] == "complete"
    assert status["totals"]["queries"] == 1
    assert status["totals"]["failed"] == 0
    assert status["totals"]["resume_skipped"] == (0 if not flags else 1)
    snapshot = screen_snapshot(output / "contig_screen")
    query_files = query_snapshot(output)
    if before is None:
        before = snapshot
        query_before = query_files
    else:
        assert snapshot == before, "Resume rewrote or replaced saved screening results"
        assert query_files == query_before, "Resume changed query archive or sidecars"
        log = log_path.read_text(encoding="utf-8")
        assert "Skipping 1 completed queries" in log
        assert "Processing query:" not in log
    (evidence / f"{name}-screen-snapshot.json").write_text(
        json.dumps(snapshot, indent=2), encoding="utf-8"
    )
    (evidence / f"{name}-query-snapshot.json").write_text(
        json.dumps(query_files, indent=2), encoding="utf-8"
    )
    print(f"Verified {name}: {expected_status}, unchanged classification fields")
