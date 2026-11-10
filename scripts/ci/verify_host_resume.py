"""Exercise two host launcher resume routes against one completed SIF run."""

import hashlib
import json
import subprocess
import sys
import tarfile
import time
from pathlib import Path


def scientific_snapshot(output: Path) -> dict[str, str]:
    """Hash summary bytes and archive file payloads, excluding gzip metadata."""
    hashes = {}
    for path in sorted(output.iterdir()):
        if path.suffix in {".tsv", ".csv"}:
            with path.open("rb") as handle:
                hashes[path.name] = hashlib.file_digest(handle, "sha256").hexdigest()
        if path.name.endswith(".tar.gz"):
            with tarfile.open(path) as archive:
                for member in archive.getmembers():
                    if member.isfile():
                        handle = archive.extractfile(member)
                        assert handle is not None, member.name
                        with handle:
                            hashes[f"{path.name}::{member.name}"] = hashlib.file_digest(
                                handle, "sha256"
                            ).hexdigest()
    return hashes


def results_digest(output: Path) -> str:
    """Hash ``results.jsonl``, the record of finished queries."""
    with (output / "results.jsonl").open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify_resume(
    output: Path,
    log: str,
    expected: dict[str, str],
    times: dict[str, str],
    baseline_json: str,
    previous_run_id: str,
) -> None:
    """Require three skips, no processing diagnostics and unchanged scientific data."""
    status = json.loads((output / "run_status.json").read_text(encoding="utf-8"))
    baseline = json.loads(baseline_json)
    assert status["status"] == "complete"
    assert status["run_id"] != previous_run_id, "Resume retained the previous run ID"
    assert status["database"] == baseline["database"]
    assert status["software_version"] == baseline["software_version"]
    assert status["settings"]["resume"] is True
    ignored = {"resume", "max_workers", "threads_per_worker"}
    assert {
        key: value for key, value in status["settings"].items() if key not in ignored
    } == {
        key: value for key, value in baseline["settings"].items() if key not in ignored
    }, "Resume changed stable settings"
    assert status["totals"]["resume_skipped"] == 3
    assert status["totals"]["queries"] == 3
    assert status["totals"]["failed"] == 0
    assert "Skipping 3 completed queries" in log
    assert "Processing 0 remaining queries" in log
    assert not any(
        marker in log
        for marker in (
            "Processing query:",
            "Running pyswrd search:",
            "Running HMM search:",
            "Optimizing genetic code for:",
            " - DEBUG - Processing marker ",
        )
    ), "Resume started scientific processing"
    assert scientific_snapshot(output) == expected, "Resume changed scientific output"
    observed_times = {
        name: query["finished_at"] for name, query in status["queries"].items()
    }
    assert observed_times == times, "Resume changed query finish times"
    assert all(query["resume_skipped"] for query in status["queries"].values())


image, output, cache, evidence = map(Path, sys.argv[1:])
status = json.loads((output / "run_status.json").read_text(encoding="utf-8"))
assert len(status["queries"]) == 3
assert (output / "gvclass_summary.tsv").is_file()
assert (output / "gvclass_summary.csv").is_file()
assert all((output / f"{name}.tar.gz").is_file() for name in status["queries"])
baseline_json = json.dumps(status)
previous_run_id = status["run_id"]
before = scientific_snapshot(output)
times = {name: query["finished_at"] for name, query in status["queries"].items()}
results_before = results_digest(output)
(evidence / "resume-before.json").write_text(
    json.dumps(
        {"scientific_sha256": before, "query_times": times, "run_status": status},
        indent=2,
    ),
    encoding="utf-8",
)
launcher = Path("gvclass").resolve()
routes = [
    ("launcher-explicit", ["--apptainer", "--image", str(image)]),
    ("launcher-implicit", [f"--image={image}"]),
]
for name, options in routes:
    command = [
        sys.executable,
        str(launcher),
        *options,
        "example",
        "-o",
        str(output),
        "--resource-cache-dir",
        str(cache),
        "-t",
        "4",
        "--mode-fast",
        "--plain-output",
        "--verbose",
        "--resume",
    ]
    log_path = evidence / f"resume-{name}.log"
    time.sleep(1.1)  # Run IDs have one-second precision.
    with log_path.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            command, stdout=handle, stderr=subprocess.STDOUT, check=False
        )
    observed_status = json.loads(
        (output / "run_status.json").read_text(encoding="utf-8")
    )
    (evidence / f"resume-{name}-observed.json").write_text(
        json.dumps(
            {
                "returncode": completed.returncode,
                "run_status": observed_status,
                "scientific_sha256": scientific_snapshot(output),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    completed.check_returncode()
    verify_resume(
        output,
        log_path.read_text(encoding="utf-8"),
        before,
        times,
        baseline_json,
        previous_run_id,
    )
    assert results_digest(output) == results_before, "Resume changed results.jsonl"
    previous_run_id = observed_status["run_id"]
    (evidence / f"resume-{name}.json").write_text(
        json.dumps(
            {
                "route": name,
                "skipped_queries": 3,
                "reprocessed_queries": 0,
                "scientific_outputs_unchanged": True,
                "query_times_unchanged": True,
            }
        ),
        encoding="utf-8",
    )
    print(f"Verified {name}: three skips, no recomputation, unchanged outputs")
