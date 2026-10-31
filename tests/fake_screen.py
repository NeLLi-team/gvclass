"""A stand-in for ``python -m src.bin.gvclass_screen`` for the step and flow tests.

The stand-in follows the command line of the real screen and writes
``contig_screen_summary.tsv`` the way it does, with the SHA-256 of the add-on
manifest it finds in the database. The mode of a
bin comes from the first word of its name: ``ok``, ``half`` (the giant-virus screen
failed), ``dead`` (both failed), ``single`` (one screen could not score the bin),
``short`` (the other one), ``crash`` (exit without a summary), ``noisy`` (5 MB on
stderr, then crash), ``hang`` (start a child process, then sleep), ``orphan`` (start a
child process, then exit with code 3), ``torn`` (a summary with empty digests, exit 1),
``wronginput`` (the digest of other bytes) and ``garbled`` (text in a number column).
Every call writes ``args.json`` in the output folder; ``hang`` and ``orphan`` also
write ``child.pid``.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

#: SHA-256 of the file that ``tests/test_contig_screen_step.py`` gives every query.
INPUT_DIGEST = hashlib.sha256(b">c1\nACGT\n").hexdigest()
FAKE_SCREEN = """
import hashlib
import json
import sys
from pathlib import Path

args = sys.argv[1:]
bin_path = Path(args[0])
folder = Path(args[args.index("-o") + 1])
database = Path(args[args.index("-d") + 1])
folder.mkdir(parents=True, exist_ok=True)
(folder / "args.json").write_text(json.dumps(args), encoding="utf-8")
mode = bin_path.stem.split("_")[0]
manifest = database / "contig_screen" / "manifest.json"
digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
HEADER = [
    "bin", "cellular_status", "giant_virus_status", "input_bp", "cellular_bp",
    "cellular_pct", "giant_virus_bp", "giant_virus_pct", "flagged_bp", "flagged_pct",
    "manifest_sha256", "input_sha256",
]
ROWS = {
    "ok": ("ok", "ok", 1000, 100, "10.00", 200, "20.00", 250, "25.00"),
    "half": ("ok", "failed", 1000, 100, "10.00", "NaN", "NaN", "NaN", "NaN"),
    "dead": ("failed", "failed", 1000, "NaN", "NaN", "NaN", "NaN", "NaN", "NaN"),
    "single": ("ok", "not_scored", 1000, 100, "10.00", "NaN", "NaN", "NaN", "NaN"),
    "short": ("not_scored", "ok", 1000, "NaN", "NaN", 200, "20.00", "NaN", "NaN"),
}
for extra in ("torn", "wronginput", "garbled"):
    ROWS[extra] = ROWS["ok"]
if mode == "crash":
    print("boom: out of memory", file=sys.stderr)
    sys.exit(3)
if mode == "noisy":
    sys.stderr.write("x" * 5_000_000 + "\\nlast line of noise\\n")
    sys.exit(3)
if mode in ("hang", "orphan"):
    import subprocess
    import time

    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
    (folder / "child.pid").write_text(str(child.pid), encoding="utf-8")
    if mode == "orphan":
        sys.exit(3)
    time.sleep(600)
input_digest = hashlib.sha256(bin_path.read_bytes()).hexdigest()
if mode == "torn":
    digest = input_digest = ""
if mode == "wronginput":
    input_digest = "0" * 64
row = [bin_path.stem, *map(str, ROWS[mode]), digest, input_digest]
if mode == "garbled":
    row[5] = "oops"
text = "\\t".join(HEADER) + "\\n" + "\\t".join(row) + "\\n"
(folder / "contig_screen_summary.tsv").write_text(text, encoding="utf-8")
sys.exit(1 if mode == "torn" or "failed" in ROWS[mode][:2] else 0)
"""
SUMMARY_HEADER = (
    "bin\tcellular_status\tgiant_virus_status\tinput_bp\tcellular_bp\tcellular_pct\t"
    "giant_virus_bp\tgiant_virus_pct\tflagged_bp\tflagged_pct\tmanifest_sha256\t"
    "input_sha256\n"
)


def install_fake_screen(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Make the step run the stand-in instead of the real screen."""
    script = tmp_path / "fake_screen.py"
    script.write_text(FAKE_SCREEN, encoding="utf-8")
    monkeypatch.setattr(
        "src.pipeline.contig_screen_step.SCREEN_COMMAND", (sys.executable, str(script))
    )


def seed_summary(
    folder: Path,
    row: str,
    digest: str | None,
    header: str = "",
    input_digest: str = INPUT_DIGEST,
) -> None:
    """Write a finished ``contig_screen_summary.tsv`` as an earlier run left it.

    ``row`` has the fields before the manifest digest; ``digest`` None writes ``row``
    as it is, under ``header``, like a row of an older run without the digests.
    """
    folder.mkdir(parents=True, exist_ok=True)
    if digest is None:
        text = header + row
    else:
        text = SUMMARY_HEADER + row.rstrip("\n") + f"\t{digest}\t{input_digest}\n"
    (folder / "contig_screen_summary.tsv").write_text(text, encoding="utf-8")
