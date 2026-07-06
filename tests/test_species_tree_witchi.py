"""Smoke test pinning the witchi CLI contract used by the species-tree supermatrix.

Section 0 of the ``--species-tree`` feature depends on the pinned ``witchi``
(v0.2.0, github ``stephkoest/witchi@22bd7da``) as a pure-python chi-squared
alignment column pruner. Its ``AlignmentPruner`` class is *not* exported in
0.2.0, so the supermatrix builder (Section 4) drives the documented
``witchi prune`` CLI instead.

This test pins the two facts Section 4 relies on:

* ``witchi prune --file <aln> --format fasta`` exits 0 on a tiny alignment.
* The pruned alignment is written *next to the input file* with the name
  ``<stem>_<algorithm>_s<top_n>_pruned.fasta`` (default algorithm
  ``wasserstein``, default ``top_n`` 1 → ``<stem>_wasserstein_s1_pruned.fasta``).

If those change in a future witchi bump, Section 4's output globbing breaks and
this test is the early warning.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

witchi = pytest.importorskip("witchi", reason="witchi not installed in this env")

_WITCHI_CLI = shutil.which("witchi")

_TINY_ALIGNMENT = """\
>seqA
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA
>seqB
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA
>seqC
MKTAYIAKQRQISFVKSHFSRQLEDRLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA
>seqD
MKQAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVAVKVKA
>seqE
MKTAYIWKQRQISFVKSHFSRQLEERLGLIEVQACILSRVGDGTQDNLSGAEKAVQVKVKA
>seqF
MKTAYIAKQRQISFVKSHFSRWLEERLGLIEVQAPILSRVGDGTPDNLSGAEKAVQVKVKA
"""


def test_witchi_version_is_pinned() -> None:
    """The pinned commit resolves to witchi 0.2.0."""
    assert witchi.__version__ == "0.2.0"


def test_witchi_alignment_pruner_unexported() -> None:
    """Document why Section 4 uses the CLI: the class is not exported in 0.2.0."""
    assert not hasattr(witchi, "AlignmentPruner")


@pytest.mark.skipif(_WITCHI_CLI is None, reason="witchi console script not on PATH")
def test_witchi_prune_cli_contract(tmp_path: Path) -> None:
    """`witchi prune` exits 0 and writes `<stem>_wasserstein_s1_pruned.fasta` next to input."""
    aln = tmp_path / "grp_smoke.fasta"
    aln.write_text(_TINY_ALIGNMENT)

    # Run from a *different* cwd to prove witchi writes next to the input, not cwd.
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    proc = subprocess.run(
        [
            _WITCHI_CLI,
            "prune",
            "--file",
            str(aln),
            "--format",
            "fasta",
            "--num_workers_chisq",
            "1",
            "--num_workers_permute",
            "1",
        ],
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert proc.returncode == 0, f"witchi prune failed: {proc.stderr}"

    expected = aln.parent / "grp_smoke_wasserstein_s1_pruned.fasta"
    assert expected.exists(), (
        "expected pruned output next to input as "
        f"{expected.name}; dir held {sorted(p.name for p in aln.parent.iterdir())}"
    )
    # cwd must stay clean (no pruned artefacts leak into the run dir).
    assert not list(run_dir.glob("*_pruned.fasta"))

    # The pruned alignment is valid FASTA preserving all six taxa.
    headers = [
        line[1:].strip()
        for line in expected.read_text().splitlines()
        if line.startswith(">")
    ]
    assert sorted(headers) == ["seqA", "seqB", "seqC", "seqD", "seqE", "seqF"]


def test_witchi_must_use_console_script_not_module() -> None:
    """`python -m witchi` is NOT supported (no __main__) — Section 4 must use the
    console script resolved via ``shutil.which('witchi')``."""
    proc = subprocess.run(
        [sys.executable, "-m", "witchi", "--help"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    # `python -m witchi` must fail (no __main__). The exact stderr wording is
    # interpreter-coupled, so we only require a nonzero exit plus a usable
    # console script — that is the real contract Section 4 depends on.
    assert proc.returncode != 0
    assert _WITCHI_CLI is not None
