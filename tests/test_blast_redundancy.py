"""Regression test locking in the single-BLAST-pass contract.

Until v1.4.3 the per-query pipeline ran BLAST twice: once inside the
engine (writing ``*.blastpout``) and again inside :class:`MarkerProcessor`
(writing ``*.m8``). Only the ``*.m8`` output was consumed downstream by
contamination scoring and tree building. The redundant pass was removed
in v1.4.3 to halve the per-query BLAST workload.

This test fails if either

* the engine grows a second BLAST call in ``_run_hmm_and_blast`` / a new
  ``_run_blast_search`` helper, or
* MarkerProcessor stops emitting ``<marker>.m8`` files (breaks
  contamination scoring).
"""

from __future__ import annotations

from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[1]
ENGINE = REPO_ROOT / "src" / "pipeline" / "query_processing_engine.py"
CONTAMINATION = REPO_ROOT / "src" / "core" / "contamination_scoring.py"
MARKER_PROCESSING = REPO_ROOT / "src" / "core" / "marker_processing.py"


def _read(path: Path) -> str:
    return path.read_text()


def test_engine_no_longer_emits_blastpout_files() -> None:
    """After v1.4.3, the engine must not write any ``*.blastpout`` files."""
    text = _read(ENGINE)
    # Must not assemble any filename ending in ``.blastpout``.
    assert not re.search(r'\.blastpout"', text), (
        "Engine is assembling a .blastpout filename again — the single-"
        "BLAST-pass contract has regressed."
    )
    # Must not re-introduce the helper that owned the duplicate pass.
    assert not re.search(r"def\s+_run_blast_search\s*\(", text), (
        "Engine re-introduced _run_blast_search — merge its logic into "
        "the canonical MarkerProcessor BLAST pass instead."
    )
    # Must not call run_blastp from the engine layer (that's MarkerProcessor's job).
    assert "run_blastp(" not in text, (
        "Engine called run_blastp directly — leave BLAST to MarkerProcessor."
    )


def test_contamination_reads_m8_exactly() -> None:
    """Contamination scoring iterates ``blast_dir.glob('*.m8')``."""
    text = _read(CONTAMINATION)
    assert 'blast_dir.glob("*.m8")' in text
    assert ".blastpout" not in text


def test_marker_processing_produces_exact_m8_path() -> None:
    """MarkerProcessor is the canonical .m8 producer."""
    text = _read(MARKER_PROCESSING)
    assert re.search(
        r'blast_out\s*=\s*self\.blast_dir\s*/\s*f"\{self\.marker\}\.m8"',
        text,
    ), (
        "MarkerProcessor no longer emits {marker}.m8 — contamination "
        "scoring will break."
    )
