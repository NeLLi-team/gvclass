"""Phase 4.3 regression test: document the BLAST-output redundancy so
future changes that remove the unused ``*.blastpout`` artefact have a
clear contract to point at.

Background
----------
``_run_blast_search`` in ``src/pipeline/query_processing_engine.py`` writes
``blastp_out/<query>.<marker>.blastpout`` files. No production consumer
reads those files today:

* ``src/core/contamination_scoring.py`` reads ``blastp_out/*.m8``.
* ``src/core/marker_processing.py::MarkerProcessor.blast_and_merge`` runs
  BLAST independently to produce ``<marker>.m8`` in the same directory.

This means BLAST is effectively run twice per query on two incompatible
naming conventions. Full unification is planned for v1.5.0; this test
ensures that if someone deletes the ``.blastpout`` emission without
also updating the contamination scorer's glob pattern, the regression
is caught.
"""

from __future__ import annotations

from pathlib import Path
import re


def _read(path: Path) -> str:
    return path.read_text()


REPO_ROOT = Path(__file__).resolve().parents[1]
ENGINE = REPO_ROOT / "src" / "pipeline" / "query_processing_engine.py"
CONTAMINATION = REPO_ROOT / "src" / "core" / "contamination_scoring.py"
MARKER_PROCESSING = REPO_ROOT / "src" / "core" / "marker_processing.py"


def test_engine_emits_blastpout_files_for_now() -> None:
    """Snapshot: the engine still writes ``.blastpout`` until v1.5.0 unifies."""
    assert ".blastpout" in _read(ENGINE)


def test_contamination_reads_m8_not_blastpout() -> None:
    """Contamination scoring globs ``*.m8``, not ``*.blastpout``."""
    text = _read(CONTAMINATION)
    assert "blast_dir.glob(\"*.m8\")" in text
    assert "*.blastpout" not in text


def test_marker_processing_still_produces_m8() -> None:
    """The .m8 artefact that contamination scoring depends on is produced
    by MarkerProcessor, not by _run_blast_search."""
    text = _read(MARKER_PROCESSING)
    assert re.search(r"\.m8", text)


def test_run_blast_search_has_v150_todo() -> None:
    """Ensure the Phase 4.3 TODO comment is present so reviewers remember
    to unify the output naming in v1.5.0."""
    assert "v1.5.0" in _read(ENGINE)
    assert "blastpout" in _read(ENGINE)
