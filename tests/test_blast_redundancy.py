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
    """Snapshot: the engine still writes ``.blastpout`` files keyed by
    ``{query_name}.{marker}`` until v1.5.0 unifies. Lock in the exact
    emission pattern so a silent rename does not skate past the test."""
    text = _read(ENGINE)
    assert re.search(
        r'blast_out\s*=\s*blast_dir\s*/\s*f"\{query_name\}\.\{marker\}\.blastpout"',
        text,
    ), "engine no longer writes {query}.{marker}.blastpout — unify BLAST naming"


def test_contamination_reads_m8_exactly() -> None:
    """Contamination scoring must iterate ``blast_dir.glob("*.m8")`` at
    the specific site in collect_contig_features. Pinning the literal
    expression ensures a glob-pattern regression is caught."""
    text = _read(CONTAMINATION)
    assert 'blast_dir.glob("*.m8")' in text
    # And must not accidentally pick up the stale blastpout naming.
    assert "*.blastpout" not in text


def test_marker_processing_produces_exact_m8_path() -> None:
    """MarkerProcessor is the canonical .m8 producer. Lock in the exact
    path pattern so any refactor that drops the .m8 suffix trips here
    rather than silently breaks contamination scoring."""
    text = _read(MARKER_PROCESSING)
    assert re.search(
        r'blast_out\s*=\s*self\.blast_dir\s*/\s*f"\{self\.marker\}\.m8"',
        text,
    ), "MarkerProcessor no longer emits {marker}.m8 — contamination scoring will break"


def test_run_blast_search_has_scoped_v150_todo() -> None:
    """The Phase 4.3 TODO must stay attached to the _run_blast_search
    function body so reviewers see it at the call site. Match the exact
    comment token rather than a bare substring."""
    text = _read(ENGINE)
    # Find the function def and ensure the TODO appears within the next
    # 30 lines of the signature (scoped assertion rather than a
    # whole-file substring match).
    match = re.search(r"def _run_blast_search\(", text)
    assert match is not None
    window = text[match.start() : match.start() + 1500]
    assert "TODO(v1.5.0)" in window, "v1.5.0 TODO not attached to _run_blast_search"
    assert "blastpout" in window
