"""M10: enforce the pyswrd top-100-per-query cap.

pyswrd's default max_alignments is 10, silently capping the documented
"top 100 per query" before run_blastp's own slice. run_blastp must pass
max_alignments explicitly.
"""

from __future__ import annotations

from pathlib import Path

from src.core import blast
from src.core.blast import run_blastp


def test_run_blastp_passes_max_alignments(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "q.faa"
    query.write_text(">q1\nMKVLAA\n")
    ref = tmp_path / "r.faa"
    ref.write_text(">r1\nMKVLAA\n")
    out = tmp_path / "out.m8"

    captured: dict = {}

    def fake_search(queries, targets, **kwargs):
        captured.update(kwargs)
        return iter(())

    monkeypatch.setattr(blast.pyswrd, "search", fake_search)
    run_blastp(str(query), str(ref), str(out))

    # Default top_per_query is 100; pyswrd must be told to return that many.
    assert captured.get("max_alignments") == 100

