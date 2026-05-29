"""M10: enforce the pyswrd top-100-per-query cap.

pyswrd's default max_alignments is 10, silently capping the documented
"top 100 per query" before run_blastp's own slice. run_blastp must pass
max_alignments explicitly, and parse_blastp's per-query cap must be honored.
"""

from __future__ import annotations

from pathlib import Path

from src.core import blast
from src.core.blast import parse_blastp, run_blastp


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


def test_parse_blastp_per_query_cap_can_exceed_ten(tmp_path: Path) -> None:
    lines = [
        f"q1\tt{i}\t90.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t{200 - i}.0"
        for i in range(15)
    ]
    out = tmp_path / "b.m8"
    out.write_text("\n".join(lines) + "\n")

    assert len(parse_blastp(str(out), max_hits_per_query=100)) == 15
    assert len(parse_blastp(str(out), max_hits_per_query=10)) == 10
