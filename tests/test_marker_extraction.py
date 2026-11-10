"""Marker sequence extraction from a query protein file."""

from __future__ import annotations

from pathlib import Path

from src.core.marker_extraction import extract_marker_sequences


def test_extract_marker_sequences_matches_by_record_id(tmp_path: Path) -> None:
    query_faa = tmp_path / "q.faa"
    query_faa.write_text(
        ">prot_1 # start=1 # end=300 # annotation\n"
        "MKLILV\n"
        ">prot_2 some descriptor\n"
        "MTKAAG\n"
    )
    output_dir = tmp_path / "markers"

    marker_hits = {"GVOGm0022": {"prot_1"}, "GVOGm0054": {"prot_2"}}

    files = extract_marker_sequences(query_faa, marker_hits, output_dir)

    assert set(files.keys()) == {"GVOGm0022", "GVOGm0054"}
    for marker in files:
        contents = files[marker].read_text()
        assert contents.startswith(">prot_"), contents
