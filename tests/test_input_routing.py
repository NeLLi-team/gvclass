"""H1: route .fasta/.fas inputs by inferred alphabet, not the literal .fna suffix.

A DNA .fasta/.fas file must take the nucleotide (gene-calling) path; a protein
.fasta and a .faa must take the protein path.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.pipeline import query_processing_engine as qpe
from src.pipeline.query_processing_engine import PreparedQueryInput


@pytest.fixture
def routing_recorder(monkeypatch):
    calls: dict[str, Path] = {}

    def fake_nuc(query_file, *args, **kwargs):
        calls["kind"] = "nucleotide"
        return PreparedQueryInput(
            protein_file=query_file,
            reformatted_file=query_file,
            best_code=11,
            is_nucleotide=True,
        )

    def fake_prot(query_file, *args, **kwargs):
        calls["kind"] = "protein"
        return PreparedQueryInput(
            protein_file=query_file,
            reformatted_file=query_file,
            best_code=None,
            is_nucleotide=False,
        )

    monkeypatch.setattr(qpe, "_prepare_nucleotide_query", fake_nuc)
    monkeypatch.setattr(qpe, "_prepare_protein_query", fake_prot)
    return calls


def _route(tmp_path: Path, name: str, content: str, calls: dict):
    qf = tmp_path / name
    qf.write_text(content)
    prepared = qpe._prepare_query_input(qf, tmp_path, tmp_path, [0, 11], 1, None)
    return calls["kind"], prepared


@pytest.mark.parametrize("ext", [".fasta", ".fas"])
def test_dna_fasta_like_routes_to_nucleotide(tmp_path, routing_recorder, ext):
    kind, prepared = _route(
        tmp_path, f"q{ext}", ">c1\nATGAAATAG\n", routing_recorder
    )
    assert kind == "nucleotide"
    assert prepared.is_nucleotide is True


@pytest.mark.parametrize("ext", [".fasta", ".fas"])
def test_protein_fasta_like_routes_to_protein(tmp_path, routing_recorder, ext):
    kind, prepared = _route(tmp_path, f"q{ext}", ">p1\nMKVLAA\n", routing_recorder)
    assert kind == "protein"
    assert prepared.is_nucleotide is False


def test_fna_routes_to_nucleotide(tmp_path, routing_recorder):
    kind, _ = _route(tmp_path, "q.fna", ">c1\nATGAAATAG\n", routing_recorder)
    assert kind == "nucleotide"


def test_faa_routes_to_protein(tmp_path, routing_recorder):
    kind, _ = _route(tmp_path, "q.faa", ">p1\nMKVLAA\n", routing_recorder)
    assert kind == "protein"


def test_copy_gate_writes_for_nucleotide_only(tmp_path):
    qout = tmp_path / "qout"
    for sub in ("query_faa", "query_fna", "query_gff"):
        (qout / sub).mkdir(parents=True)
    protein = tmp_path / "p.faa"
    protein.write_text(">p\nM\n")
    reformatted = tmp_path / "r.fna"
    reformatted.write_text(">c\nATG\n")

    # Protein input (is_nucleotide=False) -> nothing copied.
    qpe._copy_final_sequence_outputs(
        Path("q.fasta"), "q", qout, reformatted, protein, None, False, _NullLogger()
    )
    assert not (qout / "query_fna" / "q.fna").exists()

    # Nucleotide input (is_nucleotide=True) -> fna/faa copied.
    qpe._copy_final_sequence_outputs(
        Path("q.fasta"), "q", qout, reformatted, protein, 11, True, _NullLogger()
    )
    assert (qout / "query_fna" / "q.fna").exists()
    assert (qout / "query_faa" / "q.faa").exists()


class _NullLogger:
    def info(self, *a, **k):
        pass

    def warning(self, *a, **k):
        pass
