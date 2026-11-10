"""Splitting multi-contig nucleotide files into one file per contig."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.utils.contig_splitter import (
    ContigIdCollisionError,
    split_contigs,
    split_contigs_from_directory,
)
from src.utils.input_validation import InputValidator


def test_split_contigs_raises_on_colliding_ids(tmp_path: Path) -> None:
    """Two contig IDs that sanitize to the same filename must be detected
    before any files are written."""

    input_fna = tmp_path / "mixed.fna"
    input_fna.write_text(">a/b\nACGTACGTAC\n>a_b\nTTTTTTTTTT\n")
    output_dir = tmp_path / "out"

    with pytest.raises(ContigIdCollisionError, match="collide"):
        split_contigs(input_fna, output_dir=output_dir)

    # The transactional guarantee: no partial writes.
    assert not output_dir.exists() or list(output_dir.iterdir()) == []


def test_split_contigs_passes_through_unique_ids(tmp_path: Path) -> None:
    input_fna = tmp_path / "input.fna"
    input_fna.write_text(">seq1\nACGT\n>seq2\nTTTT\n")
    output_dir = tmp_path / "out"

    returned_dir, count = split_contigs(input_fna, output_dir=output_dir)

    assert returned_dir == output_dir
    assert count == 2
    assert {p.name for p in output_dir.iterdir()} == {"seq1.fna", "seq2.fna"}


def test_split_contigs_from_directory_transactional_across_files(
    tmp_path: Path,
) -> None:
    """Codex audit: a collision in a later file must NOT leave outputs
    from an earlier clean file on disk."""

    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    (input_dir / "good.fna").write_text(">g1\n" + "ACGT" * 10 + "\n")
    # Second file has two contigs that sanitize to the same filename.
    (input_dir / "bad.fna").write_text(">a/b\nACGT\n>a_b\nTTTT\n")

    output_dir = tmp_path / "out"
    with pytest.raises(ContigIdCollisionError):
        split_contigs_from_directory(input_dir, output_dir=output_dir)

    # With the cross-file plan-phase check, nothing should have been
    # written at all. (The output_dir is only created during the write
    # phase, after planning.)
    assert not output_dir.exists() or list(output_dir.iterdir()) == []


def test_fa_and_fas_extension_handling(tmp_path: Path) -> None:
    """.fa files are classified by content and .fas files are split."""
    fa_file = tmp_path / "genome.fa"
    fa_file.write_text(">contig1\n" + "ATGC" * 6000 + "\n")
    validated = InputValidator.validate_sequence_file(fa_file)
    assert validated == fa_file.resolve()

    fas_dir = tmp_path / "fas_in"
    fas_dir.mkdir()
    (fas_dir / "multi.fas").write_text(
        ">c1\n" + "A" * 100 + "\n>c2\n" + "G" * 100 + "\n"
    )
    split_out = tmp_path / "fas_out"
    out_dir, n_contigs, n_files = split_contigs_from_directory(
        fas_dir, split_out, min_length=50
    )
    assert out_dir == split_out
    assert n_contigs == 2
    assert n_files == 1
