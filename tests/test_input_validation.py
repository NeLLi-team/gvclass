"""Sequence ID reformatting, duplicate IDs, query discovery, and type inference."""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.core.reformat import calculate_stats, reformat_sequences
from src.pipeline.query_processing_engine import is_nucleotide_input
from src.utils import input_validation
from src.utils.common import validate_file_path
from src.utils.error_handling import ValidationError
from src.utils.input_validation import (
    InputValidator,
    discover_query_files,
    is_safe_query_stem,
)

LONG_DNA = "ACGT" * 6000


def _records(*ids: str) -> list[SeqRecord]:
    return [SeqRecord(Seq("ACGT"), id=sequence_id) for sequence_id in ids]


def _new_ids(*ids: str, input_file: str = "q.fna") -> list[str]:
    return [record.id for record in reformat_sequences(_records(*ids), input_file)]


def test_reformat_keeps_contigs_of_different_bins_distinct() -> None:
    assert _new_ids("binA|contig_1", "binB|contig_1") == [
        "q|binA__contig_1",
        "q|binB__contig_1",
    ]


def test_reformat_keeps_the_whole_ncbi_style_id() -> None:
    assert _new_ids("gi|9626243|ref|NC_001416.1|") == [
        "q|gi__9626243__ref__NC_001416.1"
    ]


def test_reformat_is_stable_when_its_own_output_is_fed_back() -> None:
    assert _new_ids("q|binA__contig_1", "q|ctg") == ["q|binA__contig_1", "q|ctg"]


def test_reformat_removes_only_one_leading_basename_prefix() -> None:
    assert _new_ids("q|q|ctg") == ["q|q__ctg"]


def test_reformat_replaces_colons_and_keeps_inner_underscores() -> None:
    assert _new_ids("_a:b_", "x||y") == ["q|_a__b_", "q|x____y"]


def test_reformat_clears_descriptions() -> None:
    record = SeqRecord(Seq("ACGT"), id="ctg", description="ctg some text")

    reformat_sequences([record], "q.fna")

    assert record.description == ""


@pytest.mark.parametrize("empty_id", ["q|", "||", "|"])
def test_reformat_rejects_an_id_with_no_content(empty_id: str) -> None:
    with pytest.raises(ValidationError, match="no content"):
        reformat_sequences(_records(empty_id), "q.fna")


def test_reformat_rejects_ids_that_become_equal_and_names_both() -> None:
    with pytest.raises(ValidationError) as excinfo:
        reformat_sequences(_records("ok", "x|a", "x__a"), "q.fna")

    message = str(excinfo.value)
    assert "'x|a'" in message
    assert "'x__a'" in message
    assert "q|x__a" in message


def test_reformat_rejects_the_same_original_id_twice() -> None:
    with pytest.raises(ValidationError, match="'ctg' and 'ctg'"):
        reformat_sequences(_records("ctg", "ctg"), "q.fna")


def test_validate_sequence_file_rejects_duplicate_record_ids(tmp_path: Path) -> None:
    path = tmp_path / "dup.fna"
    path.write_text(f">ctg1\n{LONG_DNA}\n>ctg2\n{LONG_DNA}\n>ctg1\n{LONG_DNA}\n")

    with pytest.raises(ValidationError, match="ctg1"):
        InputValidator.validate_sequence_file(path)


def test_discover_query_files_finds_any_case_sorted_by_name(tmp_path: Path) -> None:
    for name in ("b.fna", "A.FNA", "c.Fasta", "d.faa", "notes.txt", "e.fna.bak"):
        (tmp_path / name).write_text(">c\nACGT\n")
    (tmp_path / "sub.fna").mkdir()

    found = discover_query_files(tmp_path)

    assert [path.name for path in found] == ["A.FNA", "b.fna", "c.Fasta", "d.faa"]


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [("ACGTACGTACGTNNACGT", "fna"), ("MKWVNDHRYEFLPQ" * 3, "faa")],
)
def test_infer_fasta_type_reads_the_sequence_alphabet(
    tmp_path: Path, sequence: str, expected: str
) -> None:
    path = tmp_path / "x.fasta"
    path.write_text(f">first\n{sequence}\n")

    assert InputValidator.infer_fasta_type(path) == expected


def test_infer_fasta_type_skips_empty_records(tmp_path: Path) -> None:
    path = tmp_path / "x.fasta"
    path.write_text(">empty\n>real\nMKWVNDHRYEFLPQ\n")

    assert InputValidator.infer_fasta_type(path) == "faa"


def test_infer_fasta_type_stops_at_the_first_informative_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "x.fasta"
    path.write_text(">a\nMKWVNDHRYEFLPQ\n>b\nACGT\n")
    pulled: list[str] = []

    def lazy_parse(handle: object, fmt: str):
        for name in ("a", "b", "c"):
            pulled.append(name)
            yield SeqRecord(Seq("MKWVNDHRYEFLPQ"), id=name)

    monkeypatch.setattr(input_validation.SeqIO, "parse", lazy_parse)

    assert InputValidator.infer_fasta_type(path) == "faa"
    assert pulled == ["a"]


@pytest.mark.parametrize("stem", ["sample", "a.b", "bin-1_x", "sample.previous.abc123"])
def test_is_safe_query_stem_accepts_ordinary_names(stem: str) -> None:
    assert is_safe_query_stem(stem) is True


@pytest.mark.parametrize("stem", ["", ".", "..", ".hidden", "...", "a/b", "../x"])
def test_is_safe_query_stem_rejects_names_that_leave_the_work_directory(
    stem: str,
) -> None:
    assert is_safe_query_stem(stem) is False


def test_short_fna_rejected_by_default(tmp_path: Path) -> None:
    """Phase 3.2a: ``validate_sequence_file`` must raise on sub-20 kb
    aggregate nucleotide input unless ``allow_short=True``."""

    short_fna = tmp_path / "short.fna"
    short_fna.write_text(">c1\n" + "ACGT" * 100 + "\n")  # 400 bp << 20 kb

    with pytest.raises(ValidationError, match="below the minimum"):
        InputValidator.validate_sequence_file(short_fna)


def test_short_fna_accepted_when_allow_short_enabled(tmp_path: Path) -> None:
    short_fna = tmp_path / "short.fna"
    short_fna.write_text(">c1\n" + "ACGT" * 100 + "\n")
    result = InputValidator.validate_sequence_file(short_fna, allow_short=True)
    assert result == short_fna.resolve()


def test_short_fna_accepted_with_custom_min_sequence_length(tmp_path: Path) -> None:
    short_fna = tmp_path / "short.fna"
    short_fna.write_text(">c1\n" + "ACGT" * 100 + "\n")  # 400 bp

    result = InputValidator.validate_sequence_file(
        short_fna,
        min_sequence_length=400,
    )

    assert result == short_fna.resolve()


def test_short_fna_accepted_with_zero_min_sequence_length(tmp_path: Path) -> None:
    short_fna = tmp_path / "short.fna"
    short_fna.write_text(">c1\nACGT\n")

    result = InputValidator.validate_sequence_file(
        short_fna,
        min_sequence_length=0,
    )

    assert result == short_fna.resolve()


def test_short_fna_rejected_against_custom_min_sequence_length(tmp_path: Path) -> None:
    short_fna = tmp_path / "short.fna"
    short_fna.write_text(">c1\n" + "ACGT" * 100 + "\n")  # 400 bp

    with pytest.raises(ValidationError, match="minimum \\(500 bp\\)"):
        InputValidator.validate_sequence_file(short_fna, min_sequence_length=500)


def test_fasta_extension_content_is_validated(tmp_path: Path) -> None:
    """Phase 3.2c: .fasta files must go through DNA/protein alphabet checks."""

    bad_fasta = tmp_path / "bad.fasta"
    # Invalid characters that previously slipped through because .fasta
    # did not route to either alphabet check. The content inference now
    # classifies this as protein (mixed alphabet) and the protein-chars
    # check rejects the digits — either message is acceptable, we just
    # need the rejection to happen at all.
    bad_fasta.write_text(">c1\nATGC123\n")

    with pytest.raises(ValidationError, match="invalid (?:DNA|protein) characters"):
        InputValidator.validate_sequence_file(bad_fasta, allow_short=True)


def test_pure_dna_fasta_routes_to_dna_alphabet(tmp_path: Path) -> None:
    """A clean DNA-only .fasta must resolve to the ``fna`` branch (no
    spurious protein-alphabet rejection)."""

    good_fasta = tmp_path / "good.fasta"
    good_fasta.write_text(">c1\n" + "ACGT" * 10 + "\n")
    # No exception -> routed correctly.
    result = InputValidator.validate_sequence_file(good_fasta, allow_short=True)
    assert result == good_fasta.resolve()


def test_protein_fasta_of_iupac_letters_is_not_called_dna(tmp_path: Path) -> None:
    """Every letter of ``MKWVNDHRYSGACD`` is also an IUPAC ambiguity code, so a
    set-membership test against ``VALID_DNA_CHARS`` calls this peptide DNA."""

    peptides = tmp_path / "peptides.fasta"
    peptides.write_text(">p1\nMKWVNDHRYSGACD\n>p2\nMKWVNDHRYSGACD\n")

    records = list(SeqIO.parse(str(peptides), "fasta"))
    assert InputValidator._infer_fasta_type_from_content(records) == "faa"

    # Protein path: no DNA-alphabet rejection and no 20 kb nucleotide floor.
    assert InputValidator.validate_sequence_file(peptides) == peptides.resolve()


def test_ambiguity_rich_contig_is_still_called_dna(tmp_path: Path) -> None:
    """The other side of the same discriminator: GVClass accepts IUPAC
    ambiguity codes in nucleotide input, and ``is_nucleotide_input`` routes
    gene calling off this helper, so a misread here sends a real contig to the
    protein branch."""

    contig = tmp_path / "ambiguous.fasta"
    contig.write_text(">c1\nATGNNNRYKMSWATGCCGT\n")

    records = list(SeqIO.parse(str(contig), "fasta"))
    assert InputValidator._infer_fasta_type_from_content(records) == "fna"
    assert is_nucleotide_input(contig) is True


def test_low_complexity_protein_repeats_are_not_called_dna(tmp_path: Path) -> None:
    """Gly/Ala-rich repeats (collagen, elastin, silk fibroin) use only letters
    that are also nucleotide codes, ``S`` among them, so neither an alphabet
    test nor a base-composition test alone separates them from DNA. Base
    diversity does: real nucleotide sequence uses all four bases."""

    class _Record:
        def __init__(self, seq: str) -> None:
            self.seq = seq

    for motif in ("GAGAGSGAGAGSGAGAGS", "GAPGAPGAPGAPGAPGAP", "GGYGGGGRGGYGGGGRGG"):
        assert (
            InputValidator._infer_fasta_type_from_content([_Record(motif)]) == "faa"
        ), motif

    fibroin = tmp_path / "fibroin.fasta"
    fibroin.write_text(">p1\nGAGAGSGAGAGSGAGAGS\n")
    assert is_nucleotide_input(fibroin) is False

    # An assembly-gap N-run must not tip a contig onto the protein path either.
    gapped = tmp_path / "gapped.fasta"
    gapped.write_text(">c1\nATGAAA" + "N" * 40 + "CCGGTTACGT\n")
    assert is_nucleotide_input(gapped) is True


def test_file_size_cap_admits_assemblies_over_the_old_50mb_limit() -> None:
    """The cap was 50 MB under a constant named for sequence length in bp, but
    compared against the file size in bytes, so ``--contigs`` on a metagenome
    assembly failed with "File too large". Writing a >50 MB FASTA would slow the
    suite, so assert the cap itself: this fails if it regresses to 50 MB."""

    assert InputValidator.MAX_FILE_SIZE_BYTES > 50_000_000


def test_query_enumeration_accepts_fasta_extensions(tmp_path: Path) -> None:
    """Codex audit: a directory containing only ``.fasta`` / ``.fas``
    inputs must be enumerated as queries, not silently produce zero.
    """

    query_dir = tmp_path / "queries"
    query_dir.mkdir()
    (query_dir / "a.fasta").write_text(">c\nACGT\n")
    (query_dir / "b.fas").write_text(">c\nACGT\n")
    assert [path.name for path in discover_query_files(query_dir)] == [
        "a.fasta",
        "b.fas",
    ]


def test_validate_path_safety_allows_double_dots_in_filename(tmp_path: Path) -> None:
    """A '..' inside the file name is safe when no path component is '..'."""
    f = tmp_path / "genome..v1.fna"
    f.write_text(">c1\nATGC\n")
    resolved = validate_file_path(f)
    assert resolved == f.resolve()


def test_reformat_calculate_stats_counts_lowercase_gc() -> None:
    """calculate_stats should count both uppercase and lowercase G/C bases."""

    records = [SeqRecord(Seq("atgc"), id="c1")]
    stats_df = calculate_stats(records, "lower.fna")
    assert float(stats_df.loc[0, "GCperc"]) == 50.0
