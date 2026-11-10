"""Sequence ID reformatting, duplicate IDs, query discovery, and type inference."""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.core.reformat import reformat_sequences
from src.utils import input_validation
from src.utils.error_handling import ValidationError
from src.utils.input_validation import InputValidator, discover_query_files

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
