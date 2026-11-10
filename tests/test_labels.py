"""The label index: one loader, one genome-ID resolver, one nearest-lineage rule."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.core.labels import LabelIndex, domain_prefix

LABEL_LINES = (
    "NCLDV__Mimi_1\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|"
    "Mimiviridae|Mimivirus|Mimi_one\n"
    "NCLDV__Mimi\tNCLDV|Nucleocytoviricota|Megaviricetes|Pimascovirales|"
    "Marseilleviridae|Marseillevirus|Mimi\n"
    "PHAGE__GCA-1\tPHAGE|Uroviricota|Caudoviricetes|Crassvirales|f|g|s\n"
    "BAC__ref|contig_7\tBAC|Pseudomonadota|Gamma|Enterobacterales|Enterobacteriaceae"
    "|Escherichia|coli\n"
    "BAC__short\tBAC|Pseudomonadota\n"
)


def _write_labels(tmp_path: Path, text: str = LABEL_LINES) -> Path:
    labels_file = tmp_path / "labels.tsv"
    labels_file.write_text(text, encoding="utf-8")
    return labels_file


def test_load_builds_a_seven_rank_lineage_per_genome(tmp_path: Path) -> None:
    index = LabelIndex.load(_write_labels(tmp_path))

    assert index.lineages["NCLDV__Mimi"] == (
        "NCLDV",
        "Nucleocytoviricota",
        "Megaviricetes",
        "Pimascovirales",
        "Marseilleviridae",
        "Marseillevirus",
        "Mimi",
    )
    assert index.lineages["BAC__short"] == (
        "BAC",
        "Pseudomonadota",
        "",
        "",
        "",
        "",
        "",
    )


def test_domain_comes_from_the_id_prefix_before_double_underscore(
    tmp_path: Path,
) -> None:
    labels_file = _write_labels(tmp_path, "EUK-pEVE__x\tEUK|p|c|o|f|g|s\n")

    assert LabelIndex.load(labels_file).lineages["EUK-pEVE__x"][0] == "EUK-pEVE"
    assert domain_prefix("EUK-pEVE__x|contig_1") == "EUK-pEVE"
    assert domain_prefix("plain_id") == ""


def test_load_skips_comments_and_blank_lines(tmp_path: Path) -> None:
    labels_file = _write_labels(tmp_path, "# header\n\nNCLDV__a\tNCLDV|p|c|o|f|g|s\n")

    assert set(LabelIndex.load(labels_file).lineages) == {"NCLDV__a"}


@pytest.mark.parametrize(
    "text",
    ["malformed-label\n", "NCLDV__a\t\n", "\tNCLDV|p\n", "a\tb\tc\n"],
)
def test_load_rejects_a_malformed_record(tmp_path: Path, text: str) -> None:
    with pytest.raises(ValueError, match=r"labels\.tsv at line 1"):
        LabelIndex.load(_write_labels(tmp_path, text))


def test_load_rejects_an_empty_file(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no records"):
        LabelIndex.load(_write_labels(tmp_path, "# only a comment\n"))


def test_load_rejects_a_missing_file(tmp_path: Path) -> None:
    with pytest.raises(OSError, match=r"labels\.tsv"):
        LabelIndex.load(tmp_path / "labels.tsv")


def test_load_reads_each_file_once(tmp_path: Path) -> None:
    labels_file = _write_labels(tmp_path)

    assert LabelIndex.load(labels_file) is LabelIndex.load(labels_file)


def test_load_reads_a_changed_file_again(tmp_path: Path) -> None:
    labels_file = _write_labels(tmp_path, "NCLDV__a\tNCLDV|p|c|o|f|g|s\n")
    first = LabelIndex.load(labels_file)

    labels_file.write_text(
        "NCLDV__a\tNCLDV|p|c|o|f|g|s\nNCLDV__b\tNCLDV|p|c|o|f|g|s\n",
        encoding="utf-8",
    )

    assert set(LabelIndex.load(labels_file).lineages) == {"NCLDV__a", "NCLDV__b"}
    assert set(first.lineages) == {"NCLDV__a"}


PKV_CANONICAL_LINE = (
    "NCLDV__PkV-RF01\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|IM_19|"
    "g1787|singleton\n"
)
PKV_DUPLICATE_LINE = (
    "NCLDV__GCF_905367645.1\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|"
    "IM_19|NCLDV_unclassified|PkV-RF01\n"
)
PKV_LINEAGE = (
    "NCLDV",
    "Nucleocytoviricota",
    "Megaviricetes",
    "Imitervirales",
    "IM_19",
    "g1787",
    "singleton",
)


def test_load_gives_the_duplicate_pkv_record_the_canonical_lineage(
    tmp_path: Path,
) -> None:
    labels_file = _write_labels(tmp_path, PKV_CANONICAL_LINE + PKV_DUPLICATE_LINE)

    lineages = LabelIndex.load(labels_file).lineages

    assert lineages["NCLDV__PkV-RF01"] == PKV_LINEAGE
    assert lineages["NCLDV__GCF_905367645.1"] == PKV_LINEAGE


def test_load_keeps_the_duplicate_pkv_record_without_the_canonical_record(
    tmp_path: Path,
) -> None:
    labels_file = _write_labels(tmp_path, PKV_DUPLICATE_LINE)

    lineages = LabelIndex.load(labels_file).lineages

    assert lineages["NCLDV__GCF_905367645.1"][5:] == ("NCLDV_unclassified", "PkV-RF01")


def test_load_leaves_other_records_unchanged_beside_the_pkv_pair(
    tmp_path: Path,
) -> None:
    labels_file = _write_labels(
        tmp_path, LABEL_LINES + PKV_CANONICAL_LINE + PKV_DUPLICATE_LINE
    )

    lineages = LabelIndex.load(labels_file).lineages

    assert lineages["NCLDV__Mimi"][3:] == (
        "Pimascovirales",
        "Marseilleviridae",
        "Marseillevirus",
        "Mimi",
    )
    assert lineages["NCLDV__Mimi_1"][5:] == ("Mimivirus", "Mimi_one")
    assert set(lineages) == {
        "NCLDV__Mimi_1",
        "NCLDV__Mimi",
        "PHAGE__GCA-1",
        "BAC__ref|contig_7",
        "BAC__short",
        "NCLDV__PkV-RF01",
        "NCLDV__GCF_905367645.1",
    }


def test_index_is_immutable(tmp_path: Path) -> None:
    index = LabelIndex.load(_write_labels(tmp_path))

    with pytest.raises(AttributeError):
        index.lineages = {}  # type: ignore[misc]  # frozen dataclass
    with pytest.raises(TypeError):
        index.lineages["NCLDV__new"] = ("NCLDV",) * 7  # type: ignore[index]  # read-only view


@pytest.mark.parametrize(
    ("protein_id", "genome_id"),
    [
        # 1. an exact label key wins, even when it ends in _<digits>
        ("NCLDV__Mimi_1", "NCLDV__Mimi_1"),
        ("BAC__ref|contig_7", "BAC__ref|contig_7"),
        # 2. the ID without the trailing gene number
        ("NCLDV__Mimi_7", "NCLDV__Mimi"),
        ("PHAGE__GCA-1_7", "PHAGE__GCA-1"),
        ("BAC__ref|contig_7_42", "BAC__ref|contig_7"),
        # 3. the part before the first |
        ("BAC__short|ctg_3", "BAC__short"),
        ("NCLDV__Mimi|anything", "NCLDV__Mimi"),
        # 4. otherwise the stripped ID
        ("UNKNOWN__x_12", "UNKNOWN__x"),
        ("UNKNOWN__x|y_3", "UNKNOWN__x|y"),
    ],
)
def test_resolve_genome_id_precedence(
    tmp_path: Path, protein_id: str, genome_id: str
) -> None:
    index = LabelIndex.load(_write_labels(tmp_path))

    assert index.resolve_genome_id(protein_id) == genome_id


LINEAGE_TAXONOMY = {
    "NCLDV__a": ["NCLDV", "P", "C", "O1", "F1", "G1", "S1"],
    "NCLDV__b": ["NCLDV", "P", "C", "O2", "F2", "G2", "S2"],
    "NCLDV__twin": ["NCLDV", "P", "C", "O1", "F1", "G1", "S1"],
    "NCLDV__sibling": ["NCLDV", "P", "C", "O1", "F1", "G9", "S9"],
    "BAC__x": ["BAC", "Pseudomonadota", "Gamma", "Entero", "Ef", "Eg", "Es"],
}


def _nearest(neighbors: dict[str, float]) -> tuple[tuple[str, ...], float] | None:
    return LabelIndex.from_taxonomy(LINEAGE_TAXONOMY).nearest_lineage("m1", neighbors)


def test_nearest_lineage_of_one_neighbour_is_its_lineage() -> None:
    assert _nearest({"NCLDV__a|p_1": 0.2}) == (
        ("NCLDV", "P", "C", "O1", "F1", "G1", "S1"),
        0.2,
    )


def test_nearest_lineage_of_a_tie_keeps_the_ranks_the_neighbours_share() -> None:
    lineage, distance = _nearest({"NCLDV__a|p_1": 0.2, "NCLDV__sibling|p_1": 0.2})

    assert lineage == ("NCLDV", "P", "C", "O1", "F1", "", "")
    assert distance == 0.2


def test_nearest_lineage_of_a_tie_across_orders_stops_at_the_class() -> None:
    lineage, _ = _nearest({"NCLDV__a|p_1": 0.2, "NCLDV__b|p_1": 0.2})

    assert lineage == ("NCLDV", "P", "C", "", "", "", "")


def test_nearest_lineage_of_a_tie_across_domains_has_no_domain() -> None:
    lineage, _ = _nearest({"NCLDV__a|p_1": 0.2, "BAC__x|p_1": 0.2})

    assert lineage == ("",) * 7


def test_nearest_lineage_ignores_farther_and_unlabelled_neighbours() -> None:
    lineage, distance = _nearest(
        {"UNKNOWN__z|p_1": 0.05, "NCLDV__a|p_1": 0.2, "NCLDV__b|p_1": 0.3}
    )

    assert lineage[3] == "O1"
    assert distance == 0.2


def test_nearest_lineage_is_none_without_a_labelled_neighbour() -> None:
    assert _nearest({"UNKNOWN__z|p_1": 0.05}) is None
    assert _nearest({}) is None
