"""Tests of the best-hit selection and the per-window reference features."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.contig_screen import ref_hits

N_REFS = 40
CLASS_OF_REF = np.array(["BAC", "ARC", "EUK"] * 14)[:N_REFS]
FAMILY = np.array([f"F{i % 7}" for i in range(N_REFS)])
PREFIX = np.array([f"P{i % 5}" for i in range(N_REFS)])
SCORES = (50.0, 55.5, 68.9)
IDENTITIES = (30.0, 40.0, 55.0)


def tied_hits(unit: str, n_pieces: int, seed: int) -> pd.DataFrame:
    """Return best-per-class rows with many ties between the classes of one query."""
    rng = np.random.default_rng(seed)
    rows = []
    for piece in range(n_pieces):
        for gene in range(12):
            score = float(rng.choice(SCORES))
            identity = float(rng.choice(IDENTITIES))
            for cls in ("BAC", "ARC", "EUK"):
                if rng.random() < 0.6:
                    ref = int(rng.choice(np.flatnonzero(cls == CLASS_OF_REF)))
                    name = f"{unit}|10|{piece}|{gene}"
                    rows.append((name, cls, ref, identity, score, score / 200))
    return pd.DataFrame(rows, columns=ref_hits.BEST_COLUMNS)


def piece_index(unit: str, n_pieces: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "contig_id": [f"{unit}|10|{p}" for p in range(n_pieces)],
            "n_genes": [12] * n_pieces,
        }
    )


def features(best: pd.DataFrame, index: pd.DataFrame) -> pd.DataFrame:
    out = ref_hits.piece_features(best, index, FAMILY, PREFIX)
    return out.sort_values("contig_id").reset_index(drop=True)


def test_piece_features_ignore_row_order() -> None:
    best = tied_hits("a", 60, seed=1)
    shuffled = best.sample(frac=1.0, random_state=7)

    pd.testing.assert_frame_equal(
        features(best, piece_index("a", 60)), features(shuffled, piece_index("a", 60))
    )


def test_piece_features_ignore_unrelated_rows() -> None:
    best = tied_hits("a", 60, seed=1)
    other = tied_hits("b", 60, seed=2)

    alone = features(best, piece_index("a", 60))
    together = features(pd.concat([best, other]), piece_index("a", 60))

    pd.testing.assert_frame_equal(alone, together)


@pytest.mark.parametrize("order", [(0, 1, 2), (2, 1, 0), (1, 2, 0)])
def test_best_per_class_prefers_identity_then_lower_reference(
    order: tuple[int, ...],
) -> None:
    ref_class = np.array(["BAC"] * 4)
    chunk = pd.DataFrame(
        {
            "qseqid": ["u|10|0|0"] * 3,
            "sseqid": ["r1_x", "r2_x", "r3_x"],
            "pident": [40.0, 45.0, 45.0],
            "qlen": [200] * 3,
            "bitscore": [60.0] * 3,
        }
    ).iloc[list(order)]

    best = ref_hits.best_per_class(chunk, ref_class)

    assert best["ref"].tolist() == [2]
    assert best["pident"].tolist() == [45.0]


def test_best_per_class_drops_hits_to_references_without_a_class() -> None:
    ref_class = np.array(["", "BAC", ""], dtype=object)
    chunk = pd.DataFrame(
        {
            "qseqid": ["u|1|0|0"] * 2,
            "sseqid": ["r0_x", "r1_x"],
            "pident": [99.0, 40.0],
            "qlen": [100, 100],
            "bitscore": [200.0, 50.0],
        }
    )

    best = ref_hits.best_per_class(chunk, ref_class)

    assert best["ref"].tolist() == [1]


@pytest.mark.parametrize("flip", [False, True])
def test_equal_score_and_identity_goes_to_the_first_class_name(flip: bool) -> None:
    best = pd.DataFrame(
        {
            "qseqid": ["u|10|0|0", "u|10|0|0"],
            "cls": ["BAC", "ARC"],
            "ref": [3, 5],
            "pident": [40.0, 40.0],
            "bitscore": [68.9, 68.9],
            "spa": [0.3445, 0.3445],
        }
    )
    if flip:
        best = best.iloc[::-1]

    out = ref_hits.piece_features(best, piece_index("u", 1), FAMILY, PREFIX)

    assert out["r_share_ARC"].tolist() == [1 / 12]
    assert out["r_share_BAC"].tolist() == [0.0]


def test_best_table_keeps_one_row_per_gene_and_class_in_a_fixed_order() -> None:
    best = pd.DataFrame(
        {
            "qseqid": ["u|10|0|0", "u|10|0|0", "u|10|0|0", "u|10|0|1"],
            "cls": ["BAC", "BAC", "ARC", "EUK"],
            "ref": [3, 5, 7, 9],
            "pident": [40.0, 55.0, 30.0, 20.0],
            "bitscore": [68.9, 90.0, 50.0, 40.0],
            "spa": [0.3, 0.45, 0.25, 0.2],
        }
    )

    table = ref_hits.best_table(best.iloc[::-1])

    assert table.columns.tolist() == ref_hits.BEST_COLUMNS
    assert len(table) == 3
    assert table[table["cls"] == "BAC"]["ref"].tolist() == [5]
    assert table["qseqid"].tolist() == sorted(table["qseqid"].tolist())


def test_a_window_without_a_hit_has_zero_shares_and_no_family() -> None:
    out = ref_hits.piece_features(
        pd.DataFrame(columns=ref_hits.BEST_COLUMNS), piece_index("u", 1), FAMILY, PREFIX
    )

    assert out["r_frac_hit"].tolist() == [0.0]
    assert out["r_entropy"].tolist() == [0.0]
    assert out["r_top_family_share"].isna().all()


def test_family_token_is_the_last_lineage_token_that_names_a_family() -> None:
    lineage = "d__Bacteria; p__X; f__Foo; g__Bar"

    assert ref_hits.family_token(lineage) == "f__Foo"
    assert ref_hits.family_token("Eukaryota; Fungi") == ""
    assert ref_hits.normalize_token("f__Foo") == "Foo"


def write_reference(path: Path, rows: list[tuple[int, str]]) -> None:
    lines = ["index\tclass\tlineage"] + [f"{i}\t{c}\tf__F{i}" for i, c in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_load_reference_empties_the_skipped_classes_and_appends_the_extras(
    tmp_path: Path,
) -> None:
    write_reference(tmp_path / "main.tsv", [(0, "NCLDV"), (1, "BAC")])
    write_reference(tmp_path / "extra.tsv", [(2, "NCLDV")])

    ref = ref_hits.load_reference(
        tmp_path / "main.tsv", tmp_path / "extra.tsv", ["NCLDV"]
    )

    assert ref["class"].fillna("").tolist() == ["", "BAC", "NCLDV"]


def test_read_hits_yields_nothing_for_an_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "hits.tsv"
    path.write_text("")

    assert list(ref_hits.read_hits(path)) == []
    assert ref_hits.collect_best([path], np.array(["BAC"], dtype=object)).empty


def test_reference_features_join_the_hit_files_of_both_databases(
    tmp_path: Path,
) -> None:
    write_reference(tmp_path / "ref.tsv", [(0, "BAC"), (1, "EUK")])
    write_reference(tmp_path / "extra.tsv", [(2, "NCLDV")])
    ref = ref_hits.load_reference(tmp_path / "ref.tsv", tmp_path / "extra.tsv", [])
    row = "{q}\tr{r}_x\t{p}\t100\t100\t100\t1e-9\t{b}\n"
    (tmp_path / "h1.tsv").write_text(row.format(q="u|1|0|0", r=0, p=50.0, b=120.0))
    (tmp_path / "h2.tsv").write_text(row.format(q="u|1|0|1", r=2, p=60.0, b=150.0))
    index = pd.DataFrame({"contig_id": ["u|1|0"], "n_genes": [2]})

    features, best = ref_hits.reference_features(
        ref, [tmp_path / "h1.tsv", tmp_path / "h2.tsv"], index
    )

    assert features["r_share_BAC"].tolist() == [0.5]
    assert features["r_share_NCLDV"].tolist() == [0.5]
    assert best["cls"].tolist() == ["BAC", "NCLDV"]
