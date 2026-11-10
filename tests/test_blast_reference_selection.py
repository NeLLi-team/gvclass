"""Reference selection gives every query protein of a marker its own best hits."""

import pytest

from src.core.blast import BlastHit, select_references


def _hit(
    subject: str,
    score: float,
    *,
    evalue: float = 1e-30,
    identity: float = 80.0,
) -> BlastHit:
    return {
        "query_id": "query",
        "target_id": subject,
        "identity": identity,
        "aln_len": 100,
        "mismatches": 0,
        "gap_opens": 0,
        "q_start": 1,
        "q_end": 100,
        "s_start": 1,
        "s_end": 100,
        "evalue": evalue,
        "score": score,
    }


def test_select_references_gives_a_weak_second_copy_its_own_neighbours() -> None:
    long_copy = [_hit(f"long_{i:03d}", 1000 - i) for i in range(150)]
    short_copy = [_hit(f"short_{i:02d}", 100 - i) for i in range(20)]

    selected = select_references({"long": long_copy, "short": short_copy}, max_refs=100)

    assert len(selected) == 100
    assert selected[:4] == ["long_000", "short_00", "long_001", "short_01"]
    assert len([ref for ref in selected if ref.startswith("short_")]) == 20
    assert len([ref for ref in selected if ref.startswith("long_")]) == 80


def test_select_references_with_one_protein_keeps_score_order() -> None:
    hits = [_hit("a", 50), _hit("b", 90), _hit("c", 70)]

    assert select_references({"query": hits}, max_refs=100) == ["b", "c", "a"]


def test_select_references_with_one_protein_cuts_at_max_refs() -> None:
    hits = [_hit(f"ref_{i:03d}", 500 - i) for i in range(150)]

    selected = select_references({"query": hits}, max_refs=100)

    assert len(selected) == 100
    assert selected[0] == "ref_000"
    assert selected[-1] == "ref_099"


def test_select_references_orders_one_proteins_hits_by_score_then_evalue_then_id() -> (
    None
):
    hits = [
        _hit("d", 80, evalue=1e-20),
        _hit("c", 80, evalue=1e-30),
        _hit("b", 90),
        _hit("a", 80, evalue=1e-20),
    ]

    assert select_references({"query": hits}, max_refs=100) == ["b", "c", "a", "d"]


@pytest.mark.parametrize("subject_order", [["c", "a", "b"], ["b", "c", "a"]])
def test_select_references_breaks_ties_by_subject_id(
    subject_order: list[str],
) -> None:
    hits = [_hit(subject, 70) for subject in subject_order]

    assert select_references({"query": hits}, max_refs=100) == ["a", "b", "c"]


def test_select_references_ignores_the_order_of_query_proteins() -> None:
    first = [_hit("x", 90), _hit("y", 80)]
    second = [_hit("z", 90), _hit("w", 80)]

    forward = select_references({"p1": first, "p2": second}, max_refs=100)
    backward = select_references({"p2": second, "p1": first}, max_refs=100)

    assert forward == backward == ["x", "z", "y", "w"]


def test_select_references_skips_subjects_that_are_already_chosen() -> None:
    first = [_hit("shared", 90), _hit("only_first", 80)]
    second = [_hit("shared", 95), _hit("only_second", 60)]

    selected = select_references({"p1": first, "p2": second}, max_refs=100)

    assert selected == ["shared", "only_first", "only_second"]


def test_select_references_stops_at_max_refs_inside_a_rank() -> None:
    hits_by_query = {
        "p1": [_hit("a", 90)],
        "p2": [_hit("b", 90)],
        "p3": [_hit("c", 90)],
    }

    assert select_references(hits_by_query, max_refs=2) == ["a", "b"]


def test_select_references_drops_hits_below_the_identity_floor() -> None:
    hits = [_hit("close", 90, identity=30.0), _hit("far", 95, identity=29.9)]

    assert select_references({"query": hits}, max_refs=100) == ["close"]


def test_select_references_returns_every_hit_when_there_are_fewer_than_max_refs() -> (
    None
):
    hits = [_hit(f"ref_{i:02d}", 100 - i) for i in range(15)]

    assert len(select_references({"query": hits}, max_refs=100)) == 15


def test_select_references_without_hits_is_empty() -> None:
    assert select_references({}, max_refs=100) == []
