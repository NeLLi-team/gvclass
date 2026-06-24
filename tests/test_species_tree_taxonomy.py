"""Unit tests for species-tree Section 5: tree-placement taxonomy.

Deterministic synthetic newicks exercise: patristic precompute parity vs ete3,
solo nearest-reference placement, shared multi-query clades (id/taxonomy/mean/
clade-ref distance), maximality, query-nested-among-refs = solo, one-row-per-query
/ no reference rows, lexicographic tie-breaks, and rooting invariance.
"""

from __future__ import annotations

from pathlib import Path

from ete3 import Tree

from src.core.species_tree.taxonomy_placement import (
    SpeciesTreeClassifier,
    _Patristic,
    format_lineage,
)
from src.core.species_tree.tree_inference import infer_species_tree


def _labels(tmp_path: Path) -> Path:
    path = tmp_path / "labels.tsv"
    path.write_text(
        "NCLDV__R1\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|Mimi|Mega|sp1\n"
        "NCLDV__R2\tNCLDV|Nucleocytoviricota|Megaviricetes|Algavirales|Phyco|Pra|sp2\n"
        "NCLDV__R3\tNCLDV|Nucleocytoviricota|Pokkesviricetes|Asfuvirales|Asf|Asfi|sp3\n"
    )
    return path


def _classifier(tmp_path: Path) -> SpeciesTreeClassifier:
    return SpeciesTreeClassifier(_labels(tmp_path))


def _write(tmp_path: Path, name: str, newick: str) -> Path:
    path = tmp_path / name
    path.write_text(newick.strip() + "\n")
    return path


def _row_key(row):
    return (
        row.query,
        row.nn_genome,
        row.nn_taxonomy,
        round(row.nn_distance, 9),
        row.clade_id,
        row.clade_members,
        None if row.clade_mean_dist is None else round(row.clade_mean_dist, 9),
        None if row.clade_ref_dist is None else round(row.clade_ref_dist, 9),
    )


def test_format_lineage() -> None:
    assert format_lineage("NCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|Mimi|Mega|sp1") == (
        "d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_Mimi;g_Mega;s_sp1"
    )


def test_patristic_parity_vs_ete3() -> None:
    tree = Tree(
        "((QUERY__q1:0.1,NCLDV__R1:0.2):0.3,(NCLDV__R2:0.15,NCLDV__R3:0.25):0.1);"
    )
    pat = _Patristic(tree)
    leaves = tree.get_leaves()
    for a in leaves:
        for b in leaves:
            assert abs(pat.dist(a, b) - tree.get_distance(a, b)) < 1e-9


def test_solo_query_nearest_reference(tmp_path: Path) -> None:
    tree = _write(
        tmp_path,
        "solo.treefile",
        "((QUERY__q1:0.1,NCLDV__R1:0.1):0.2,(NCLDV__R2:0.3,NCLDV__R3:0.3):0.1);",
    )
    rows = _classifier(tmp_path).classify(tree)
    assert len(rows) == 1
    row = rows[0]
    assert row.query == "q1"
    assert row.nn_genome == "NCLDV__R1"  # sibling, nearest
    assert row.nn_taxonomy.startswith("d_NCLDV;p_Nucleocytoviricota")
    assert "o_Imitervirales" in row.nn_taxonomy
    # Solo -> clade fields blank.
    assert row.clade_id == "" and row.clade_members == ""
    assert row.clade_mean_dist is None and row.clade_ref_dist is None


def test_two_query_clade_shares_taxonomy_and_metrics(tmp_path: Path) -> None:
    tree = _write(
        tmp_path,
        "clade.treefile",
        "((QUERY__q1:0.1,QUERY__q2:0.1):0.2,(NCLDV__R1:0.1,NCLDV__R2:0.2):0.1);",
    )
    rows = _classifier(tmp_path).classify(tree)
    assert {r.query for r in rows} == {"q1", "q2"}
    # Both members share clade id, members, taxonomy, and clade metrics.
    assert {r.clade_id for r in rows} == {"q1"}
    assert {r.clade_members for r in rows} == {"q1;q2"}
    assert len({r.nn_taxonomy for r in rows}) == 1
    assert len({round(r.clade_mean_dist, 9) for r in rows}) == 1
    assert len({round(r.clade_ref_dist, 9) for r in rows}) == 1
    # mean pairwise distance of the two queries = 0.1 + 0.1.
    assert abs(rows[0].clade_mean_dist - 0.2) < 1e-9


def test_maximal_clade_not_subclades(tmp_path: Path) -> None:
    tree = _write(
        tmp_path,
        "maximal.treefile",
        "(((QUERY__q1:0.1,QUERY__q2:0.1):0.1,QUERY__q3:0.2):0.3,"
        "(NCLDV__R1:0.1,NCLDV__R2:0.1):0.1);",
    )
    rows = _classifier(tmp_path).classify(tree)
    # One maximal clade of all three queries, not nested sub-clades.
    assert {r.query for r in rows} == {"q1", "q2", "q3"}
    assert {r.clade_members for r in rows} == {"q1;q2;q3"}
    assert {r.clade_id for r in rows} == {"q1"}


def test_query_nested_among_references_is_solo(tmp_path: Path) -> None:
    # q1 sits between references -> its maximal all-query subtree is itself.
    tree = _write(
        tmp_path,
        "nested.treefile",
        "((NCLDV__R1:0.1,QUERY__q1:0.1):0.1,(NCLDV__R2:0.1,NCLDV__R3:0.1):0.1);",
    )
    rows = _classifier(tmp_path).classify(tree)
    assert len(rows) == 1
    assert rows[0].query == "q1"
    assert rows[0].clade_id == ""  # solo


def test_one_row_per_query_no_reference_rows(tmp_path: Path) -> None:
    tree = _write(
        tmp_path,
        "mix.treefile",
        "((QUERY__q1:0.1,QUERY__q2:0.1):0.2,((NCLDV__R1:0.1,NCLDV__R2:0.1):0.1,QUERY__q3:0.3):0.1);",
    )
    rows = _classifier(tmp_path).classify(tree)
    queries = {r.query for r in rows}
    # Exactly the three queries, one row each; references never produce a row.
    assert queries == {"q1", "q2", "q3"}
    assert len(rows) == 3


def test_no_references_returns_empty(tmp_path: Path) -> None:
    tree = _write(
        tmp_path,
        "norefs.treefile",
        "(QUERY__q1:0.1,(QUERY__q2:0.1,QUERY__q3:0.1):0.1);",
    )
    assert _classifier(tmp_path).classify(tree) == []


def test_single_reference_places_all_queries(tmp_path: Path) -> None:
    tree = _write(
        tmp_path,
        "oneref.treefile",
        "((QUERY__q1:0.1,QUERY__q2:0.1):0.2,NCLDV__R1:0.3);",
    )
    rows = _classifier(tmp_path).classify(tree)
    assert {r.query for r in rows} == {"q1", "q2"}
    assert all(r.nn_genome == "NCLDV__R1" for r in rows)


def test_set_outgroup_failure_fails_closed(tmp_path: Path, monkeypatch) -> None:
    tree = _write(
        tmp_path,
        "failreroot.treefile",
        "((QUERY__q1:0.1,NCLDV__R1:0.1):0.2,(NCLDV__R2:0.1,NCLDV__R3:0.1):0.1);",
    )

    def boom(self, *args, **kwargs):
        raise RuntimeError("forced reroot failure")

    monkeypatch.setattr(Tree, "set_outgroup", boom)
    # Fail closed: no placement rather than one on an arbitrary rooting.
    assert _classifier(tmp_path).classify(tree) == []


def test_placement_is_rooting_invariant(tmp_path: Path) -> None:
    newick = (
        "((QUERY__q1:0.1,QUERY__q2:0.1):0.2,"
        "(NCLDV__R1:0.1,(NCLDV__R2:0.2,NCLDV__R3:0.2):0.1):0.1);"
    )
    t1 = Tree(newick)
    file1 = tmp_path / "root1.treefile"
    t1.write(outfile=str(file1))

    # A different rooting of the SAME topology.
    t2 = Tree(newick)
    t2.set_outgroup(t2.search_nodes(name="NCLDV__R3")[0])
    file2 = tmp_path / "root2.treefile"
    t2.write(outfile=str(file2))

    clf = _classifier(tmp_path)
    rows1 = [_row_key(r) for r in clf.classify(file1)]
    rows2 = [_row_key(r) for r in clf.classify(file2)]
    assert rows1 == rows2


def test_infer_species_tree_parses(tmp_path: Path) -> None:
    faa = tmp_path / "sm.supermatrix.faa"
    faa.write_text(
        ">QUERY__q1\nMKTAYIAKQRQISFVKSHFSR\n"
        ">NCLDV__R1\nMKTAYIAKQRQISFVKSHFSR\n"
        ">NCLDV__R2\nMKQAYIAKQRQISFVKSHFSR\n"
        ">NCLDV__R3\nMKTAYIWKQRQISFVKSHFSR\n"
    )
    out = tmp_path / "sm.treefile"
    infer_species_tree(faa, out, threads=1)
    tree = Tree(str(out))
    assert {leaf.name for leaf in tree.get_leaves()} == {
        "QUERY__q1",
        "NCLDV__R1",
        "NCLDV__R2",
        "NCLDV__R3",
    }


def test_nearest_reference_tie_break_is_rooting_invariant(tmp_path: Path) -> None:
    """q1 is equidistant (0.207491) from R1 and R2; the 9-dp-rounded tie-break must
    pick the smaller id (R1) regardless of the input tree's rooting (guards the FP
    non-associativity the §5 review demonstrated)."""
    newick = (
        "((QUERY__q1:0.0,(NCLDV__R1:0.207491,NCLDV__R2:0.207491):0.0):0.2,"
        "NCLDV__R3:0.3);"
    )
    t1 = Tree(newick)
    file1 = tmp_path / "ti1.treefile"
    t1.write(outfile=str(file1))

    t2 = Tree(newick)
    t2.set_outgroup(t2.search_nodes(name="NCLDV__R3")[0])
    file2 = tmp_path / "ti2.treefile"
    t2.write(outfile=str(file2))

    clf = _classifier(tmp_path)
    nn1 = clf.classify(file1)[0].nn_genome
    nn2 = clf.classify(file2)[0].nn_genome
    assert nn1 == nn2 == "NCLDV__R1"


def test_zero_length_branches_and_lexicographic_ties(tmp_path: Path) -> None:
    # q1 equidistant (0.0) from R1 and R2 via zero-length branches -> tie broken
    # by smallest reference id (R1 < R2).
    tree = _write(
        tmp_path,
        "ties.treefile",
        "((QUERY__q1:0.0,(NCLDV__R1:0.0,NCLDV__R2:0.0):0.0):0.1,NCLDV__R3:0.2);",
    )
    rows = _classifier(tmp_path).classify(tree)
    assert len(rows) == 1
    assert rows[0].nn_genome == "NCLDV__R1"  # lexicographic tie-break
    assert rows[0].nn_distance == 0.0
