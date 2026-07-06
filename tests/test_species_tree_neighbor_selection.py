"""Unit tests for species-tree Section 2: top-K NCLDV neighbor selection.

Covers ``TreeAnalyzer.get_top_k_neighbors`` (tree-based, prefix-filtered,
genome-deduplicated, deterministic) and the neighbor_selection primitives
(``select_ncldv_neighbors_from_trees``, ``dereplicate_reference_genomes``,
``build_ncldv_reference_subset``) using deterministic synthetic newick trees and
tiny FASTA fixtures. The dedicated search + gene-tree builders are exercised with
a small real pyswrd/FAMSA/VeryFastTree run.
"""

from __future__ import annotations

from pathlib import Path

from src.core.tree_analysis import TreeAnalyzer
from src.core.species_tree.neighbor_selection import (
    build_ncldv_reference_subset,
    build_species_tree_gene_tree,
    gather_ncldv_neighbor_proteins,
    select_ncldv_neighbors_from_trees,
)


def _labels_file(tmp_path: Path) -> Path:
    """A minimal labels.tsv mapping the synthetic NCLDV genomes to lineages."""
    path = tmp_path / "labels.tsv"
    lines = []
    for g in ("G1", "G2", "G3", "G4"):
        lines.append(
            f"NCLDV__{g}\tNCLDV|Nucleocytoviricota|Imitervirales|fam|gen|sp_{g}"
        )
    lines.append("EUK-pEVE__E1\t" "EUK-pEVE|Discosea-pEVE|Flabellinia-pEVE|o|f|g|s")
    lines.append("BAC__B1\tBAC|Pseudomonadota|Gamma|Ent|Esc|coli")
    path.write_text("\n".join(lines) + "\n")
    return path


def _analyzer(tmp_path: Path) -> TreeAnalyzer:
    return TreeAnalyzer(_labels_file(tmp_path))


def _write_tree(tmp_path: Path, name: str, newick: str) -> Path:
    path = tmp_path / name
    path.write_text(newick.strip() + "\n")
    return path


# Query QUERY_1 distances: NCLDV__G1_1=0.3, NCLDV__G2_1=0.4, NCLDV__G3_1=0.9,
# BAC__B1_1 filtered. NCLDV__G1_2 is a far paralog of G1 (dist 0.9; dedup target).
_MIXED_TREE = (
    "((QUERY_1:0.1,(NCLDV__G1_1:0.1,NCLDV__G2_1:0.2):0.1):0.1,"
    "((NCLDV__G3_1:0.5,BAC__B1_1:0.3):0.1,NCLDV__G1_2:0.6):0.1);"
)


def test_top_k_ordering_and_ncldv_only_filter(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(tmp_path, "mixed.treefile", _MIXED_TREE)

    hits = analyzer.get_top_k_neighbors(tree, "QUERY", k=50, keep_prefix="NCLDV__")
    genomes = [genome_id for _, genome_id, _ in hits]

    # NCLDV-only (BAC excluded), ordered by ascending distance, genome-deduped.
    assert genomes == ["NCLDV__G1", "NCLDV__G2", "NCLDV__G3"]
    # Distances strictly increasing in the documented order.
    dists = [d for _, _, d in hits]
    assert dists == sorted(dists)
    # The kept G1 protein is the nearer paralog (G1_1, not the far G1_2).
    g1_protein = next(p for p, g, _ in hits if g == "NCLDV__G1")
    assert g1_protein == "NCLDV__G1_1"


def test_top_k_truncates_to_k(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(tmp_path, "mixed.treefile", _MIXED_TREE)

    hits = analyzer.get_top_k_neighbors(tree, "QUERY", k=2, keep_prefix="NCLDV__")
    assert [g for _, g, _ in hits] == ["NCLDV__G1", "NCLDV__G2"]


def test_top_k_accepts_auxiliary_peve_prefix(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(
        tmp_path,
        "peve.treefile",
        "((QUERY_1:0.1,EUK-pEVE__E1|p_1:0.1):0.1,"
        "(NCLDV__G1_1:0.3,EUK__ORDINARY|p_1:0.1):0.1);",
    )

    hits = analyzer.get_top_k_neighbors(
        tree,
        "QUERY",
        k=50,
        keep_prefix=("NCLDV__", "EUK-pEVE__"),
    )

    genomes = [genome_id for _, genome_id, _ in hits]
    assert genomes == ["EUK-pEVE__E1", "NCLDV__G1"]
    assert "EUK__ORDINARY" not in genomes


def test_fewer_than_k_returns_all(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(tmp_path, "mixed.treefile", _MIXED_TREE)
    hits = analyzer.get_top_k_neighbors(tree, "QUERY", k=50, keep_prefix="NCLDV__")
    assert len(hits) == 3  # only three distinct NCLDV genomes present


def test_no_ncldv_neighbors_returns_empty(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(
        tmp_path,
        "noncldv.treefile",
        "(QUERY_1:0.1,(BAC__B1_1:0.2,BAC__B2_1:0.3):0.1);",
    )
    assert analyzer.get_top_k_neighbors(tree, "QUERY", keep_prefix="NCLDV__") == []


def test_no_query_leaf_returns_empty(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(
        tmp_path,
        "noquery.treefile",
        "(NCLDV__G1_1:0.1,(NCLDV__G2_1:0.2,NCLDV__G3_1:0.3):0.1);",
    )
    assert analyzer.get_top_k_neighbors(tree, "QUERY", keep_prefix="NCLDV__") == []


def test_missing_and_skipped_tree_return_empty(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    assert analyzer.get_top_k_neighbors(tmp_path / "absent.treefile", "QUERY") == []

    present = _write_tree(tmp_path, "skip.treefile", _MIXED_TREE)
    present.with_suffix(".SKIPPED").write_text("skipped\n")
    assert analyzer.get_top_k_neighbors(present, "QUERY") == []


def test_query_prefixed_query_is_not_a_neighbor(tmp_path: Path) -> None:
    """Drift #4: a query stem may itself start with the keep_prefix; query leaves
    must still never be returned as neighbors."""
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(
        tmp_path,
        "qprefix.treefile",
        "((NCLDV__QUERY_1:0.1,NCLDV__G1_1:0.1):0.1,NCLDV__G2_1:0.3);",
    )
    hits = analyzer.get_top_k_neighbors(
        tree, "NCLDV__QUERY", k=50, keep_prefix="NCLDV__"
    )
    genomes = [g for _, g, _ in hits]
    assert "NCLDV__QUERY" not in genomes
    assert set(genomes) == {"NCLDV__G1", "NCLDV__G2"}


def test_select_across_trees(tmp_path: Path) -> None:
    analyzer = _analyzer(tmp_path)
    # Two group trees that share genome G2; G2 nearer in tree B.
    tree_a = _write_tree(
        tmp_path,
        "grpA.treefile",
        "((QUERY_1:0.1,NCLDV__G1_1:0.1):0.1,NCLDV__G2_1:0.5);",
    )
    tree_b = _write_tree(
        tmp_path,
        "grpB.treefile",
        "((QUERY_1:0.1,NCLDV__G2_1:0.1):0.1,NCLDV__G4_1:0.4);",
    )
    hits = select_ncldv_neighbors_from_trees(
        {"grpA": tree_a, "grpB": tree_b}, analyzer, "QUERY", k=50
    )

    # Group attribution preserved across trees.
    assert {h.group for h in hits} == {"grpA", "grpB"}
    by_group_genomes = {(h.group, h.genome_id) for h in hits}
    assert ("grpA", "NCLDV__G1") in by_group_genomes
    assert ("grpB", "NCLDV__G4") in by_group_genomes


def test_top_k_resolves_genome_from_pipe_headers(tmp_path: Path) -> None:
    """Real reference headers are '>NCLDV__<genome>|<protein>'; the genome must
    resolve via the '|'-split path of _extract_genome_id (synthetic '_'-only
    headers exercise a different branch)."""
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "NCLDV__IMGM3300009076_BIN259\tNCLDV|Nucleocytoviricota|Imitervirales|f|g|s1\n"
        "NCLDV__bin_6581\tNCLDV|Nucleocytoviricota|Pimascovirales|f|g|s2\n"
    )
    analyzer = TreeAnalyzer(labels)
    # Two proteins of the SAME genome BIN259 (paralogs) + one of bin_6581.
    tree = _write_tree(
        tmp_path,
        "pipe.treefile",
        "((PkV-RF01|contig_1_5:0.1,"
        "NCLDV__IMGM3300009076_BIN259|Ga0115550_1009924_4:0.1):0.1,"
        "(NCLDV__IMGM3300009076_BIN259|Ga0115550_1009924_9:0.6,"
        "NCLDV__bin_6581|Ga0590176_00006195_36:0.3):0.1);",
    )
    hits = analyzer.get_top_k_neighbors(tree, "PkV-RF01", k=50, keep_prefix="NCLDV__")
    genomes = [g for _, g, _ in hits]
    # Both BIN259 paralogs collapse to one genome; bin_6581 distinct.
    assert genomes.count("NCLDV__IMGM3300009076_BIN259") == 1
    assert set(genomes) == {"NCLDV__IMGM3300009076_BIN259", "NCLDV__bin_6581"}
    # The nearer BIN259 paralog (the _4 protein at dist 0.3) is kept.
    bin259_protein = next(p for p, g, _ in hits if g == "NCLDV__IMGM3300009076_BIN259")
    assert bin259_protein.endswith("_4")


def test_build_ncldv_reference_subset_filters_prefix(tmp_path: Path) -> None:
    group_faa = tmp_path / "grp.faa"
    group_faa.write_text(
        ">NCLDV__G1|c_1\nMKTAYIAK\n"
        ">BAC__B1|c_1\nMKQAYIAK\n"
        ">NCLDV__G2|c_1\nMKTAYIWK\n"
        ">PPV__P1|c_1\nMKTAYIAR\n"
    )
    out = tmp_path / "grp.ncldv.faa"
    n = build_ncldv_reference_subset(group_faa, out)
    assert n == 2
    headers = [
        ln[1:].strip() for ln in out.read_text().splitlines() if ln.startswith(">")
    ]
    assert all(h.startswith("NCLDV__") for h in headers)
    assert len(headers) == 2


def test_build_reference_subset_filters_multiple_prefixes(tmp_path: Path) -> None:
    group_faa = tmp_path / "grp.faa"
    group_faa.write_text(
        ">NCLDV__G1|c_1\nMKTAYIAK\n"
        ">EUK-pEVE__E1|c_1\nMKQAYIAK\n"
        ">EUK__ORDINARY|c_1\nMKTAYIWK\n"
        ">PPV__P1|c_1\nMKTAYIAR\n"
    )
    out = tmp_path / "grp.panel.faa"

    n = build_ncldv_reference_subset(
        group_faa,
        out,
        prefix=("NCLDV__", "EUK-pEVE__"),
    )

    assert n == 2
    headers = [
        ln[1:].strip() for ln in out.read_text().splitlines() if ln.startswith(">")
    ]
    assert headers == ["NCLDV__G1|c_1", "EUK-pEVE__E1|c_1"]


def test_dedicated_search_and_tree_build(tmp_path: Path) -> None:
    """Light end-to-end of the breadth mechanism: NCLDV-only BLAST -> dedicated
    gene tree -> top-K, on tiny real inputs (pyswrd + FAMSA + VeryFastTree)."""
    # One query protein clearly closest to G1; G2/G3 diverged; B1 is non-NCLDV.
    base = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA"
    g1 = base
    g2 = base.replace("EERLGL", "DDRLGL").replace("VGDGT", "VGEGT")
    g3 = base.replace("AKQRQ", "AKQSQ").replace("SGAEK", "SGADK")
    query = tmp_path / "query.faa"
    query.write_text(f">QUERY_1\n{base}\n")
    ncldv_ref = tmp_path / "grp.ncldv.faa"
    ncldv_ref.write_text(
        f">NCLDV__G1_1\n{g1}\n>NCLDV__G2_1\n{g2}\n>NCLDV__G3_1\n{g3}\n"
    )

    candidates = gather_ncldv_neighbor_proteins(
        query, ncldv_ref, tmp_path / "q.m8", max_candidates=50, threads=1
    )
    assert "NCLDV__G1_1" in candidates

    tree = build_species_tree_gene_tree(
        query, candidates, ncldv_ref, tmp_path / "grp.treefile", threads=1
    )
    assert tree is not None and tree.exists()

    analyzer = _analyzer(tmp_path)
    hits = analyzer.get_top_k_neighbors(tree, "QUERY", k=50, keep_prefix="NCLDV__")
    assert hits, "expected at least one NCLDV neighbor from the dedicated tree"
    # G1 (identical to query) is the nearest neighbor.
    assert hits[0][1] == "NCLDV__G1"


def test_build_species_tree_gene_tree_too_few_sequences(tmp_path: Path) -> None:
    query = tmp_path / "query.faa"
    query.write_text(">QUERY_1\nMKTAYIAK\n")
    ncldv_ref = tmp_path / "grp.ncldv.faa"
    ncldv_ref.write_text(">NCLDV__G1_1\nMKTAYIAK\n")
    # Only 2 sequences total (< min 3) -> None.
    out = build_species_tree_gene_tree(
        query, ["NCLDV__G1_1"], ncldv_ref, tmp_path / "g.treefile", threads=1
    )
    assert out is None


def test_exact_distance_tie_breaks_by_genome_id(tmp_path: Path) -> None:
    """Two NCLDV genomes equidistant from the query resolve deterministically by
    genome id (the reason for the (distance, genome_id) sort key)."""
    analyzer = _analyzer(tmp_path)
    tree = _write_tree(
        tmp_path,
        "tie.treefile",
        "(QUERY_1:0.1,(NCLDV__G1_1:0.2,NCLDV__G2_1:0.2):0.1);",
    )
    hits = analyzer.get_top_k_neighbors(tree, "QUERY", k=50, keep_prefix="NCLDV__")
    # Identical patristic distance (0.4) for both -> ordered G1 before G2 by id.
    assert [g for _, g, _ in hits] == ["NCLDV__G1", "NCLDV__G2"]
    assert hits[0][2] == hits[1][2]  # genuinely tied distances


# 60-aa protein families for the orchestration end-to-end test (pyswrd needs
# enough length/identity to clear the 1e-5 e-value gate).
_FAM_A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA"
_FAM_A2 = "MKTAYIAKQRQISFVKSHFSRQLDERLGLIEVQAPILSRVGDGTQDNLSGAEKAVAVKVKA"
_FAM_B = "MSEQNNTSGFLGKKVDLSSLTGKKVAVDASHALYQFLIAVRQEGGQLTNEAGETTSHLMGM"
_FAM_B2 = "MSEQNNTSGFLGKKVDLSSLTGKKVAVDASHALYQFLIAVRQEGGQLTNEAGDTTSHLLGM"


def test_build_query_group_trees_and_select_orchestration(
    tmp_path, monkeypatch
) -> None:
    """End-to-end: subset -> dedicated search -> tree -> top-K across two groups,
    with shared subset_cache reuse and skip handling for empty/missing groups.
    Exercises the live primitives the per-query hook composes."""
    import src.core.species_tree.neighbor_selection as ns

    query_hits = tmp_path / "query_hits_faa"
    query_hits.mkdir()
    ref_dir = tmp_path / "faa"
    ref_dir.mkdir()

    # grpA / grpB: query hit + NCLDV refs (one identical, one diverged) + a BAC.
    (query_hits / "grpA.faa").write_text(f">QUERYGENOME_1\n{_FAM_A}\n")
    (ref_dir / "grpA.faa").write_text(
        f">NCLDV__A1_1\n{_FAM_A}\n>NCLDV__A2_1\n{_FAM_A2}\n>BAC__BA_1\n{_FAM_B}\n"
    )
    (query_hits / "grpB.faa").write_text(f">QUERYGENOME_1\n{_FAM_B}\n")
    (ref_dir / "grpB.faa").write_text(
        f">NCLDV__B1_1\n{_FAM_B}\n>NCLDV__B2_1\n{_FAM_B2}\n>BAC__BB_1\n{_FAM_A}\n"
    )
    # grpC: query hit present but NO reference faa -> skipped.
    (query_hits / "grpC.faa").write_text(f">QUERYGENOME_1\n{_FAM_A}\n")
    # grpEmpty: empty query faa -> skipped.
    (query_hits / "grpEmpty.faa").write_text("")
    (ref_dir / "grpEmpty.faa").write_text(f">NCLDV__E1_1\n{_FAM_A}\n")

    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "\n".join(
            f"NCLDV__{g}\tNCLDV|Nucleocytoviricota|Imitervirales|f|g|s_{g}"
            for g in ("A1", "A2", "B1", "B2", "E1")
        )
        + "\n"
    )
    analyzer = TreeAnalyzer(labels)

    # Count how often the NCLDV subset is actually (re)built.
    calls = {"n": 0}
    real_subset = ns.build_ncldv_reference_subset

    def _counting_subset(group_faa, out_faa, prefix=ns.NCLDV_PREFIX):
        calls["n"] += 1
        return real_subset(group_faa, out_faa, prefix=prefix)

    monkeypatch.setattr(ns, "build_ncldv_reference_subset", _counting_subset)

    groups = ["grpA", "grpB", "grpC", "grpEmpty"]
    cache: dict = {}
    work = tmp_path / "work"

    group_to_tree = ns.build_query_group_trees(
        query_hits,
        ref_dir,
        groups,
        "QUERYGENOME",
        work,
        threads=1,
        subset_cache=cache,
    )
    hits = ns.select_ncldv_neighbors_from_trees(
        group_to_tree,
        analyzer,
        "QUERYGENOME",
        k=50,
    )

    # Aggregated across the two viable groups; grpC (no ref) and grpEmpty skipped.
    assert {h.group for h in hits} == {"grpA", "grpB"}
    assert {h.genome_id for h in hits} >= {"NCLDV__A1", "NCLDV__B1"}
    # Subset built once per viable group; cache populated.
    assert {key[0] for key in cache} == {"grpA", "grpB"}
    assert calls["n"] == 2

    # Second query reuses the cached subsets -> no extra subset builds.
    (query_hits / "grpA.faa").write_text(f">QUERYGENOME2_1\n{_FAM_A2}\n")
    group_to_tree2 = ns.build_query_group_trees(
        query_hits,
        ref_dir,
        groups,
        "QUERYGENOME2",
        work,
        threads=1,
        subset_cache=cache,
    )
    hits2 = ns.select_ncldv_neighbors_from_trees(
        group_to_tree2,
        analyzer,
        "QUERYGENOME2",
        k=50,
    )
    assert calls["n"] == 2  # unchanged: subsets reused from cache
    assert any(h.group == "grpA" for h in hits2)
