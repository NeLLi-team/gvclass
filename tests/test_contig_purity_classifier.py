"""Unit tests for the v1.4.3 Phase 2 per-contig taxonomic-purity classifier.

Covers the six classifier-level cases listed in
docs/plans/active/v1_4_3_contig_purity_contamination.md under the
"Tests" section.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.config.marker_sets import BUSCO_MODELS, UNI56_MODELS
from src.core.contamination_scoring import _marker_category
from src.core.summarize.paths_stats import QueryEvidence
from tests.conftest import skip_if_no_runtime_resources

pytestmark = skip_if_no_runtime_resources()


def _make_scorer(tmp_path: Path):
    """Build a ContaminationScorer with a synthetic labels table covering
    the cellular lineages the tests reference."""
    from src.core.contamination_scoring import ContaminationScorer
    from tests.conftest import stage_db_resources

    db = tmp_path / "db"
    stage_db_resources(db, markers=False)
    labels = db / "labels.tsv"
    # LabelIndex expects genome_id<TAB>domain|phylum|class|order|family|genus|species.
    # The domain of a lineage is the ID prefix before "__".
    labels.write_text(
        "EUK__Nannochloropsis_A\tEuk|Stramenopiles|Eustigmatophyceae|Eustigmatales|Monodopsidaceae|Nannochloropsis|Nannochloropsis_gaditana\n"
        "EUK__Nannochloropsis_B\tEuk|Stramenopiles|Eustigmatophyceae|Eustigmatales|Monodopsidaceae|Nannochloropsis|Nannochloropsis_oceanica\n"
        "EUK__Acanthocystis_A\tEuk|Centroheliozoa|Haptista|Acanthocystidae|Acanthocystidae|Acanthocystis|NA\n"
        "EUK__NA_A\tEuk|NA|NA|NA|NA|NA|NA\n"
        "EUK__Other_A\tEuk|NA|NA|OtherOrder|NA|NA|NA\n"
        "EUK__NoOrder_A\tEuk|Stramenopiles|Eustigmatophyceae||Monodopsidaceae|Nannochloropsis|x\n"
        "EUK-pEVE__Host_A\tEuk-pEVE|NA|NA|NA|NA|NA|NA\n"
        "NCLDV__Mimi_1\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|Mimiviridae|Mimivirus|Acanthamoeba_polyphaga_mimivirus\n"
        "BAC__Eco_A\tBac|Pseudomonadota|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|Escherichia_coli\n"
    )
    # Minimal completeness table (not used by the classifier).
    (db / "order_completeness.tab").write_text(
        "Order\tOrthogroups\tAverage_Percent\tStd_Percent\n"
    )
    return ContaminationScorer(db)


def _tree_nn(marker: str, protein: str, neighbor: str, distance: float = 0.5):
    """Build a single-row tree_nn_results slice."""
    return {marker: {protein: {neighbor: distance}}}


def _merge(*slices):
    out: dict = {}
    for slice_ in slices:
        for marker, proteins in slice_.items():
            out.setdefault(marker, {}).update(proteins)
    return out


def _blast_hit(
    query_id: str,
    subject_domain: str,
    subject_order: str = "",
    identity: float = 80.0,
    score: float = 120.0,
):
    from src.core.contamination_scoring import BlastHit

    return BlastHit(
        query_id=query_id,
        subject_id=f"{subject_domain}__synthetic",
        identity=identity,
        score=score,
        subject_domain=subject_domain,
        subject_order=subject_order,
        subject_family="",
    )


# ---------------------------------------------------------------------------
# Six classifier cases from the plan
# ---------------------------------------------------------------------------


def test_few_cellular_genes_on_viral_contig_are_not_removal_evidence(
    tmp_path: Path,
) -> None:
    """Two cellular-like genes stay below the coherent-contamination floor."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M_NCLDV", "p_viral", "NCLDV__Mimi_1"),
        _tree_nn("M_EUK_1", "p_euk_1", "EUK__Nannochloropsis_A"),
        _tree_nn("M_EUK_2", "p_euk_2", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        "p_viral": _blast_hit("p_viral", "NCLDV", "Imitervirales"),
        "p_euk_1": _blast_hit("p_euk_1", "EUK", "Eustigmatales", identity=90.0),
        "p_euk_2": _blast_hit("p_euk_2", "EUK", "Eustigmatales", identity=90.0),
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_1",
        contig_proteins=["p_viral", "p_euk_1", "p_euk_2"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "viral_bearing"
    assert result["n_viral_proteins"] == 1
    assert result["cellular_evidence_coherent"] is False


def test_coherent_cellular_and_viral_evidence_stays_mixed(tmp_path: Path) -> None:
    """Five coherent cellular hits plus one viral hit retain both signals."""
    scorer = _make_scorer(tmp_path)
    cellular = [f"p_cell_{i}" for i in range(5)]
    tree = _merge(
        *[
            _tree_nn(f"M_cell_{i}", protein, "EUK__Nannochloropsis_A")
            for i, protein in enumerate(cellular)
        ],
        _tree_nn("M_viral", "p_viral", "NCLDV__Mimi_1"),
    )
    best_hits = {
        **{
            protein: _blast_hit(protein, "EUK", "Eustigmatales", identity=90.0)
            for protein in cellular
        },
        "p_viral": _blast_hit("p_viral", "NCLDV", "Imitervirales", identity=90.0),
    }
    proteins = [*cellular, "p_viral"]

    result = scorer._classify_contig_purity(
        contig_id="contig_mixed",
        contig_proteins=proteins,
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "mixed_evidence"
    assert result["n_cellular_proteins"] == 5
    assert result["n_viral_proteins"] == 1

    features = scorer.compute_contig_purity_features(
        protein_to_contig={protein: "contig_mixed" for protein in proteins},
        contig_lengths={"contig_mixed": 1000},
        tree_nn_results=tree,
        best_hits=best_hits,
        attribution_mode="fna_gene_calling",
    ).features
    assert features["cellular_coherent_contig_count"] == 1
    assert features["viral_bearing_contig_count"] == 1
    assert features["mixed_evidence_contig_count"] == 1


def test_peve_hits_are_reported_without_diluting_cellular_fraction(
    tmp_path: Path,
) -> None:
    """pEVE evidence is ambiguous and is not part of the classified denominator."""
    scorer = _make_scorer(tmp_path)
    best_hits = {
        **{
            f"p_cell_{i}": _blast_hit(
                f"p_cell_{i}", "EUK", "Eustigmatales", identity=90.0
            )
            for i in range(3)
        },
        **{
            f"p_eve_{i}": _blast_hit(f"p_eve_{i}", "EUK-pEVE", identity=90.0)
            for i in range(5)
        },
    }
    features = scorer.compute_blast_features(best_hits, "", "")
    assert features["nonviral_best_hit_fraction"] == 100.0
    assert features["ambiguous_eve_best_hit_count"] == 5
    assert features["ambiguous_eve_best_hit_fraction"] == 62.5

    eve_proteins = [f"p_eve_{i}" for i in range(5)]
    tree = _merge(
        *[
            _tree_nn(f"M_eve_{i}", protein, "EUK-pEVE__Host_A")
            for i, protein in enumerate(eve_proteins)
        ]
    )
    classification = scorer._classify_contig_purity(
        contig_id="contig_eve",
        contig_proteins=eve_proteins,
        tree_nn_results=tree,
        best_hits={protein: best_hits[protein] for protein in eve_proteins},
    )
    assert classification["classification"] == "ambiguous"
    assert classification["n_ambiguous_eve_proteins"] == 5


def test_cellular_coherent_contig_with_consistent_order(tmp_path: Path) -> None:
    """Contig with 4 markers all resolving to EUK__Nannochloropsis passes."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Nannochloropsis_B"),
        _tree_nn("M3", "p3", "EUK__Nannochloropsis_A"),
        _tree_nn("M4", "p4", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=85.0)
        for pid in ("p1", "p2", "p3", "p4")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_cell",
        contig_proteins=["p1", "p2", "p3", "p4"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "cellular_coherent"
    assert result["n_cellular_proteins"] == 4
    # All four proteins agree on (EUK, Eustigmatales).
    assert result["cellular_lineage_purity"] == pytest.approx(1.0)


def test_cellular_scatter_fails_purity(tmp_path: Path) -> None:
    """Four EUK markers scattered across Nannochloropsis / Acanthocystis /
    NA / Other -> purity 0.25, classifier returns ambiguous."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Acanthocystis_A"),
        _tree_nn("M3", "p3", "EUK__NA_A"),
        _tree_nn("M4", "p4", "EUK__Other_A"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "various", identity=85.0)
        for pid in ("p1", "p2", "p3", "p4")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_scatter",
        contig_proteins=["p1", "p2", "p3", "p4"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "ambiguous"
    assert result["cellular_lineage_purity"] < 0.6


def test_below_marker_floor_is_ambiguous(tmp_path: Path) -> None:
    """Two cellular markers fall under the 3-marker floor -> ambiguous."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=90.0)
        for pid in ("p1", "p2")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_short",
        contig_proteins=["p1", "p2"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "ambiguous"


def test_low_identity_rejects_cellular_coherent(tmp_path: Path) -> None:
    """Four markers all tree-resolving to one EUK order but median identity
    45% -> HGT divergent, classifier must refuse the cellular call."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Nannochloropsis_A"),
        _tree_nn("M3", "p3", "EUK__Nannochloropsis_B"),
        _tree_nn("M4", "p4", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=45.0)
        for pid in ("p1", "p2", "p3", "p4")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_divergent",
        contig_proteins=["p1", "p2", "p3", "p4"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "ambiguous"


def test_protein_fraction_aggregates_correctly(tmp_path: Path) -> None:
    """Full bin-level aggregation: 3 contigs, one cellular_coherent, two
    viral_bearing. Check protein fraction over marker-bearing proteins."""
    scorer = _make_scorer(tmp_path)

    protein_to_contig = {
        # cellular contig (4 proteins, all EUK Nannochloropsis)
        "pC1": "contig_cell",
        "pC2": "contig_cell",
        "pC3": "contig_cell",
        "pC4": "contig_cell",
        # viral contig 1 (3 proteins)
        "pV1_a": "contig_viral_1",
        "pV1_b": "contig_viral_1",
        "pV1_c": "contig_viral_1",
        # viral contig 2 (3 proteins)
        "pV2_a": "contig_viral_2",
        "pV2_b": "contig_viral_2",
        "pV2_c": "contig_viral_2",
    }
    contig_lengths = {
        "contig_cell": 12000,
        "contig_viral_1": 30000,
        "contig_viral_2": 20000,
    }
    tree = _merge(
        *[
            _tree_nn(f"M_cell_{i}", pid, "EUK__Nannochloropsis_A")
            for i, pid in enumerate(("pC1", "pC2", "pC3", "pC4"))
        ],
        *[
            _tree_nn(f"M_viral_{i}", pid, "NCLDV__Mimi_1")
            for i, pid in enumerate(
                ("pV1_a", "pV1_b", "pV1_c", "pV2_a", "pV2_b", "pV2_c")
            )
        ],
    )
    best_hits = {
        **{
            pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=85.0)
            for pid in ("pC1", "pC2", "pC3", "pC4")
        },
        **{
            pid: _blast_hit(pid, "NCLDV", "Imitervirales", identity=90.0)
            for pid in (
                "pV1_a",
                "pV1_b",
                "pV1_c",
                "pV2_a",
                "pV2_b",
                "pV2_c",
            )
        },
    }

    out = scorer.compute_contig_purity_features(
        protein_to_contig=protein_to_contig,
        contig_lengths=contig_lengths,
        tree_nn_results=tree,
        best_hits=best_hits,
        attribution_mode="fna_gene_calling",
    ).features
    assert out["per_contig_classifier_active"] is True
    assert out["cellular_coherent_contig_count"] == 1
    assert out["viral_bearing_contig_count"] == 2
    # 4 of 10 marker-bearing proteins on the cellular coherent contig.
    assert out["cellular_coherent_protein_fraction"] == pytest.approx(40.0)


def test_bp_fraction_aggregates_correctly(tmp_path: Path) -> None:
    """Diagnostic bp fraction: 12kb cellular of 62kb total -> ~19.35%."""
    scorer = _make_scorer(tmp_path)
    protein_to_contig = {
        "pC1": "contig_cell",
        "pC2": "contig_cell",
        "pC3": "contig_cell",
        "pV1_a": "contig_viral_1",
        "pV1_b": "contig_viral_1",
    }
    contig_lengths = {"contig_cell": 12000, "contig_viral_1": 50000}
    tree = _merge(
        *[
            _tree_nn(f"M_cell_{i}", pid, "EUK__Nannochloropsis_A")
            for i, pid in enumerate(("pC1", "pC2", "pC3"))
        ],
        _tree_nn("M_viral_0", "pV1_a", "NCLDV__Mimi_1"),
        _tree_nn("M_viral_1", "pV1_b", "NCLDV__Mimi_1"),
    )
    best_hits = {
        **{
            pid: _blast_hit(pid, "EUK", "Eustigmatales", 85.0)
            for pid in ("pC1", "pC2", "pC3")
        },
        **{
            pid: _blast_hit(pid, "NCLDV", "Imitervirales", 90.0)
            for pid in ("pV1_a", "pV1_b")
        },
    }
    out = scorer.compute_contig_purity_features(
        protein_to_contig=protein_to_contig,
        contig_lengths=contig_lengths,
        tree_nn_results=tree,
        best_hits=best_hits,
        attribution_mode="fna_gene_calling",
    ).features
    expected_fraction = round((12000 / 62000) * 100.0, 2)
    assert out["cellular_coherent_bp_fraction"] == pytest.approx(
        expected_fraction, abs=0.01
    )


@pytest.mark.parametrize("prefix", ["", "sample|"])
def test_dna_lengths_are_exact_and_mixed_contigs_are_review_only(
    tmp_path: Path,
    prefix: str,
) -> None:
    """Protein counts cannot inflate DNA fractions or make mixed removal candidates."""
    scorer = _make_scorer(tmp_path)
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(
        f">{prefix}contig_a\n{'A' * 1000}\n>{prefix}contig_b\n{'C' * 1000}\n"
    )

    contig_a_proteins = [f"{prefix}contig_a_{i}" for i in range(1, 4)]
    contig_b_cellular = [f"{prefix}contig_b_{i}" for i in range(1, 6)]
    contig_b_viral = f"{prefix}contig_b_6"
    query_faa = tmp_path / "query.faa"
    query_faa.write_text(
        "".join(
            f">{protein}\n{'M' * length}\n"
            for protein, length in [
                *((protein, 50) for protein in contig_a_proteins),
                *((protein, 100) for protein in contig_b_cellular),
                (contig_b_viral, 100),
            ]
        )
    )

    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    cellular_proteins = [*contig_a_proteins, *contig_b_cellular]
    blast_lines = [
        "\t".join(
            [
                protein,
                "EUK__Nannochloropsis_A",
                "90",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "120",
            ]
        )
        for protein in cellular_proteins
    ]
    blast_lines.append(
        "\t".join(
            [
                contig_b_viral,
                "NCLDV__Mimi_1",
                "90",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "120",
            ]
        )
    )
    (blast_dir / "M1.m8").write_text("\n".join(blast_lines) + "\n")

    tree = _merge(
        *[
            _tree_nn(f"M_cell_{i}", protein, "EUK__Nannochloropsis_A")
            for i, protein in enumerate(cellular_proteins)
        ],
        _tree_nn("M_viral", contig_b_viral, "NCLDV__Mimi_1"),
    )
    features = scorer.collect_contig_features(
        query_fna=query_fna,
        query_faa=query_faa,
        scan=scorer.scan_blast_dir(blast_dir),
        evidence=QueryEvidence({}, {}, {}),
        primary_order="Imitervirales",
        tree_nn_results=tree,
    )

    assert features["contig_length_basis"] == "dna_exact"
    assert features["protein_estimated_contig_count"] == 0
    assert features["suspicious_bp_fraction"] == 50.0
    assert len(features["suspicious_contigs"]) == 1
    candidate = features["suspicious_contigs"][0]
    assert candidate["contig_id"] == f"{prefix}contig_a"
    assert candidate["length_bp"] == 1000
    assert candidate["length_basis"] == "dna_exact"
    assert features["mixed_evidence_contig_count"] == 1
    assert [record["contig_id"] for record in features["review_contigs"]] == [
        f"{prefix}contig_b"
    ]
    assert features["review_contigs"][0]["recommended_action"] == "review"


# ---------------------------------------------------------------------------
# Empty orders, tied neighbours, unassessed contigs and mixed marker groups
# ---------------------------------------------------------------------------


def _cellular_contig(
    scorer, neighbors_by_protein: dict[str, str], *, identity: float = 90.0
) -> dict:
    tree = _merge(
        *[
            _tree_nn(f"M_{protein}", protein, neighbor)
            for protein, neighbor in neighbors_by_protein.items()
        ]
    )
    best_hits = {
        protein: _blast_hit(protein, "EUK", "Eustigmatales", identity=identity)
        for protein in neighbors_by_protein
    }
    return scorer._classify_contig_purity(
        contig_id="contig_1",
        contig_proteins=list(neighbors_by_protein),
        tree_nn_results=tree,
        best_hits=best_hits,
    )


def test_proteins_without_an_order_stay_in_the_denominator_and_never_agree(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)

    result = _cellular_contig(
        scorer,
        {
            "p1": "EUK__Nannochloropsis_A",
            "p2": "EUK__Nannochloropsis_B",
            "p3": "EUK__NoOrder_A",
            "p4": "EUK__NoOrder_A",
            "p5": "EUK__NoOrder_A",
        },
    )

    assert result["n_cellular_proteins"] == 5
    assert result["cellular_lineage_purity"] == pytest.approx(0.4)
    assert result["classification"] == "ambiguous"


def test_proteins_that_all_lack_an_order_have_no_purity(tmp_path: Path) -> None:
    scorer = _make_scorer(tmp_path)

    result = _cellular_contig(scorer, {f"p{i}": "EUK__NoOrder_A" for i in range(4)})

    assert result["cellular_lineage_purity"] == 0.0
    assert result["classification"] == "ambiguous"


def test_tied_neighbours_of_one_order_cast_one_vote(tmp_path: Path) -> None:
    scorer = _make_scorer(tmp_path)
    tree = {
        # Two references of one order at the same distance: one vote.
        "M1": {
            "p1": {
                "EUK__Nannochloropsis_A|g_1": 0.5,
                "EUK__Nannochloropsis_B|g_1": 0.5,
            }
        },
        "M2": {"p1": {"EUK__Acanthocystis_A|g_1": 0.4}},
    }

    lineage = scorer._per_protein_dominant_lineage("p1", tree)

    # One vote each: the lexicographically smallest triplet wins the tie.
    assert lineage == {
        "domain": "EUK",
        "order": "Acanthocystidae",
        "genus": "Acanthocystis",
    }


def test_tied_neighbours_of_two_orders_vote_for_their_shared_lineage(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)
    tree = {
        "M1": {
            "p1": {"EUK__Nannochloropsis_A|g_1": 0.5, "EUK__Acanthocystis_A|g_1": 0.5}
        }
    }

    lineage = scorer._per_protein_dominant_lineage("p1", tree)

    assert lineage == {"domain": "EUK", "order": "", "genus": ""}


def test_tied_neighbours_of_two_domains_give_the_protein_no_domain(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)
    tree = {"M1": {"p1": {"EUK__Nannochloropsis_A|g_1": 0.5, "BAC__Eco_A|g_1": 0.5}}}

    assert scorer._per_protein_dominant_lineage("p1", tree) is None


def test_farther_neighbours_do_not_vote(tmp_path: Path) -> None:
    scorer = _make_scorer(tmp_path)
    tree = {"M1": {"p1": {"EUK__Nannochloropsis_A|g_1": 0.1, "BAC__Eco_A|g_1": 0.9}}}

    lineage = scorer._per_protein_dominant_lineage("p1", tree)

    assert lineage["domain"] == "EUK"


def _m8_line(protein: str, subject: str = "EUK__Nannochloropsis_A") -> str:
    fields = [protein, subject, "90", "0", "0", "0", "0", "0", "0", "0", "0", "120"]
    return "\t".join(fields)


def _contig_features_with_cellular_contig(
    scorer,
    tmp_path: Path,
    cellular_markers: list[str],
    cellular_neighbors: dict[str, float] | None = None,
) -> dict:
    """Return the features of a query whose classifier is active.

    ``contig_v`` holds five tree-placed viral proteins. ``contig_a`` holds one
    protein per cellular marker. Each of them has the tree neighbours
    ``cellular_neighbors`` (neighbour name to distance), or none.
    """
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(f">contig_v\n{'A' * 1000}\n>contig_a\n{'C' * 1000}\n")
    viral = [f"contig_v_{i}" for i in range(1, 6)]
    cellular = [f"contig_a_{i}" for i, _ in enumerate(cellular_markers, start=1)]
    query_faa = tmp_path / "query.faa"
    query_faa.write_text("".join(f">{p}\nMKLV\n" for p in [*viral, *cellular]))
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    (blast_dir / "M_viral.m8").write_text(
        "\n".join(_m8_line(p, "NCLDV__Mimi_1") for p in viral) + "\n"
    )
    for protein, marker in zip(cellular, cellular_markers, strict=True):
        (blast_dir / f"{marker}.m8").write_text(_m8_line(protein) + "\n")
    tree = _merge(*[_tree_nn(f"M_{p}", p, "NCLDV__Mimi_1") for p in viral])
    if cellular_neighbors:
        tree.update({f"M_{p}": {p: cellular_neighbors} for p in cellular})
    return scorer.collect_contig_features(
        query_fna=query_fna,
        query_faa=query_faa,
        scan=scorer.scan_blast_dir(blast_dir),
        evidence=QueryEvidence({}, {}, {}),
        primary_order="Imitervirales",
        tree_nn_results=tree,
    )


def test_unplaced_contig_gets_the_cellular_marker_rule_while_classifier_is_active(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)

    features = _contig_features_with_cellular_contig(
        scorer, tmp_path, [BUSCO_MODELS[0], UNI56_MODELS[0]]
    )

    assert features["per_contig_classifier_active"] is True
    assert [(c["contig_id"], c["reason"]) for c in features["suspicious_contigs"]] == [
        ("contig_a", "cellular_markers")
    ]


def test_unplaced_contig_with_one_cellular_marker_is_not_suspicious(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)

    features = _contig_features_with_cellular_contig(
        scorer, tmp_path, [BUSCO_MODELS[0]]
    )

    assert features["suspicious_contigs"] == []


def test_contig_with_tied_domain_placements_is_assessed_and_not_flagged(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)

    features = _contig_features_with_cellular_contig(
        scorer,
        tmp_path,
        [BUSCO_MODELS[0], UNI56_MODELS[0]],
        cellular_neighbors={"BAC__Eco_A|g_1": 0.5, "NCLDV__Mimi_1|g_1": 0.5},
    )

    assert features["per_contig_classifier_active"] is True
    assert features["suspicious_contigs"] == []
    assert features["suspicious_bp_fraction"] == 0.0


def test_contig_whose_neighbours_have_no_label_still_gets_the_marker_rule(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)

    features = _contig_features_with_cellular_contig(
        scorer,
        tmp_path,
        [BUSCO_MODELS[0], UNI56_MODELS[0]],
        cellular_neighbors={"UNLABELLED__ref|g_1": 0.5},
    )

    assert [(c["contig_id"], c["reason"]) for c in features["suspicious_contigs"]] == [
        ("contig_a", "cellular_markers")
    ]


@pytest.mark.parametrize(
    ("marker", "scores", "category"),
    [
        # grp_COG0085 holds the cellular COG0085 and the viral GVOGm0022.
        ("grp_COG0085", {"COG0085": 500.0, "GVOGm0022": 300.0}, "cellular"),
        ("grp_COG0085", {"COG0085": 300.0, "GVOGm0022": 500.0}, "viral"),
        ("grp_COG0085", {"COG0085": 400.0, "GVOGm0022": 400.0}, "viral"),
        ("grp_COG0085", {"COG0085": 400.0}, "cellular"),
        ("grp_COG0085", {"OTHER": 900.0}, "viral"),
        ("grp_COG0085", {}, "viral"),
        # grp_OG114 holds the viral OG114 and the phage GENOMAD.029205.VV.
        ("grp_OG114", {"OG114": 40.0, "GENOMAD.029205.VV": 90.0}, "phage"),
        ("grp_OG114", {"OG114": 90.0, "GENOMAD.029205.VV": 40.0}, "viral"),
        # A group whose members share one category keeps it.
        ("grp_1275837at2759", {}, "cellular"),
        ("not_a_group", {"COG0085": 500.0}, "viral"),
    ],
)
def test_marker_category_of_a_mixed_group_follows_the_best_scoring_model(
    marker: str, scores: dict[str, float], category: str
) -> None:
    assert _marker_category(marker, scores) == category


def _two_group_contig_features(
    scorer, tmp_path: Path, scores: dict[str, dict[str, float]]
) -> dict:
    """Return the features of one contig whose two proteins hit mixed groups."""
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(f">contig_a\n{'C' * 1000}\n")
    proteins = list(scores)
    query_faa = tmp_path / "query.faa"
    query_faa.write_text("".join(f">{p}\nMKLV\n" for p in proteins))
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    for protein, group in zip(proteins, ["grp_COG0085", "grp_COG0086"], strict=True):
        (blast_dir / f"{group}.m8").write_text(_m8_line(protein) + "\n")
    return scorer.collect_contig_features(
        query_fna=query_fna,
        query_faa=query_faa,
        scan=scorer.scan_blast_dir(blast_dir),
        evidence=QueryEvidence({}, {}, scores),
        primary_order="Imitervirales",
        tree_nn_results={},
    )


def test_cellular_models_in_mixed_groups_count_as_cellular_markers(
    tmp_path: Path,
) -> None:
    scorer = _make_scorer(tmp_path)

    features = _two_group_contig_features(
        scorer,
        tmp_path,
        {
            "contig_a_1": {"COG0085": 500.0, "GVOGm0022": 300.0},
            "contig_a_2": {"COG0086": 500.0, "GVOGm0023": 300.0},
        },
    )

    candidate = features["suspicious_contigs"][0]
    assert candidate["reason"] == "cellular_markers"
    assert candidate["cellular_marker_count"] == 2
    assert candidate["viral_marker_count"] == 0


def test_viral_models_in_mixed_groups_count_as_viral_markers(tmp_path: Path) -> None:
    scorer = _make_scorer(tmp_path)

    features = _two_group_contig_features(
        scorer,
        tmp_path,
        {
            "contig_a_1": {"COG0085": 300.0, "GVOGm0022": 500.0},
            "contig_a_2": {"COG0086": 300.0, "GVOGm0023": 500.0},
        },
    )

    assert features["suspicious_contigs"] == []
