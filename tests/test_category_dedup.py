from collections import Counter
from pathlib import Path

from src.config.marker_sets import BUSCO_MODELS, PHAGE_MODELS, UNI56_MODELS
from src.core.contamination_scoring import ContaminationScorer
from src.core.summarize_full import FullSummarizer


def test_add_marker_metrics_deduplicates_proteins_within_categories():
    summarizer = FullSummarizer.__new__(FullSummarizer)
    marker_hits = {
        "OG1352": {"ncldv_p1"},
        "OG484": {"ncldv_p1", "ncldv_p2", "ncldv_p4"},
        "GVOGm0003": {"ncldv_p1", "ncldv_p3"},
        "VP_MCP_alpha": {"vp_p1", "vp_p2"},
        "VP_MCP_beta": {"vp_p1"},
        "VP_Penton_alpha": {"vp_p2"},
        "VP_ATPase_alpha": {"vp_p3"},
        "PLV_marker_a": {"plv_p1"},
        "PLV_marker_b": {"plv_p1", "plv_p2"},
        "Mirus_MCP": {"mirus_p1"},
        "Mirus_JellyRoll": {"mirus_p1"},
        "Mirus_Portal": {"mirus_p2"},
        "VLTF3": {"mrya_p1"},
        "ATPase": {"mrya_p1", "mrya_p2"},
        PHAGE_MODELS[0]: {"phage_p1"},
        PHAGE_MODELS[1]: {"phage_p1", "phage_p2"},
        BUSCO_MODELS[0]: {"cell_p1"},
        UNI56_MODELS[0]: {"cell_p1", "cell_p2"},
    }
    marker_counts = {marker: len(proteins) for marker, proteins in marker_hits.items()}

    result = {}
    summarizer._add_marker_metrics(result, marker_counts, marker_hits)

    assert result["ncldv_mcp_total"] == 4
    assert result["vp_completeness"] == "3/4"
    assert result["vp_mcp"] == 2
    assert result["plv"] == 2
    assert result["vp_df"] == 1.0
    assert result["mirus_completeness"] == "2/4"
    assert result["mirus_df"] == 0.5
    assert result["mrya_completeness"] == "2/6"
    assert result["mrya_dup"] == 1.5
    assert result["phage_completeness"] == "2/20"
    assert result["phage_dup"] == 1.5
    assert result["busco_completeness"] == "1/255"
    assert result["busco_dup"] == 1.0
    assert result["cog_completeness"] == "1/56"
    assert result["cog_dup"] == 2.0


def test_mrya_consolidation_excludes_grouped_nonmrya_winners():
    """A grouped MRYA model (VLTF3/ATPase/...) co-occurs with GVOG/NCLDV markers
    on ordinary NCLDV proteins. It must only count toward MRYA when the MRYA model
    is the protein's best hit within its marker group; standalone MRYA models
    (HUH) always count."""
    summarizer = FullSummarizer.__new__(FullSummarizer)
    marker_scores = {
        # VLTF3 (MRYA) loses to co-grouped GVOGm0890 (GVOG8) -> not MRYA.
        "p1": {"VLTF3": 100.0, "GVOGm0890": 200.0},
        # ATPase (MRYA) beats co-grouped GVOGm0760 (GVOG8) -> genuine MRYA.
        "p2": {"ATPase": 300.0, "GVOGm0760": 250.0},
        # HUH is a standalone MRYA marker -> always counts.
        "p3": {"HUH": 50.0},
    }
    marker_counts = {m: 1 for scores in marker_scores.values() for m in scores}

    result = {}
    summarizer._add_marker_metrics(
        result, marker_counts, marker_hits=None, marker_scores=marker_scores
    )

    assert result["mrya_completeness"] == "2/6"  # ATPase + HUH, not VLTF3
    assert result["mrya_dup"] == 1.0


def test_mrya_falls_back_to_raw_panel_without_scores():
    """Counts-only path (no per-protein scores) keeps the raw panel behaviour."""
    summarizer = FullSummarizer.__new__(FullSummarizer)
    marker_counts = {"VLTF3": 1, "ATPase": 2, "GVOGm0890": 1, "GVOGm0760": 1}

    result = {}
    summarizer._add_marker_metrics(result, marker_counts, marker_hits=None)

    assert result["mrya_completeness"] == "2/6"  # VLTF3 + ATPase, no consolidation


def test_cog_consolidation_excludes_grouped_rna_pol():
    """COG0085/COG0086 (RNA pol) share marker groups with the NCLDV RNA pol
    markers GVOGm0022/GVOGm0023; an NCLDV protein must not inflate the cellular
    cog panel. Standalone COGs still count."""
    summarizer = FullSummarizer.__new__(FullSummarizer)
    marker_scores = {
        # NCLDV RNA pol: GVOGm0022 out-scores co-grouped COG0085 -> not cog.
        "p1": {"COG0085": 500.0, "GVOGm0022": 1500.0},
        # standalone COG -> counts.
        "p2": {"COG0013": 100.0},
    }
    marker_counts = {m: 1 for scores in marker_scores.values() for m in scores}

    result = {}
    summarizer._add_marker_metrics(
        result, marker_counts, marker_hits=None, marker_scores=marker_scores
    )

    assert result["cog_completeness"] == "1/56"  # COG0013 only, COG0085 dropped


def test_vp_completeness_is_tree_bucketed():
    """A VP marker on a protein that places with NCLDV is not counted toward the
    PPV panel; a VP marker placing with PPV is."""
    summarizer = FullSummarizer.__new__(FullSummarizer)
    marker_hits = {"VP_ATPase_1": {"p_ncldv"}, "VP_MCP_a": {"p_ppv"}}
    protein_domain = {"p_ncldv": "NCLDV", "p_ppv": "PPV"}

    comp, _vp_mcp, _plv, _vp_df = summarizer.calculate_vp_metrics(
        {}, marker_hits, tree_nn_results=None, protein_domain=protein_domain
    )
    assert comp == "1/4"  # MCP (PPV) kept, ATPase (NCLDV) dropped

    # Without tree placements, fall back to raw category presence.
    comp_raw, _, _, _ = summarizer.calculate_vp_metrics({}, marker_hits)
    assert comp_raw == "2/4"


def test_mirus_completeness_is_tree_bucketed():
    """A Mirus marker on a protein placing with NCLDV is excluded; one placing
    with MIRUS is kept."""
    summarizer = FullSummarizer.__new__(FullSummarizer)
    marker_hits = {"Mirus_Terminase_ATPase": {"p_ncldv"}, "Mirus_Portal": {"p_mirus"}}
    protein_domain = {"p_ncldv": "NCLDV", "p_mirus": "MIRUS"}

    comp, _df = summarizer.calculate_mirus_completeness(
        {}, marker_hits, protein_domain=protein_domain
    )
    assert comp == "1/4"  # Portal (MIRUS) kept, ATPase (NCLDV) dropped


def test_protein_nearest_domain_min_distance_and_label_lookup():
    """Each protein is bucketed to the domain of its single nearest tree neighbour
    across all trees; neighbours absent from labels are skipped."""
    summarizer = FullSummarizer.__new__(FullSummarizer)
    summarizer.labels_dict = {"NCLDV__g1": ["NCLDV"], "PPV__g2": ["PPV"]}
    tree_nn = {
        # nearest neighbour (0.05) is unlabeled -> skipped; PPV at 0.20 then wins
        # over the NCLDV placement at 0.40 in another tree.
        "g_mcp": {"p1": {"UNKNOWN__x": 0.05, "PPV__g2": 0.20}},
        "g_atp": {"p1": {"NCLDV__g1": 0.40}},
    }
    domains = summarizer._protein_nearest_domain(tree_nn)
    assert domains["p1"] == "PPV"

    # No tree placements -> empty mapping.
    assert summarizer._protein_nearest_domain({}) == {}


def test_rule_based_contamination_uses_unique_proteins_for_category_totals(tmp_path):
    scorer = ContaminationScorer.__new__(ContaminationScorer)
    scorer.collect_best_blast_hits = lambda blast_dir: {}
    scorer.compute_blast_features = lambda best_hits, primary_order, primary_family: {
        "strong_cellular_hit_count": 0,
        "nonviral_best_hit_fraction": 0.0,
        "strong_phage_hit_count": 0,
        "foreign_viral_order_hit_count": 0,
        "foreign_viral_family_hit_count": 0,
        "dominant_nonviral_lineage_fraction": 0.0,
        "dominant_foreign_viral_order_fraction": 0.0,
    }
    scorer.compute_tree_features = lambda tax_counters: {
        "order_majority_fraction": 0.0,
        "family_majority_fraction": 0.0,
        "order_secondary_fraction": 0.0,
        "family_secondary_fraction": 0.0,
    }
    scorer._extract_tax_token = lambda counter: ""

    query_output_dir = Path(tmp_path)
    hmmout_dir = query_output_dir / "hmmout"
    hmmout_dir.mkdir()
    (hmmout_dir / "models.out.filtered").write_text(
        "\n".join(
            [
                f"cell_p1\t-\t-\t{BUSCO_MODELS[0]}",
                f"cell_p1\t-\t-\t{UNI56_MODELS[0]}",
                f"cell_p2\t-\t-\t{UNI56_MODELS[0]}",
                f"phage_p1\t-\t-\t{PHAGE_MODELS[0]}",
                f"phage_p1\t-\t-\t{PHAGE_MODELS[1]}",
                f"phage_p2\t-\t-\t{PHAGE_MODELS[1]}",
            ]
        )
        + "\n"
    )

    marker_counts = {
        BUSCO_MODELS[0]: 1,
        UNI56_MODELS[0]: 2,
        PHAGE_MODELS[0]: 1,
        PHAGE_MODELS[1]: 2,
    }

    result = scorer.score_rule_based(
        result={"order_dup": 1.0, "gvog8_dup": 1.0},
        marker_counts=marker_counts,
        tax_counters={"order": Counter(), "family": Counter()},
        query_output_dir=query_output_dir,
    )

    assert result["contamination_cellular_signal_v1"] == 9.0
    assert result["contamination_phage_signal_v1"] == 10.0
    assert result["contamination_score_v1"] == 4.62
