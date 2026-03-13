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
    assert result["mcp_total"] == 5
    assert result["vp_completeness"] == "3/4"
    assert result["vp_mcp"] == 2
    assert result["plv"] == 2
    assert result["vp_df"] == 1.0
    assert result["mirus_completeness"] == "2/4"
    assert result["mirus_df"] == 0.5
    assert result["mrya_total"] == 2
    assert result["phage_total"] == 2
    assert result["cellular_total"] == 2
    assert result["cellular_dup"] == 1.0


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
