"""Marker panels and marker groups must name models that the database holds.

A panel name that the database does not hold never matches a hit. The panel then
reports a wrong size and a constant completeness, and no test sees it. These tests
compare every panel constant, every prefix panel and every marker group member with
``docs/reference/marker-annotations.tsv``, which lists the database models.
"""

import csv
from pathlib import Path

import pytest

from src.config import marker_sets
from src.config.marker_sets import (
    BUSCO_MODELS,
    PANEL_MODELS,
    PANEL_PREFIXES,
    PHAGE_MODELS,
    load_marker_groups,
    missing_panel_models,
)
from src.core import contamination_scoring
from src.core.summarize import marker_panels
from src.core.summarize.paths_stats import QueryEvidence

REPO_ROOT = Path(__file__).resolve().parents[1]
ANNOTATIONS_TSV = REPO_ROOT / "docs" / "reference" / "marker-annotations.tsv"
MARKER_GROUPS_TSV = REPO_ROOT / "src" / "config" / "marker_groups.tsv"

# Group members that are in the database HMM file but not in the annotation table:
# the order-level OG models without an annotation, the PLV_MCP_* models and the
# capscan models. A new name outside this set is drift and fails the test.
GROUP_MEMBERS_WITHOUT_ANNOTATION = frozenset(
    """
    OG10 OG1046 OG114 OG12902 OG13 OG130 OG1357 OG1419 OG1421 OG14872 OG15 OG1515
    OG15424 OG1552 OG1590 OG1608 OG1683 OG1688 OG17 OG17782 OG1859 OG19209 OG1983
    OG2081 OG2155 OG21616 OG2196 OG23 OG249 OG25874 OG25875 OG25876 OG26 OG2605
    OG267 OG27 OG270 OG28775 OG299 OG30 OG33 OG34 OG3465 OG3587 OG36 OG37 OG3859
    OG4181 OG4547 OG4640 OG516 OG5703 OG5880 OG6034 OG6097 OG6330 OG6520 OG6561
    OG667 OG771 OG806 OG8381 OG876 OG9378
    PLV_MCP_1 PLV_MCP_10 PLV_MCP_2 PLV_MCP_3 PLV_MCP_4 PLV_MCP_5 PLV_MCP_6
    PLV_MCP_8 PLV_MCP_9
    plv_mcp_caps_Dino_Aquinto plv_mcp_caps_GKS1_Gossevirus_Mavpol
    plv_mcp_caps_GKS2_Alpenseevirus_NCV_like plv_mcp_caps_GKS3_CHI_Aquinto
    plv_mcp_caps_Mavpol2 plv_mcp_caps_Nova6_VC35_Mavpol2
    plv_mcp_caps_Nova9_Aquinto_Trimcap_related plv_mcp_caps_PaM2_NCV_like
    plv_mcp_caps_PgVV_Aquinto plv_mcp_caps_SP67_Aquinto plv_mcp_caps_SP70_Aquinto
    plv_mcp_caps_SP71_Aquinto plv_mcp_caps_SP_Aquinto
    plv_mcp_caps_Trimcap_A1_Aquinto plv_mcp_caps_Trimcap_A2_Aquinto
    plv_mcp_caps_Trimcap_A3_74_Aquinto plv_mcp_caps_Trimcap_A3_Aquinto
    plv_mcp_caps_Trimcap_A4_Aquinto plv_mcp_caps_Trimcap_B1_Aquinto
    plv_mcp_caps_Trimcap_B2_Aquinto plv_mcp_caps_Trimcap_B3_Aquinto
    plv_mcp_caps_VC104_Aquinto plv_mcp_caps_VC124_VC134_NCV_like
    plv_mcp_caps_VC147_GKS3_Aquinto plv_mcp_caps_VC224_NA
    plv_mcp_caps_VC40_NCV_like plv_mcp_caps_Zephyrvirus4_NCV_like
    plv_mcp_caps_Zephyrvirus5_NCV_like plv_mcp_caps_Zoravirus_VC102_NCV_like
    """.split()  # noqa: SIM905  # 100 names read better as wrapped words
)

MODEL_LIST_CONSTANTS = sorted(
    name
    for name, value in vars(marker_sets).items()
    if name.endswith("_MODELS") and isinstance(value, list)
)


def _database_models() -> set[str]:
    with ANNOTATIONS_TSV.open(encoding="utf-8", newline="") as handle:
        return {row["marker_id"] for row in csv.DictReader(handle, delimiter="\t")}


def _marker_group_members() -> set[str]:
    members: set[str] = set()
    with MARKER_GROUPS_TSV.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("#", "group_name\t")):
                continue
            members.update(m for m in line.rstrip("\n").split("\t")[7].split(";") if m)
    return members


@pytest.mark.parametrize("constant", MODEL_LIST_CONSTANTS)
def test_model_list_constant_names_only_database_models(constant):
    models = getattr(marker_sets, constant)

    absent = sorted(set(models) - _database_models())

    assert absent == []


def test_model_list_constants_include_the_panels():
    assert {"BUSCO_MODELS", "PHAGE_MODELS", "UNI56_MODELS", "MCP_MODELS"} <= set(
        MODEL_LIST_CONSTANTS
    )


@pytest.mark.parametrize("panel", sorted(PANEL_PREFIXES))
def test_prefix_panel_matches_a_database_model(panel):
    prefix = PANEL_PREFIXES[panel]

    matches = [m for m in _database_models() if m.startswith(prefix)]

    assert matches != []


def test_missing_panel_models_is_empty_for_the_database_models():
    assert missing_panel_models(_database_models()) == {}


def test_missing_panel_models_lists_absent_models_per_panel():
    dropped = {PHAGE_MODELS[0], BUSCO_MODELS[0]}
    available = {
        m
        for m in _database_models() - dropped
        if not m.startswith(PANEL_PREFIXES["vp_protease"])
    }

    missing = missing_panel_models(available)

    assert missing == {
        "phage": [PHAGE_MODELS[0]],
        "busco": [BUSCO_MODELS[0]],
        "vp_protease": ["VP_PRO*"],
    }


def test_missing_panel_models_reports_every_panel_when_nothing_is_available():
    missing = missing_panel_models(set())

    assert set(missing) == set(PANEL_MODELS) | set(PANEL_PREFIXES)
    assert missing["busco"] == BUSCO_MODELS


def test_busco_panel_holds_the_models_the_database_ships():
    assert len(BUSCO_MODELS) == 69
    assert len(set(BUSCO_MODELS)) == 69


def test_phage_panel_uses_database_names():
    assert PHAGE_MODELS[:2] == ["GENOMAD.000023.VV", "GENOMAD.000026.VV"]
    assert len(PHAGE_MODELS) == 20


def test_marker_group_members_are_database_models():
    unknown = _marker_group_members() - _database_models()

    assert unknown == GROUP_MEMBERS_WITHOUT_ANNOTATION


def test_add_marker_metrics_counts_database_phage_models():
    marker_counts = {"GENOMAD.000023.VV": 2, "GENOMAD.000026.VV": 1}

    result = marker_panels._add_marker_metrics(
        QueryEvidence(marker_counts, {}, {}), None, 0, ""
    )

    assert result["phage_completeness"] == "2/20"
    assert result["phage_unique"] == 2
    assert result["phage_total"] == 3


def test_group_inherits_a_category_only_when_all_members_are_in_it(monkeypatch):
    groups = {"grp_all": ["a", "b"], "grp_mixed": ["a", "x"], "grp_none": ["x", "y"]}
    monkeypatch.setattr(contamination_scoring, "GROUP_TO_MODELS", groups)

    names = contamination_scoring._category_group_names({"a", "b"})

    assert names == {"grp_all"}


def test_phage_group_set_excludes_the_mixed_group():
    assert "grp_GENOMAD.000023.VV" in contamination_scoring.PHAGE_MODELS_SET
    assert "grp_OG114" not in contamination_scoring.PHAGE_MODELS_SET


def test_cellular_group_set_excludes_groups_that_hold_a_viral_member():
    cellular = contamination_scoring.CELLULAR_MODELS

    assert "grp_1275837at2759" in cellular
    assert "grp_COG0085" not in cellular
    assert "grp_COG0086" not in cellular


def test_load_marker_groups_raises_for_a_missing_table(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_marker_groups(tmp_path / "missing.tsv")


def test_load_marker_groups_raises_for_a_short_row(tmp_path):
    table = tmp_path / "groups.tsv"
    table.write_text("group_name\tn_members\ngrp_a\t2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="line 2"):
        load_marker_groups(table)


def test_load_marker_groups_raises_for_a_table_without_groups(tmp_path):
    table = tmp_path / "groups.tsv"
    table.write_text("# comment only\n", encoding="utf-8")

    with pytest.raises(ValueError, match="no marker groups"):
        load_marker_groups(table)


def test_load_marker_groups_reads_combinable_groups(tmp_path):
    table = tmp_path / "groups.tsv"
    table.write_text(
        "group_name\tn_members\ttotal_refs\tannotation_coverage\tcoherence\t"
        "dom_annotation\tflag\tmembers\n"
        "grp_a\t2\t10\t1.00\t1.00\tx\tok\tm1;m2\n"
        "grp_b\t2\t10\t1.00\t0.20\tx\tREVIEW_mixed\tm3;m4\n",
        encoding="utf-8",
    )

    group_to_models, model_to_group = load_marker_groups(table)

    assert group_to_models == {"grp_a": ["m1", "m2"]}
    assert model_to_group == {"m1": "grp_a", "m2": "grp_a"}
