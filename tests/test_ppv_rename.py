"""Regression tests for the PPV (Preplasmiviricota) domain rename.

Pins two contracts that must not silently drift:
  1. the contamination scorer recognises the PPV domain as phage-like (a skew here
     would silently weaken contamination scoring for every Preplasmiviricota ref);
  2. the migration transform converts exactly the domain prefix / leading token and
     never touches marker names or single-underscore lifestyle tags.
"""
import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


def _load_migration():
    path = REPO / "scripts" / "relabel_ppv_preplasmiviricota.py"
    spec = importlib.util.spec_from_file_location("ppv_migration", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_ppv_is_phage_domain_not_cellular_or_giant():
    from src.core.contamination_scoring import (
        PHAGE_PREFIXES,
        CELLULAR_PREFIXES,
        GIANT_VIRUS_PREFIXES,
    )

    assert "PPV" in PHAGE_PREFIXES
    assert "PPV" not in CELLULAR_PREFIXES
    assert "PPV" not in GIANT_VIRUS_PREFIXES
    # legacy codes kept transitionally so a pre-PPV bundle still scores correctly
    assert {"PLV", "VP"} <= PHAGE_PREFIXES


@pytest.mark.parametrize(
    "src,expected",
    [
        ("PLV__EP00001", "PPV__EP00001"),
        ("VP__Dishui_lake_virophage__1", "PPV__Dishui_lake_virophage__1"),
        ("PLV__GCA-015143345-1__vpdup", "PPV__GCA-015143345-1__vpdup"),
        ("PLV", "PPV"),
        ("VP", "PPV"),
    ],
)
def test_transform_changes_domain_ids(src, expected):
    assert _load_migration().remap_token(src) == expected


@pytest.mark.parametrize(
    "tok",
    [
        "PLV_MCP_1", "VP_MCP_2", "VP_ATPase_3", "VP_Penton_1", "VP_PRO_1",
        "plv_mcp_caps_PgVV_Aquinto", "gamadvirusMCP",
        "PLV_unclassified", "VP_unclassified",
        "Polintoviricetes", "Virophaviricetes", "Preplasmiviricota",
        "NCLDV__GCA_x", "PHAGE__x", "BAC",
    ],
)
def test_transform_preserves_markers_and_lifestyle_tags(tok):
    assert _load_migration().remap_token(tok) == tok


def test_lineage_preserves_lifestyle_but_changes_domain_and_pos7():
    m = _load_migration()
    src = ("PLV|Preplasmiviricota|Virophaviricetes|VP_unclassified|"
           "VP_unclassified|VP_unclassified|VP__Dishui_lake_virophage__1")
    expected = ("PPV|Preplasmiviricota|Virophaviricetes|VP_unclassified|"
                "VP_unclassified|VP_unclassified|PPV__Dishui_lake_virophage__1")
    assert m.remap_lineage(src) == expected


def test_faa_header_protein_suffix_untouched():
    m = _load_migration()
    assert m.remap_header("PLV__EMALE07_E4-10__001|prot_PLV_keepme") == \
        "PPV__EMALE07_E4-10__001|prot_PLV_keepme"
