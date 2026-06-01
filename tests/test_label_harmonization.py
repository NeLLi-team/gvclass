"""Focused tests for reusable resource-label harmonization rules."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "harmonize_labels.py"
SPEC = importlib.util.spec_from_file_location("harmonize_labels", SCRIPT)
assert SPEC is not None
harmonize_labels = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["harmonize_labels"] = harmonize_labels
SPEC.loader.exec_module(harmonize_labels)


def _record(label_id: str, taxonomy: str):
    return harmonize_labels.LabelRecord(
        label_id=label_id,
        ranks=tuple(taxonomy.split("|")),
        line_no=1,
    )


def _evidence() -> harmonize_labels.Evidence:
    evidence = harmonize_labels.Evidence()
    evidence.add("marker.faa")
    return evidence


def test_candidate_label_ids_include_current_faa_resolution_forms() -> None:
    candidates = harmonize_labels.candidate_label_ids(
        "VP__IMGVR_UViG_3300050337|002562_1"
    )

    assert "VP__IMGVR_UViG_3300050337" in candidates
    assert "VP__IMGVR_UViG_3300050337|002562" in candidates
    assert "VP__IMGVR_UViG_3300050337_002562" in candidates


def test_mirus_placeholder_is_lifted_to_mirusviricota_with_domain_placeholders() -> None:
    ranks, contexts = harmonize_labels.harmonize_taxonomy(
        "MIRUS__TARA_PSW_NCLDV_00055",
        tuple("MIRUS|MIRUS|MIRUS|MIRUS|MIRUS|MIRUS|TARA_PSW_NCLDV_00055".split("|")),
        {},
    )

    assert contexts == []
    assert ranks == (
        "MIRUS",
        "Mirusviricota",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "TARA_PSW_NCLDV_00055",
    )


def test_mirus_mcp_host_context_is_removed_from_taxonomy() -> None:
    ranks, contexts = harmonize_labels.harmonize_taxonomy(
        "MIRUS__EP00178",
        tuple(
            "MIRUS|MCP in Eurhodophytina|MCP in Florideophyceae|"
            "MCP in Gracilariales|MCP in Rhodymeniophycidae|"
            "MCP in Gracilaria|MCP in Gracilaria changii".split("|")
        ),
        {},
    )

    assert ranks == (
        "MIRUS",
        "Mirusviricota",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "EP00178",
    )
    assert ("phylum", "MCP in Eurhodophytina") in contexts
    assert ("class", "MCP in Florideophyceae") in contexts
    assert ("species", "MCP in Gracilaria changii") in contexts


def test_duplicate_metadata_label_aliases_to_unique_faa_backed_canonical() -> None:
    taxonomy = (
        "MIRUS|Mirusviricota|Class_01|Cryptic_ORD10|"
        "MIRUS_unclassified|MIRUS_unclassified|MIRUS__Delmont2025__MIRUS_G_0050"
    )
    labels = {
        "MIRUS__Delmont2025_MIRUS_G_0050": _record(
            "MIRUS__Delmont2025_MIRUS_G_0050", taxonomy
        ),
        "MIRUS__MIRUS_G_0050": _record("MIRUS__MIRUS_G_0050", taxonomy),
    }
    evidence = {"MIRUS__MIRUS_G_0050": _evidence()}

    active, aliases, _contexts, inactive, _excluded = harmonize_labels.harmonize_bundle(
        labels, evidence, {}
    )

    assert "MIRUS__MIRUS_G_0050" in active
    assert aliases["MIRUS__Delmont2025_MIRUS_G_0050"] == (
        "MIRUS__MIRUS_G_0050",
        "duplicate_taxonomy_without_faa",
    )
    assert inactive == []


def test_phage_virophage_marker_label_becomes_canonical_vp_with_aliases() -> None:
    labels = {
        "PHAGE__GCA-000915035-1": _record(
            "PHAGE__GCA-000915035-1",
            "PHAGE|PHAGE__Zamilon virus|PHAGE__Zamilon virus|"
            "PHAGE__Zamilon virus|PHAGE__Zamilon virus|PHAGE__Zamilon virus|"
            "PHAGE__Zamilon virus",
        ),
        "PLV__GCA-000915035-1": _record(
            "PLV__GCA-000915035-1",
            "PLV|Preplasmiviricota|Maveriviricetes|Priklausovirales|"
            "Lavidaviridae|Mimivirus-dependent virus Zamilon|Zamilon virus",
        ),
    }
    evidence = {"PHAGE__GCA-000915035-1": _evidence()}

    active, aliases, contexts, inactive, _excluded = harmonize_labels.harmonize_bundle(
        labels, evidence, {}
    )

    assert active["VP__GCA-000915035-1"] == (
        "VP",
        "Preplasmiviricota",
        "Maveriviricetes",
        "Priklausovirales",
        "Lavidaviridae",
        "Mimivirus-dependent virus Zamilon",
        "Zamilon virus",
    )
    assert aliases["PHAGE__GCA-000915035-1"] == (
        "VP__GCA-000915035-1",
        "phage_virophage_relabel",
    )
    assert aliases["PLV__GCA-000915035-1"] == (
        "VP__GCA-000915035-1",
        "plv_virophage_alias",
    )
    assert {row["context_type"] for row in contexts} == {
        "original_phage_taxonomy",
        "original_plv_taxonomy",
    }
    assert inactive == []


def test_rewrite_header_token_preserves_gene_suffix_and_payload() -> None:
    alias_targets = {"PHAGE__GCA-000915035-1": "VP__GCA-000915035-1"}

    assert harmonize_labels.rewrite_header_token(
        "PHAGE__GCA-000915035-1_7", alias_targets
    ) == (
        "VP__GCA-000915035-1_7",
        "PHAGE__GCA-000915035-1",
        "VP__GCA-000915035-1",
    )
    assert harmonize_labels.rewrite_header_token(
        "PHAGE__GCA-000915035-1|contig_7", alias_targets
    ) == (
        "VP__GCA-000915035-1|contig_7",
        "PHAGE__GCA-000915035-1",
        "VP__GCA-000915035-1",
    )


def test_rewrite_header_token_ignores_alias_that_already_resolves_by_base_id() -> None:
    alias_targets = {
        "VP__IMGVR_UViG_3300037332_000371": "VP__IMGVR_UViG_3300037332"
    }

    assert harmonize_labels.rewrite_header_token(
        "VP__IMGVR_UViG_3300037332|000371_3", alias_targets
    ) == (
        "VP__IMGVR_UViG_3300037332|000371_3",
        None,
        None,
    )


def test_prefix_colliding_active_labels_are_not_mangled_without_alias() -> None:
    active = {"VP__Sputnik_virophage", "VP__Sputnik_virophage_3"}

    assert harmonize_labels.rewrite_ambiguous_header_token(
        "VP__Sputnik_virophage|3_2", active
    ) == (
        "VP__Sputnik_virophage|3_2",
        None,
        None,
    )
    assert harmonize_labels.rewrite_ambiguous_header_token(
        "VP__Sputnik_virophage_3", active
    ) == (
        "VP__Sputnik_virophage_3",
        None,
        None,
    )


def test_euk_taxonkit_lineage_replaces_rank_shifted_context() -> None:
    lineage = {
        "Rhincodon typus": {
            "phylum": "Chordata",
            "class": "Chondrichthyes",
            "order": "Orectolobiformes",
            "family": "Rhincodontidae",
            "genus": "Rhincodon",
            "species": "Rhincodon typus",
            "taxid": "259920",
            "terminal_rank": "species",
        }
    }

    ranks, contexts = harmonize_labels.harmonize_taxonomy(
        "EUK__GCA-001642345-1",
        tuple(
            "EUK|Metazoa|Vertebrata|unclassified|Rhincodontidae|"
            "Rhincodon|Rhincodon typus".split("|")
        ),
        lineage,
    )

    assert ranks == (
        "EUK",
        "Chordata",
        "Chondrichthyes",
        "Orectolobiformes",
        "Rhincodontidae",
        "Rhincodon",
        "Rhincodon typus",
    )
    assert ("ncbi_taxid", "259920") in contexts
    assert any(context_type == "original_taxonomy" for context_type, _ in contexts)
