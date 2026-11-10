"""Marker sets configuration for the GVClass pipeline.

Based on docs/markers_gvclassV1.tsv

Every model name below must exist in the reference database;
``missing_panel_models`` lists the names a database lacks.

Marker Categories:
------------------
NCLDV (Nucleocytoviricota):
  - MCP: gamadvirusMCP, GVOGm0003, OG1352, OG484

VP (Virophage) - 4-gene core for completeness:
  - MCP: VP_MCP_* (prefix match)
  - Penton: VP_Penton_* (prefix match)
  - ATPase: VP_ATPase_* (prefix match)
  - Protease: VP_PRO_* (prefix match)

PLV (Polinton-like virus):
  - Single marker with PLV_ prefix
  - PLVs share VP markers but have PLV marker in addition

Mirus (Mirusviricota) - 4-category core for completeness:
  - MCP: Mirus_MCP, Mirus_JellyRoll
  - ATPase: Mirus_Terminase_ATPase, Mirus_Terminase_merged
  - Portal: Mirus_Portal
  - Triplex: Mirus_Triplex1, Mirus_Triplex2
"""

from pathlib import Path

# =============================================================================
# VP (Virophage) marker categories - prefix-based matching
# =============================================================================
# VP completeness = unique categories present / 4
VP_CATEGORY_PREFIXES = {
    "MCP": "VP_MCP",
    "Penton": "VP_Penton",
    "ATPase": "VP_ATPase",
    "Protease": "VP_PRO",
}

# PLV (Polinton-like virus) marker prefix
# PLVs share VP markers but have PLV-specific marker(s) in addition
PLV_PREFIX = "PLV_"

# =============================================================================
# Mirus (Mirusviricota) marker categories
# =============================================================================
# Mirus completeness = unique categories present / 4
MIRUS_CATEGORY_MODELS = {
    "MCP": ["Mirus_MCP", "Mirus_JellyRoll"],
    "ATPase": ["Mirus_Terminase_ATPase", "Mirus_Terminase_merged"],
    "Portal": ["Mirus_Portal"],
    "Triplex": ["Mirus_Triplex1", "Mirus_Triplex2"],
}

# All Mirus models (flat list)
MIRUS_MODELS = [model for models in MIRUS_CATEGORY_MODELS.values() for model in models]

# BUSCO markers: the 69 BUSCO models that the database ships
BUSCO_MODELS = [
    "1003258at2759",
    "1019762at2759",
    "1032689at2759",
    "1038775at2759",
    "1057950at2759",
    "1079130at2759",
    "1079827at2759",
    "1090038at2759",
    "1106766at2759",
    "1111142at2759",
    "1157302at2759",
    "1161199at2759",
    "1178688at2759",
    "1193442at2759",
    "1194691at2759",
    "1197019at2759",
    "1200489at2759",
    "1220881at2759",
    "1222562at2759",
    "1236198at2759",
    "1249159at2759",
    "1260807at2759",
    "1275837at2759",
    "1284731at2759",
    "1291729at2759",
    "1309031at2759",
    "1322299at2759",
    "1338131at2759",
    "1346165at2759",
    "1346432at2759",
    "1348942at2759",
    "1380409at2759",
    "1382842at2759",
    "1423847at2759",
    "1430056at2759",
    "1432328at2759",
    "1434266at2759",
    "1455730at2759",
    "1488436at2759",
    "1504863at2759",
    "1530008at2759",
    "1552706at2759",
    "1567796at2759",
    "1576404at2759",
    "1588798at2759",
    "1617752at2759",
    "1623701at2759",
    "1645187at2759",
    "176625at2759",
    "257318at2759",
    "296129at2759",
    "320059at2759",
    "383503at2759",
    "426305at2759",
    "543764at2759",
    "582756at2759",
    "625387at2759",
    "674160at2759",
    "692986at2759",
    "734666at2759",
    "736068at2759",
    "777920at2759",
    "834694at2759",
    "855834at2759",
    "869548at2759",
    "923657at2759",
    "939345at2759",
    "974865at2759",
    "976469at2759",
]

# Phage markers (genomad)
PHAGE_MODELS = [
    "GENOMAD.000023.VV",
    "GENOMAD.000026.VV",
    "GENOMAD.000172.VV",
    "GENOMAD.000217.VV",
    "GENOMAD.000254.VV",
    "GENOMAD.000616.VV",
    "GENOMAD.000979.VV",
    "GENOMAD.001596.VV",
    "GENOMAD.001683.VV",
    "GENOMAD.006053.VV",
    "GENOMAD.008072.VV",
    "GENOMAD.016854.VV",
    "GENOMAD.020215.VV",
    "GENOMAD.020681.VV",
    "GENOMAD.029205.VV",
    "GENOMAD.051630.VV",
    "GENOMAD.057391.VV",
    "GENOMAD.066074.VV",
    "GENOMAD.090963.VV",
    "GENOMAD.132751.VV",
]

# GVOG markers - 4 marker set
GVOG4M_MODELS = [
    "GVOGm0461",  # DNA topoisomerase II
    "GVOGm0022",  # DNA-directed RNA polymerase subunit beta
    "GVOGm0023",  # DNA-directed RNA polymerase subunit alpha
    "GVOGm0054",  # DNA polymerase elongation subunit family B
]

# GVOG markers - 8 marker set
GVOG8M_MODELS = [
    "GVOGm0013",  # DEAD/SNF2-like helicase
    "GVOGm0022",  # DNA-directed RNA polymerase subunit beta
    "GVOGm0023",  # DNA-directed RNA polymerase subunit alpha
    "GVOGm0054",  # DNA polymerase elongation subunit family B
    "GVOGm0172",  # putative transcription initiation factor IIB
    "GVOGm0461",  # DNA topoisomerase II
    "GVOGm0760",  # packaging ATPase
    "GVOGm0890",  # Poxvirus Late Transcription Factor VLTF3 like
]

# UNI56 markers
UNI56_MODELS = [
    "COG0013",  # alanyl-tRNA_synthetase
    "COG0016",  # phenylalanyl-tRNA_synthetase_alpha_subunit
    "COG0018",  # arginyl-tRNA_synthetase
    "COG0048",  # ribosomal_protein_S12
    "COG0049",  # ribosomal_protein_S7
    "COG0051",  # ribosomal_protein_S10
    "COG0052",  # ribosomal_protein_S2
    "COG0060",  # isoleucyl-tRNA_synthetase
    "COG0072",  # phenylalanyl-tRNA_synthetase_beta_subunit
    "COG0080",  # ribosomal_protein_L11
    "COG0081",  # ribosomal_protein_L1
    "COG0085",  # DNA-directed_RNA_polymerase_beta_subunit (RpoB)
    "COG0086",  # DNA-directed_RNA_polymerase_beta_prime_subunit (RpoC)
    "COG0087",  # ribosomal_protein_L3
    "COG0088",  # ribosomal_protein_L4
    "COG0089",  # ribosomal_protein_L23
    "COG0090",  # ribosomal_protein_L2
    "COG0091",  # ribosomal_protein_L22
    "COG0092",  # ribosomal_protein_S3
    "COG0093",  # ribosomal_protein_L14
    "COG0094",  # ribosomal_protein_L5
    "COG0096",  # ribosomal_protein_S8
    "COG0097",  # ribosomal_protein_L6P
    "COG0098",  # ribosomal_protein_S5
    "COG0099",  # ribosomal_protein_S13
    "COG0100",  # ribosomal_protein_S11
    "COG0102",  # ribosomal_protein_L13
    "COG0103",  # ribosomal_protein_S9
    "COG0127",  # Xanthosine_triphosphate_pyrophosphatase
    "COG0130",  # Pseudouridine_synthase
    "COG0164",  # ribonuclease_HII
    "COG0172",  # seryl-tRNA_synthetase
    "COG0184",  # ribosomal_protein_S15P
    "COG0185",  # ribosomal_protein_S19
    "COG0186",  # ribosomal_protein_S17
    "COG0193",  # peptidyl-tRNA_hydrolase
    "COG0197",  # ribosomal_protein_L16
    "COG0198",  # ribosomal_protein_L24
    "COG0200",  # ribosomal_protein_L15
    "COG0201",  # preprotein_translocase_subunit_SecY
    "COG0202",  # DNA-directed_RNA_polymerase_alpha_subunit (RpoA)
    "COG0216",  # protein_chain_release_factor_A
    "COG0233",  # ribosome_recycling_factor
    "COG0244",  # ribosomal_protein_L10
    "COG0255",  # ribosomal_protein_L29
    "COG0256",  # ribosomal_protein_L18
    "COG0343",  # queuine/archaeosine_tRNA-ribosyltransferase
    "COG0481",  # membrane_GTPase_LepA
    "COG0495",  # leucyl-tRNA_synthetase
    "COG0504",  # CTP_synthase
    "COG0519",  # GMP_synthase_PP-ATPase_domain
    "COG0532",  # translation_initiation_factor_2
    "COG0533",  # metal-dependent_proteases_with_possible_chaperone_activity
    "COG0541",  # signal_recognition_particle_GTPase
    "COG0691",  # tmRNA-binding_protein
    "COG0858",  # ribosome-binding_factor_A
]

# =============================================================================
# MCP markers (Major Capsid Protein) - NCLDV specific
# =============================================================================
NCLDV_MCP_MODELS = [
    "gamadvirusMCP",  # NCLDV major capsid protein (Gamadvirus)
    "GVOGm0003",  # NCLDV major capsid protein (GVOG)
    "OG1352",  # NCLDV major capsid protein
    "OG484",  # NCLDV major capsid protein
]

# All MCP markers (NCLDV + Mirusviricota)
MCP_MODELS = NCLDV_MCP_MODELS + MIRUS_CATEGORY_MODELS["MCP"]

# MRYA markers (Metagenomic Russian Yokohama-Asfarviridae)
MRYA_MODELS = [
    "HUH",  # HUH
    "HUHlong",  # HUHlong
    "VLTF2",  # Viral Late Transcription Factor VLTF2 like
    "VLTF3",  # Viral Late Transcription Factor VLTF3 like
    "ATPase",  # ATPase
    "gamadvirusMCP",  # NCLDV major capsid protein (Gamadvirus)
]

# capscan caps marker (plv_mcp_caps_*) -> Bellas&Sommaruga 2026 group, for the
# best-hit capscan_group field (avoids shared-MCP tree-vote fragmentation).
CAPS_MARKER_GROUP = {
    "plv_mcp_caps_Dino_Aquinto": "Dino_cluster",
    "plv_mcp_caps_GKS1_Gossevirus_Mavpol": "Gossevirus",
    "plv_mcp_caps_GKS2_Alpenseevirus_NCV_like": "Alpenseevirus",
    "plv_mcp_caps_GKS3_CHI_Aquinto": "GKS3",
    "plv_mcp_caps_Mavpol2": "Maverick-Polinton",
    "plv_mcp_caps_MMM_PC12_NA": "Metamonada_PLV",
    "plv_mcp_caps_MMM_PC40_NA": "Metamonada_PLV",
    "plv_mcp_caps_MMM_PC40b_NA": "Metamonada_PLV",
    "plv_mcp_caps_MMM_PC4_NA": "unclassified",
    "plv_mcp_caps_Trimcap_A4_Aquinto": "Trimcap_cluster_1",
    "plv_mcp_caps_Trimcap_B1_Aquinto": "Trimcap_cluster_2",
    "plv_mcp_caps_Trimcap_B2_Aquinto": "Trimcap_cluster_2",
    "plv_mcp_caps_Trimcap_B3_Aquinto": "Trimcap_cluster_2",
    "plv_mcp_caps_Zephyrvirus4_NCV_like": "Zephyrvirus",
    "plv_mcp_caps_Zephyrvirus5_NCV_like": "Zephyrvirus",
    "plv_mcp_caps_Nova6_VC35_Mavpol2": "Maverick-Polinton_2",
    "plv_mcp_caps_Nova7_Aquinto": "unclassified",
    "plv_mcp_caps_Nova8_NCV_like": "unclassified",
    "plv_mcp_caps_Nova9_Aquinto_Trimcap_related": "Trimcap_cluster_1",
    "plv_mcp_caps_PaM1_Aquinto": "Pam1",
    "plv_mcp_caps_PaM2_NCV_like": "Pam2",
    "plv_mcp_caps_VC124_VC134_NCV_like": "VC124-134",
    "plv_mcp_caps_PC6_20_Aquinto": "unclassified",
    "plv_mcp_caps_PgVV_Aquinto": "PgVV",
    "plv_mcp_caps_SP67_Aquinto": "SP67",
    "plv_mcp_caps_SP70_Aquinto": "SP70",
    "plv_mcp_caps_SP71_Aquinto": "SP71",
    "plv_mcp_caps_SP_Aquinto": "SP",
    "plv_mcp_caps_Trimcap_A1_Aquinto": "Trimcap_cluster_1",
    "plv_mcp_caps_Trimcap_A2_Aquinto": "Trimcap_cluster_1",
    "plv_mcp_caps_Trimcap_A3_74_Aquinto": "Trimcap_cluster_1",
    "plv_mcp_caps_Trimcap_A3_Aquinto": "Trimcap_cluster_1",
    "plv_mcp_caps_VC40_NCV_like": "VC40",
    "plv_mcp_caps_VC224_NA": "Mavpol2",
    "plv_mcp_caps_VC100_Mavpol1_like": "VC100-Mavpol1-like",
    "plv_mcp_caps_VC147_GKS3_Aquinto": "GKS3-CHI",
    "plv_mcp_caps_Cluster_187_Aquinto": "unclassified",
    "plv_mcp_caps_VC104_Aquinto": "Animal_PLV_group",
    "plv_mcp_caps_Cluster_276_NA": "unclassified",
    "plv_mcp_caps_Cluster_217_NA": "unclassified",
    "gamadvirusmcp_caps_VC289_Mriyavirus_NCV_like": "Mriyavirus",
    "plv_mcp_caps_VC218_Pam2_NCV_like": "unclassified",
    "plv_mcp_caps_Cluster_261_NA": "unclassified",
    "plv_mcp_caps_Cluster_152_Aquinto": "Trimcap_cluster_2",
    "plv_mcp_caps_Cluster_280_NA": "unclassified",
    "plv_mcp_caps_Zoravirus_VC102_NCV_like": "Zoravirus_VC102",
    "gamadvirusmcp_caps_Mriya_Nova10_NCV_like": "Mriyavirus",
    "gamadvirusmcp_caps_gamadvirus_MCP_NCV_like": "Mriyavirus",
    "gamadvirusmcp_caps_yaravirus_MCP_NCV_like": "Mriyavirus",
}

# =============================================================================
# Panel registry: every model name that a panel metric reads.
# =============================================================================
PANEL_MODELS: dict[str, list[str]] = {
    "gvog4": GVOG4M_MODELS,
    "gvog8": GVOG8M_MODELS,
    "busco": BUSCO_MODELS,
    "cog": UNI56_MODELS,
    "mrya": MRYA_MODELS,
    "phage": PHAGE_MODELS,
    "ncldv_mcp": NCLDV_MCP_MODELS,
    "mirus": MIRUS_MODELS,
}

# Panels that match every model name with a given prefix.
PANEL_PREFIXES: dict[str, str] = {
    **{
        f"vp_{category.lower()}": prefix
        for category, prefix in VP_CATEGORY_PREFIXES.items()
    },
    "plv": PLV_PREFIX,
}


def missing_panel_models(available_models: set[str]) -> dict[str, list[str]]:
    """Return, per panel, the panel models that ``available_models`` lacks.

    A prefix panel is missing, reported as ``"<prefix>*"``, when no available
    model starts with its prefix. The result is empty when every panel is present.
    """
    missing: dict[str, list[str]] = {}
    for panel, models in PANEL_MODELS.items():
        absent = [model for model in models if model not in available_models]
        if absent:
            missing[panel] = absent
    for panel, prefix in PANEL_PREFIXES.items():
        if not any(model.startswith(prefix) for model in available_models):
            missing[panel] = [f"{prefix}*"]
    return missing


# =============================================================================
# Marker groups (data-driven gene families sharing reference proteins).
# scripts/build_marker_groups.py builds the lookup src/config/marker_groups.tsv.
# Grouped models are combined into one {group}.faa -> one alignment -> one tree
# -> one nearest-neighbour vote, removing the per-model vote over-counting that
# otherwise fragments classification (e.g. an MCP cross-hitting 8 profiles -> 8 votes).
# =============================================================================
_MARKER_GROUPS_TSV = Path(__file__).with_name("marker_groups.tsv")
# Flags whose groups are safe to combine into a shared tree. REVIEW_mixed groups
# (possible over-merge) are excluded so their members stay one-model-per-tree.
_COMBINE_FLAGS = frozenset({"ok", "no_annotation"})
_GROUP_TABLE_FIELDS = 8


def load_marker_groups(
    path: Path = _MARKER_GROUPS_TSV,
    combine_flags: frozenset[str] = _COMBINE_FLAGS,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Return (group_to_models, model_to_group) for groups that are safe to combine.

    Reads the lookup written by ``scripts/build_marker_groups.py``. Only groups whose
    flag is in ``combine_flags`` are returned; REVIEW_mixed groups are skipped.

    Raises:
        FileNotFoundError: The table does not exist.
        ValueError: A row has too few fields, or the table holds no combinable group.
    """
    group_to_models: dict[str, list[str]] = {}
    model_to_group: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if line.startswith(("#", "group_name\t")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < _GROUP_TABLE_FIELDS:
                raise ValueError(
                    f"{path}: line {line_no} has {len(parts)} fields, "
                    f"expected {_GROUP_TABLE_FIELDS}"
                )
            name, flag = parts[0], parts[6]
            members = [m for m in parts[7].split(";") if m]
            if flag not in combine_flags or len(members) < 2:
                continue
            group_to_models[name] = members
            for model in members:
                model_to_group[model] = name
    if not group_to_models:
        raise ValueError(f"no marker groups found in {path}")
    return group_to_models, model_to_group


GROUP_TO_MODELS, MODEL_TO_GROUP = load_marker_groups()
