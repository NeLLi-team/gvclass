"""
Marker sets configuration for GVClass pipeline.
Based on docs/markers_gvclassV1.tsv
"""

# MIRUS markers (Mirusviricota)
MIRUS_MODELS = [
    "mOG0000014",  # MCP (Mirusviricota)
    "mOG0000019",  # Portal (Mirusviricota)
    "mOG0000020",  # Triplex1 (Mirusviricota)
    "mOG0000030",  # Triplex2 (Mirusviricota)
    "mOG0000040",  # Maturation (Mirusviricota)
]

# BUSCO markers (only those present in the HMM database)
BUSCO_MODELS = [
    "1003258at2759", "1019762at2759", "1032689at2759", "1038775at2759", "1057950at2759",
    "1079130at2759", "1079827at2759", "1090038at2759", "1106766at2759", "1111142at2759",
    "1157302at2759", "1161199at2759", "1178688at2759", "1193442at2759", "1194691at2759",
    "1197019at2759", "1200489at2759", "1220881at2759", "1222562at2759", "1236198at2759",
    "1249159at2759", "1260807at2759", "1275837at2759", "1284731at2759", "1291729at2759",
    "1309031at2759", "1322299at2759", "1338131at2759", "1346165at2759", "1346432at2759",
    "1348942at2759", "1380409at2759", "1382842at2759", "1423847at2759", "1430056at2759",
    "1432328at2759", "1434266at2759", "1455730at2759", "1488436at2759", "1504863at2759",
    "1530008at2759", "1552706at2759", "1567796at2759", "1576404at2759", "1588798at2759",
    "1617752at2759", "1623701at2759", "1645187at2759", "176625at2759", "257318at2759",
    "296129at2759", "320059at2759", "383503at2759", "426305at2759", "543764at2759",
    "582756at2759", "625387at2759", "674160at2759", "692986at2759", "734666at2759",
    "736068at2759", "777920at2759", "834694at2759", "855834at2759", "869548at2759",
    "923657at2759", "939345at2759", "974865at2759", "976469at2759",
]

# Phage markers (genomad)
PHAGE_MODELS = [
    "genomad000023", "genomad000026", "genomad000172", "genomad000217", "genomad000254",
    "genomad000616", "genomad000979", "genomad001596", "genomad001683", "genomad006053",
    "genomad008072", "genomad016854", "genomad020215", "genomad020681", "genomad029205",
    "genomad051630", "genomad057391", "genomad066074", "genomad090963", "genomad132751",
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

# MCP markers (Major Capsid Protein)
MCP_MODELS = [
    "gamadvirusMCP",  # NCLDV major capsid protein (Gamadvirus)
    "yaravirusMCP",   # NCLDV major capsid protein (Yaravirus)
    "PoxMCP",         # NCLDV major capsid protein (Poxviruses)
    "GVOGm0003",      # NCLDV major capsid protein
    "mOG0000014",     # MCP (Mirusviricota)
]

# MRYA markers (Metagenomic Russian Yokohama-Asfarviridae)
MRYA_MODELS = [
    "HUH",            # HUH
    "HUHlong",        # HUHlong
    "VLTF2",          # Viral Late Transcription Factor VLTF2 like
    "VLTF3",          # Viral Late Transcription Factor VLTF3 like
    "ATPase",         # ATPase
    "gamadvirusMCP",  # NCLDV major capsid protein (Gamadvirus)
]