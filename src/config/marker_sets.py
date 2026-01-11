"""
Marker sets configuration for GVClass pipeline.
Based on docs/markers_gvclassV1.tsv

Marker Categories:
------------------
NCLDV (Nucleocytoviricota):
  - MCP: gamadvirusMCP, yaravirusMCP, PoxMCP, GVOGm0003, OG1352, OG484

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

# Legacy MIRUS models (v1.1.x database naming)
MIRUS_MODELS_LEGACY = [
    "mOG0000014",  # MCP (Mirusviricota)
    "mOG0000019",  # Portal (Mirusviricota)
    "mOG0000020",  # Triplex1 (Mirusviricota)
    "mOG0000030",  # Triplex2 (Mirusviricota)
    "mOG0000040",  # Maturation (Mirusviricota)
]

# All Mirus models (flat list for backward compatibility)
MIRUS_MODELS = [
    model for models in MIRUS_CATEGORY_MODELS.values() for model in models
] + MIRUS_MODELS_LEGACY

# BUSCO markers
BUSCO_MODELS = [
    "1001705at2759",
    "1003258at2759",
    "100698at2759",
    "1010730at2759",
    "1014314at2759",
    "1018517at2759",
    "1019762at2759",
    "1025450at2759",
    "1030907at2759",
    "1032689at2759",
    "1038775at2759",
    "1041560at2759",
    "1049599at2759",
    "1051021at2759",
    "1053181at2759",
    "1057950at2759",
    "1065019at2759",
    "1076134at2759",
    "1079130at2759",
    "1079827at2759",
    "1085752at2759",
    "1087488at2759",
    "1090038at2759",
    "1094121at2759",
    "1096688at2759",
    "1106766at2759",
    "1107630at2759",
    "1108845at2759",
    "1111142at2759",
    "1112002at2759",
    "1115196at2759",
    "1128607at2759",
    "1129824at2759",
    "1138059at2759",
    "1157302at2759",
    "1161199at2759",
    "1173229at2759",
    "1178688at2759",
    "1182451at2759",
    "1193442at2759",
    "1194691at2759",
    "1194797at2759",
    "1197019at2759",
    "1200489at2759",
    "1217666at2759",
    "1220881at2759",
    "1222562at2759",
    "1223488at2759",
    "1228942at2759",
    "1233814at2759",
    "1236198at2759",
    "1247641at2759",
    "1248958at2759",
    "1249159at2759",
    "1251252at2759",
    "1257440at2759",
    "1258856at2759",
    "1259741at2759",
    "1260807at2759",
    "1264469at2759",
    "1266231at2759",
    "1269244at2759",
    "1269649at2759",
    "1275837at2759",
    "1284731at2759",
    "1287094at2759",
    "1287401at2759",
    "1291729at2759",
    "1304061at2759",
    "1309031at2759",
    "1312453at2759",
    "1314980at2759",
    "1322299at2759",
    "1322642at2759",
    "1323575at2759",
    "1324510at2759",
    "1338131at2759",
    "1339553at2759",
    "1342242at2759",
    "1346165at2759",
    "1346432at2759",
    "1348942at2759",
    "1355894at2759",
    "1358374at2759",
    "1364586at2759",
    "1370285at2759",
    "1370304at2759",
    "1377237at2759",
    "1379600at2759",
    "1380409at2759",
    "1382842at2759",
    "1395649at2759",
    "1398309at2759",
    "1404162at2759",
    "1405073at2759",
    "1405146at2759",
    "1407446at2759",
    "1421503at2759",
    "1423847at2759",
    "142542at2759",
    "1426075at2759",
    "1428265at2759",
    "1428854at2759",
    "1430056at2759",
    "1432328at2759",
    "1434266at2759",
    "1440961at2759",
    "1441854at2759",
    "1442062at2759",
    "1445102at2759",
    "1450538at2759",
    "1454155at2759",
    "1455730at2759",
    "1459797at2759",
    "1474000at2759",
    "1479417at2759",
    "1488235at2759",
    "1488436at2759",
    "1490022at2759",
    "1504863at2759",
    "1513531at2759",
    "1525971at2759",
    "1526152at2759",
    "1530008at2759",
    "1538526at2759",
    "1545004at2759",
    "1552706at2759",
    "1558822at2759",
    "156083at2759",
    "1563319at2759",
    "1567796at2759",
    "1575179at2759",
    "1576404at2759",
    "1588798at2759",
    "1593937at2759",
    "160593at2759",
    "1617752at2759",
    "1620056at2759",
    "1623701at2759",
    "1626636at2759",
    "1633672at2759",
    "1645187at2759",
    "166920at2759",
    "176625at2759",
    "179362at2759",
    "245208at2759",
    "257318at2759",
    "261328at2759",
    "261419at2759",
    "270107at2759",
    "271586at2759",
    "290630at2759",
    "293315at2759",
    "296129at2759",
    "299766at2759",
    "306227at2759",
    "320059at2759",
    "324863at2759",
    "325552at2759",
    "330169at2759",
    "331411at2759",
    "341721at2759",
    "345441at2759",
    "355408at2759",
    "369837at2759",
    "375960at2759",
    "383503at2759",
    "388820at2759",
    "390348at2759",
    "392369at2759",
    "39650at2759",
    "396755at2759",
    "418107at2759",
    "426305at2759",
    "450058at2759",
    "453044at2759",
    "457861at2759",
    "464990at2759",
    "491869at2759",
    "513979at2759",
    "543764at2759",
    "549762at2759",
    "551907at2759",
    "570797at2759",
    "582756at2759",
    "598949at2759",
    "604979at2759",
    "621827at2759",
    "625387at2759",
    "655400at2759",
    "664730at2759",
    "671536at2759",
    "671846at2759",
    "673132at2759",
    "674160at2759",
    "674169at2759",
    "679187at2759",
    "679771at2759",
    "687505at2759",
    "692986at2759",
    "708105at2759",
    "719531at2759",
    "720340at2759",
    "721605at2759",
    "722805at2759",
    "734341at2759",
    "734666at2759",
    "736068at2759",
    "751335at2759",
    "759498at2759",
    "761109at2759",
    "768809at2759",
    "774318at2759",
    "777920at2759",
    "779909at2759",
    "801857at2759",
    "814241at2759",
    "817008at2759",
    "834694at2759",
    "836599at2759",
    "83779at2759",
    "844748at2759",
    "855834at2759",
    "858842at2759",
    "865202at2759",
    "866359at2759",
    "869548at2759",
    "878143at2759",
    "87842at2759",
    "887370at2759",
    "891751at2759",
    "898782at2759",
    "901894at2759",
    "905026at2759",
    "911863at2759",
    "917326at2759",
    "918816at2759",
    "919955at2759",
    "923657at2759",
    "924753at2759",
    "931188at2759",
    "937275at2759",
    "937686at2759",
    "939345at2759",
    "944899at2759",
    "946128at2759",
    "956854at2759",
    "97116at2759",
    "973442at2759",
    "974865at2759",
    "975158at2759",
    "975557at2759",
    "976469at2759",
    "981902at2759",
    "996662at2759",
]

# Phage markers (genomad)
PHAGE_MODELS = [
    "genomad000023",
    "genomad000026",
    "genomad000172",
    "genomad000217",
    "genomad000254",
    "genomad000616",
    "genomad000979",
    "genomad001596",
    "genomad001683",
    "genomad006053",
    "genomad008072",
    "genomad016854",
    "genomad020215",
    "genomad020681",
    "genomad029205",
    "genomad051630",
    "genomad057391",
    "genomad066074",
    "genomad090963",
    "genomad132751",
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
    "yaravirusMCP",  # NCLDV major capsid protein (Yaravirus)
    "PoxMCP",  # NCLDV major capsid protein (Poxviruses)
    "GVOGm0003",  # NCLDV major capsid protein (GVOG)
    "OG1352",  # NCLDV major capsid protein
    "OG484",  # NCLDV major capsid protein
]

# All MCP markers (NCLDV + Mirusviricota)
MCP_MODELS = NCLDV_MCP_MODELS + [
    "mOG0000014",  # MCP (Mirusviricota) - legacy naming
    "Mirus_MCP",  # MCP (Mirusviricota) - v1.2 naming
    "Mirus_JellyRoll",  # MCP (Mirusviricota) - v1.2 naming
]

# MRYA markers (Metagenomic Russian Yokohama-Asfarviridae)
MRYA_MODELS = [
    "HUH",  # HUH
    "HUHlong",  # HUHlong
    "VLTF2",  # Viral Late Transcription Factor VLTF2 like
    "VLTF3",  # Viral Late Transcription Factor VLTF3 like
    "ATPase",  # ATPase
    "gamadvirusMCP",  # NCLDV major capsid protein (Gamadvirus)
]
