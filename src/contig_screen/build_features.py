"""Turn the extraction tables of a bin into one feature row per window.

Feature groups (column prefix): ``s_`` gene structure and RBS, ``c_`` composition,
``h_`` reference hits. A bin is scored against every reference subject (reference
hold-out ``none``), so one row per window is built.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HOLDOUT = "none"
GROUPS = ("BAC", "ARC", "EUK", "EVE", "ORG", "NCLDV", "MIRUS", "PHAGE")
GROUP_OF = {
    "BAC": "BAC",
    "ARC": "ARC",
    "EUK": "EUK",
    "EUK-pEVE": "EVE",
    "MITO": "ORG",
    "PLASTID": "ORG",
    "NCLDV": "NCLDV",
    "MIRUS": "MIRUS",
    "PHAGE": "PHAGE",
    "PPV": "PHAGE",
}
CELLULAR = ("BAC", "ARC", "EUK", "ORG")
ID_BINS = (0.0, 35.0, 50.0, 70.0, 90.0, 100.1)
ID_NAMES = ("lt35", "35_50", "50_70", "70_90", "ge90")
SPA_BINS = (0.0, 1.0, 1.8, 2.6, 3.6, 100.0)
SPA_NAMES = ("lt1", "1_1.8", "1.8_2.6", "2.6_3.6", "ge3.6")


@dataclass(frozen=True, slots=True)
class BinSpec:
    """A value column of the best hit, its bin edges and the output column prefix."""

    column: str
    edges: tuple[float, ...]
    names: tuple[str, ...]
    prefix: str


IDENTITY_SPEC = BinSpec("best_identity", ID_BINS, ID_NAMES, "h_icell")
BITSCORE_SPEC = BinSpec("best_spa", SPA_BINS, SPA_NAMES, "h_bcell")
RANK_NAMES = ("tax_domain", "phylum", "class", "order", "family", "genus", "species")
MODEL_CLASSES = {
    "og": re.compile(r"^(OG|GVOG|VLTF|Mirus|MRYA|NCLDV)", re.IGNORECASE),
    "cog": re.compile(r"^COG"),
    "phage": re.compile(r"^(GENOMAD|PLV|plv|VP_|gamad)"),
}
CODONS = tuple(a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_ORDER = tuple(a + b + c for a in BASES for b in BASES for c in BASES)
AA_OF = dict(zip(CODON_ORDER, AMINO, strict=True))
AA_LETTERS = tuple(sorted({aa for aa in AMINO if aa != "*"}))
RBS_MOTIFS = (
    "GGA/GAG/AGG",
    "GGAG/GAGG",
    "AGGAG(G)/GGAGG",
    "GGAGG",
    "AGGAG",
    "AGGA",
    "GGAG",
    "AGxAG",
    "3Base/5BMM",
    "4Base/6BMM",
    "GGxGG",
    "None",
)
START_TYPES = ("ATG", "GTG", "TTG", "Edge")


def load_labels(labels_file: Path, label_domain: Path) -> pd.DataFrame:
    """Return label id to domain group and the taxon name at each rank."""
    domain = pd.read_parquet(label_domain).drop_duplicates("label_id")
    table = pd.read_parquet(labels_file)
    ranks = table["taxonomy"].str.split("|", expand=True).reindex(columns=range(7))
    ranks.columns = list(RANK_NAMES)
    ranks = ranks.fillna("")
    ranks["label_id"] = table["label_id"].to_numpy()
    merged = domain.merge(ranks, on="label_id", how="left").fillna("")
    merged["group"] = merged["domain"].map(GROUP_OF).fillna("OTHER")
    return merged.set_index("label_id")


def unit_features(
    extract_dir: Path, genome_id: str, genome_class: str, labels: pd.DataFrame
) -> pd.DataFrame:
    """Return the feature rows of one bin at reference hold-out ``none``."""
    stem = extract_dir / genome_id
    contigs = pd.read_parquet(f"{stem}.contigs.parquet")
    genes = pd.read_parquet(f"{stem}.genes.parquet")
    hits = pd.read_parquet(f"{stem}.hits.parquet")
    base = structure_and_composition(contigs, genes)
    base["class"] = genome_class
    hits = annotate_hits(hits, genes, labels)
    frame = base.join(hit_features(contigs, genes, hits))
    frame["holdout"] = HOLDOUT
    return frame.reset_index()


def annotate_hits(
    hits: pd.DataFrame, genes: pd.DataFrame, labels: pd.DataFrame
) -> pd.DataFrame:
    """Add label, group, lineage names and coverage to the hit table."""
    out = hits.copy()
    out["label"] = out["subject"].str.split("|", n=1).str[0]
    aa_len = genes.drop_duplicates("prot_hash").set_index("prot_hash")["aa_len"]
    out["coverage"] = out["aln_len"] / out["prot_hash"].map(aa_len)
    out["spa"] = out["score"] / out["prot_hash"].map(aa_len).clip(lower=1)
    info = labels.reindex(out["label"])
    for column in ("group", "phylum", "class", "order"):
        out[column] = info[column].to_numpy()
    out["group"] = out["group"].fillna("OTHER")
    return out


def structure_and_composition(
    contigs: pd.DataFrame, genes: pd.DataFrame
) -> pd.DataFrame:
    """Return gene structure (``s_``) and composition (``c_``) features per contig."""
    frame = contigs.set_index("contig_id")
    out = pd.DataFrame(index=frame.index)
    out["genome_id"] = frame["genome_id"]
    out["class"] = frame["class"]
    out["length_kb"] = frame["length_kb"]
    out["length_bp"] = frame["length_bp"]
    out["source_contig"] = frame["source_contig"]
    out["source_start"] = frame["source_start"]
    out["s_n_genes"] = frame["n_genes"]
    out["s_genes_per_kb"] = frame["n_genes"] / frame["length_bp"] * 1000
    out["s_coding_density"] = frame["coding_bp"] / frame["length_bp"]
    out["c_gc"] = frame["gc"]
    codon = frame[[f"cod_{c}" for c in CODONS]].to_numpy(dtype=float)
    out = out.join(composition(codon, frame.index))
    return out.join(gene_structure(genes))


def composition(codon: np.ndarray, index: pd.Index) -> pd.DataFrame:
    """Return codon, amino acid and codon-position GC fractions from codon counts."""
    total = codon.sum(axis=1)
    safe = np.where(total > 0, total, np.nan)
    data: dict[str, np.ndarray] = {}
    for position in range(3):
        gc = [i for i, c in enumerate(CODONS) if c[position] in "GC"]
        data[f"c_gc{position + 1}"] = codon[:, gc].sum(axis=1) / safe
    for aa in AA_LETTERS:
        cols = [i for i, c in enumerate(CODONS) if AA_OF[c] == aa]
        data[f"c_aa_{aa}"] = codon[:, cols].sum(axis=1) / safe
    data["c_stop_frac"] = (
        codon[:, [i for i, c in enumerate(CODONS) if AA_OF[c] == "*"]].sum(axis=1)
        / safe
    )
    for i, codon_name in enumerate(CODONS):
        if AA_OF[codon_name] != "*":
            data[f"c_cod_{codon_name}"] = codon[:, i] / safe
    return pd.DataFrame(data, index=index)


def gene_structure(genes: pd.DataFrame) -> pd.DataFrame:
    """Return per-contig statistics of strand, spacing, starts, RBS and scores."""
    ordered = genes.sort_values(["contig_id", "begin"])
    grouped = ordered.groupby("contig_id", sort=False)
    prev_end = grouped["end"].shift(1)
    prev_strand = grouped["strand"].shift(1)
    ordered = ordered.assign(
        gap=ordered["begin"] - prev_end - 1,
        switch=(ordered["strand"] != prev_strand) & prev_strand.notna(),
        length_nt=ordered["end"] - ordered["begin"] + 1,
        has_prev=prev_end.notna(),
    )
    ordered["overlap"] = (ordered["gap"] < 0) & ordered["has_prev"]
    ordered["gap_pos"] = ordered["gap"].where(
        ordered["has_prev"] & (ordered["gap"] >= 0)
    )
    ordered["score_per_aa"] = ordered["score"] / ordered["aa_len"].clip(lower=1)
    grouped = ordered.groupby("contig_id", sort=False)
    out = pd.DataFrame(
        {
            "s_gene_len_mean": grouped["aa_len"].mean(),
            "s_gene_len_sd": grouped["aa_len"].std(),
            "s_gene_len_median": grouped["aa_len"].median(),
            "s_strand_switch_rate": grouped["switch"].mean(),
            "s_overlap_frac": grouped["overlap"].mean(),
            "s_gap_mean": grouped["gap_pos"].mean(),
            "s_gap_median": grouped["gap_pos"].median(),
            "s_partial_frac": grouped["partial"].apply(lambda s: (s > 0).mean()),
            "s_minus_frac": grouped["strand"].apply(lambda s: (s < 0).mean()),
            "s_gene_gc_mean": grouped["gc_cont"].mean(),
            "s_gene_gc_sd": grouped["gc_cont"].std(),
            "s_conf_mean": grouped["conf"].mean(),
            "s_score_per_aa_mean": grouped["score_per_aa"].mean(),
            "s_rscore_mean": grouped["rscore"].mean(),
            "s_sscore_mean": grouped["sscore"].mean(),
            "s_tscore_mean": grouped["tscore"].mean(),
            "s_uscore_mean": grouped["uscore"].mean(),
            "s_cscore_mean": grouped["cscore"].mean(),
        }
    )
    for start in START_TYPES:
        out[f"s_start_{start}"] = grouped["start_type"].apply(
            lambda s, t=start: (s == t).mean()
        )
    for motif in RBS_MOTIFS:
        out[f"s_rbs_{motif}"] = grouped["rbs_motif"].apply(
            lambda s, m=motif: (s == m).mean()
        )
    out["s_rbs_spacer_none"] = grouped["rbs_spacer"].apply(
        lambda s: (s == "None").mean()
    )
    return out


def hit_features(
    contigs: pd.DataFrame, genes: pd.DataFrame, hits: pd.DataFrame
) -> pd.DataFrame:
    """Return hit-based (``h_``) features per contig for the kept hits."""
    index = pd.Index(contigs["contig_id"], name="contig_id")
    gene = genes[["contig_id", "prot_hash", "best_model", "aa_len", "gc_cont"]].copy()
    best = best_hits(hits)
    gene = gene.join(best, on="prot_hash")
    gene["has_hmm"] = gene["best_model"].notna()
    gene["has_ref"] = gene["best_group"].notna()
    out = pd.DataFrame(index=index)
    grouped = gene.groupby("contig_id")
    n_genes = grouped.size().reindex(index).clip(lower=1)
    out["h_frac_hmm"] = grouped["has_hmm"].sum().reindex(index) / n_genes
    out["h_frac_ref"] = grouped["has_ref"].sum().reindex(index) / n_genes
    hmm_only = gene["has_hmm"] & ~gene["has_ref"]
    out["h_frac_hmm_only"] = (
        hmm_only.groupby(gene["contig_id"]).sum().reindex(index) / n_genes
    )
    out["h_frac_nohit"] = 1.0 - out["h_frac_hmm"]
    nohit = gene[~gene["has_hmm"]]
    out["h_nohit_aa_len"] = nohit.groupby("contig_id")["aa_len"].mean().reindex(index)
    out["h_nohit_gc"] = nohit.groupby("contig_id")["gc_cont"].mean().reindex(index)
    out = out.join(model_class_fractions(gene, index, n_genes))
    out = out.join(group_shares(gene, index, n_genes))
    out = out.join(group_bin_matrix(gene, index, n_genes, IDENTITY_SPEC))
    out = out.join(group_bin_matrix(gene, index, n_genes, BITSCORE_SPEC))
    out = out.join(
        best_value_by_group(gene, hits, index, column="identity", prefix="h_idmax")
    )
    out = out.join(
        best_value_by_group(gene, hits, index, column="spa", prefix="h_bsmax")
    )
    out = out.join(coherence(gene, index))
    out["h_identity_mean"] = (
        gene.groupby("contig_id")["best_identity"].mean().reindex(index)
    )
    out["h_score_per_aa_mean"] = (
        (gene["best_score"] / gene["aa_len"].clip(lower=1))
        .groupby(gene["contig_id"])
        .mean()
        .reindex(index)
    )
    return out


def best_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """Return the best kept hit of every protein."""
    columns = [
        "best_group",
        "best_identity",
        "best_score",
        "best_spa",
        "best_phylum",
        "best_class",
        "best_order",
    ]
    if hits.empty:
        return pd.DataFrame(columns=columns)
    top = hits.sort_values(["prot_hash", "score"], ascending=[True, False])
    top = top.drop_duplicates("prot_hash").set_index("prot_hash")
    return top.rename(
        columns={
            "group": "best_group",
            "identity": "best_identity",
            "score": "best_score",
            "spa": "best_spa",
            "phylum": "best_phylum",
            "class": "best_class",
            "order": "best_order",
        }
    )[columns]


def model_class_fractions(
    gene: pd.DataFrame, index: pd.Index, n_genes: pd.Series
) -> pd.DataFrame:
    """Return the fraction of genes whose best catalog model is of each class."""
    out = pd.DataFrame(index=index)
    names = gene["best_model"].fillna("")
    for key, pattern in MODEL_CLASSES.items():
        flag = names.map(lambda n, p=pattern: bool(p.match(n)))
        counts = flag.groupby(gene["contig_id"]).sum().reindex(index)
        out[f"h_model_{key}"] = counts / n_genes
    return out


def group_bin_matrix(
    gene: pd.DataFrame, index: pd.Index, n_genes: pd.Series, spec: BinSpec
) -> pd.DataFrame:
    """Return the share of genes per best-hit group and value bin."""
    have = gene[gene["has_ref"]].copy()
    have["bin"] = pd.cut(
        have[spec.column], bins=spec.edges, labels=spec.names, right=False
    ).astype(str)
    have["cell"] = have["best_group"].astype(str) + "_" + have["bin"]
    table = pd.crosstab(have["contig_id"], have["cell"]).reindex(index, fill_value=0)
    wanted = [f"{g}_{b}" for g, b in itertools.product(GROUPS, spec.names)]
    table = table.reindex(columns=wanted, fill_value=0).div(n_genes, axis=0)
    table.columns = [f"{spec.prefix}_{c}" for c in table.columns]
    return table


def group_shares(
    gene: pd.DataFrame, index: pd.Index, n_genes: pd.Series
) -> pd.DataFrame:
    """Return the share of genes whose best hit is in each group."""
    counts = pd.crosstab(
        gene.loc[gene["has_ref"], "contig_id"], gene.loc[gene["has_ref"], "best_group"]
    )
    counts = counts.reindex(index, fill_value=0).reindex(columns=GROUPS, fill_value=0)
    shares = counts.div(n_genes, axis=0)
    shares.columns = [f"h_share_{g}" for g in GROUPS]
    shares["h_share_cellular"] = shares[[f"h_share_{g}" for g in CELLULAR]].sum(axis=1)
    shares["h_share_giant"] = shares[["h_share_NCLDV", "h_share_MIRUS"]].sum(axis=1)
    return shares


def best_value_by_group(
    gene: pd.DataFrame,
    hits: pd.DataFrame,
    index: pd.Index,
    *,
    column: str,
    prefix: str,
) -> pd.DataFrame:
    """Return, per contig, the mean over genes of the best ``column`` to each group."""
    out = pd.DataFrame(
        index=index, columns=[f"{prefix}_{g}" for g in GROUPS], dtype=float
    )
    if hits.empty:
        return out
    per = hits.pivot_table(
        index="prot_hash", columns="group", values=column, aggfunc="max"
    )
    joined = gene[["contig_id", "prot_hash"]].join(per, on="prot_hash")
    for group in GROUPS:
        if group in per.columns:
            out[f"{prefix}_{group}"] = (
                joined.groupby("contig_id")[group].mean().reindex(index)
            )
    return out


def coherence(gene: pd.DataFrame, index: pd.Index) -> pd.DataFrame:
    """Return phylum and class level agreement of the cellular best hits."""
    out = pd.DataFrame(index=index)
    cell = gene[gene["best_group"].isin(CELLULAR)]
    for rank in ("phylum", "class"):
        column = f"best_{rank}"
        grouped = cell.groupby("contig_id")[column]
        top = grouped.agg(lambda s: s.value_counts().iloc[0] / len(s))
        out[f"h_cell_top_{rank}_share"] = top.reindex(index)
        out[f"h_cell_n_{rank}"] = grouped.nunique().reindex(index)
    shannon = (
        gene[gene["has_ref"]]
        .groupby("contig_id")["best_group"]
        .agg(lambda s: entropy(s.value_counts().to_numpy()))
    )
    out["h_group_entropy"] = shannon.reindex(index)
    return out


def entropy(counts: np.ndarray) -> float:
    """Return the Shannon entropy in bits of a count vector."""
    p = counts / counts.sum()
    return float(-(p * np.log2(p)).sum())
