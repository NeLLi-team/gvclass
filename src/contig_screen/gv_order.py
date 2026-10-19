"""Gene order and copy number of catalog markers on the contigs of a giant-virus bin.

A clean bin holds one genome, and most order-level markers (``OG``, ``GVOGm``) occur
once in it. A second copy of a marker, and above all a second copy of a short
stretch of markers in the same order, points to a second genome in the bin. This
module turns the gene table of ``gv_genes.py`` into per-contig features of that
evidence. The synthetic bins and the fit of the priors are not ported; the priors
come from the model bundle.

Per contig the features are: copies of its markers on other contigs, weighted by how
rarely the marker is duplicated in single genomes (priors from reference genomes),
the amino-acid identity of those copies, shared neighbour pairs of markers (rough
order: pairs up to ``SKIP`` markers apart), and the composition distance of the
contig to the rest of the bin.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pyswrd

VIRAL_PREFIXES = ("OG", "GVOGm", "Mirus", "MRYA")
SKIP = 3
MIN_FRAGMENT_BP = 5000
FLOOR = 0.01
SINGLE_COPY_P = 0.1
TAXON_CLASSES = ("NCLDV", "MIRUS", "OTHERV")
TAXON_LEVELS = ("family", "genus")
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def canonical_table() -> tuple[np.ndarray, int]:
    """Return the 4-mer code (0 to 255) to canonical class table and the class count."""
    pair = {"A": "T", "C": "G", "G": "C", "T": "A"}
    index: dict[str, int] = {}
    table = np.zeros(256, dtype=np.int16)
    for code in range(256):
        kmer = "".join("ACGT"[(code >> (2 * (3 - i))) & 3] for i in range(4))
        reverse = "".join(pair[b] for b in reversed(kmer))
        table[code] = index.setdefault(min(kmer, reverse), len(index))
    return table, len(index)


CANONICAL, N_CLASSES = canonical_table()
LOOKUP = np.full(256, 4, dtype=np.uint8)
for _letter, _value in BASE.items():
    LOOKUP[ord(_letter)] = _value


def composition(sequence: str) -> np.ndarray:
    """Return canonical 4-mer frequencies followed by the GC fraction."""
    codes = LOOKUP[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    out = np.zeros(N_CLASSES + 1)
    if len(codes) < 4:
        return out
    valid = codes < 4
    c = codes.astype(np.int32)
    word = c[:-3] * 64 + c[1:-2] * 16 + c[2:-1] * 4 + c[3:]
    keep = valid[:-3] & valid[1:-2] & valid[2:-1] & valid[3:]
    counts = np.bincount(CANONICAL[word[keep]], minlength=N_CLASSES)
    if counts.sum() > 0:
        out[:N_CLASSES] = counts / counts.sum()
    gc = int(((codes == 1) | (codes == 2)).sum())
    out[N_CLASSES] = gc / max(int(valid.sum()), 1)
    return out


CODON_COLUMNS = [f"cod_{i}" for i in range(64)]
THIRD_GC = np.array([(i % 4) in (1, 2) for i in range(64)])


def gene_codon_counts(sequence: str, begin: int, end: int, strand: int) -> np.ndarray:
    """Return the 64 codon counts of one gene (reverse strand complemented)."""
    codes = LOOKUP[np.frombuffer(sequence[begin - 1 : end].encode("ascii"), np.uint8)]
    if strand < 0:
        codes = np.where(codes < 4, 3 - codes, 4)[::-1]
    n = len(codes) // 3 * 3
    triple = codes[:n].reshape(-1, 3).astype(np.int32)
    ok = (triple < 4).all(axis=1)
    index = triple[ok, 0] * 16 + triple[ok, 1] * 4 + triple[ok, 2]
    return np.bincount(index, minlength=64)


def codon_vector(genes: pd.DataFrame) -> np.ndarray:
    """Return codon frequencies (64) and the GC fraction of third positions."""
    out = np.zeros(65)
    if genes.empty or CODON_COLUMNS[0] not in genes:
        return out
    counts = genes[CODON_COLUMNS].to_numpy(dtype=float).sum(axis=0)
    total = counts.sum()
    if total > 0:
        out[:64] = counts / total
        out[64] = counts[THIRD_GC].sum() / total
    return out


def structure_vector(genes: pd.DataFrame, length: int) -> np.ndarray:
    """Return coding density, mean gene length, strand switch rate, genes per kb."""
    if genes.empty:
        return np.zeros(4)
    ordered = genes.sort_values("begin")
    strand = ordered["strand"].to_numpy()
    switch = float((strand[1:] != strand[:-1]).mean()) if len(strand) > 1 else 0.0
    coding = float((ordered["end"] - ordered["begin"] + 1).sum()) / length
    return np.array(
        [
            min(coding, 1.5),
            float(ordered["aa_len"].mean()) / 100.0,
            switch,
            len(ordered) / (length / 1000.0),
        ]
    )


@dataclass
class Taxonomy:
    """Best catalog hits of the genes of one contig.

    ``n_hit`` counts genes with a hit of any class, ``n_ncldv`` those with a hit to a
    NCLDV record; ``family`` and ``genus`` count the codes of those NCLDV hits and
    ``classes`` the class of the best hit of every gene with a hit.
    """

    n_genes: int
    n_hit: int
    n_ncldv: int
    family: dict[str, int]
    genus: dict[str, int]
    classes: dict[str, int]
    pident_sum: float
    bits_sum: float


@dataclass
class Contig:
    """One contig (or fragment) of a bin with its marker genes in order."""

    name: str
    genome: str
    length: int
    models: list[str] = field(default_factory=list)
    proteins: list[str] = field(default_factory=list)
    hashes: list[str] = field(default_factory=list)
    comp: np.ndarray | None = None
    codon: np.ndarray | None = None
    struct: np.ndarray | None = None
    tax: Taxonomy | None = None


@dataclass
class Priors:
    """Single-copy priors of markers from reference genomes."""

    p_multi: dict[str, float]
    default: float
    mu: np.ndarray
    sd: np.ndarray
    extra: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)

    def single_copy(self, model: str) -> bool:
        """Return True when ``model`` is rarely duplicated in single genomes."""
        return self.p_multi.get(model, self.default) <= SINGLE_COPY_P

    def weight(self, model: str) -> float:
        """Return how surprising a second copy of ``model`` is (log10 units)."""
        return -math.log10(max(self.p_multi.get(model, self.default), FLOOR))


def is_marker(model: str) -> bool:
    """Return True for the order-level and core giant-virus models."""
    return model.startswith(VIRAL_PREFIXES)


def contig_taxonomy(genes: pd.DataFrame) -> Taxonomy | None:
    """Return the hit counts of the genes of a contig, or None without hit columns."""
    if "tx_class" not in genes:
        return None
    classes = genes["tx_class"].fillna("")
    ncldv = genes[genes["tx_pident"].notna()]
    return Taxonomy(
        n_genes=len(genes),
        n_hit=int((classes != "").sum()),
        n_ncldv=len(ncldv),
        family=dict(Counter(c for c in ncldv["tx_family"].fillna("") if c)),
        genus=dict(Counter(c for c in ncldv["tx_genus"].fillna("") if c)),
        classes=dict(Counter(c for c in classes if c)),
        pident_sum=float(ncldv["tx_pident"].sum()),
        bits_sum=float(ncldv["tx_bits"].sum()),
    )


def contig_from_genes(
    genes: pd.DataFrame, name: str, genome: str, start: int, end: int
) -> Contig:
    """Return the marker genes of one source contig between ``start`` and ``end``."""
    within = genes[(genes["begin"] > start) & (genes["end"] <= end)]
    inside = within[within["model"].map(is_marker).astype(bool)]
    inside = inside.sort_values("begin")
    return Contig(
        name=name,
        genome=genome,
        length=end - start,
        models=inside["model"].tolist(),
        proteins=inside["protein"].tolist(),
        hashes=inside["prot_hash"].tolist(),
        codon=codon_vector(within),
        struct=structure_vector(within, end - start),
        tax=contig_taxonomy(within),
    )


def top_code(counts: dict[str, int]) -> str:
    """Return the commonest code (the smaller string on ties), or an empty string."""
    positive = {code: n for code, n in counts.items() if n > 0}
    if not positive:
        return ""
    return min(positive, key=lambda code: (-positive[code], code))


def rest_code(total: Counter[str], own: dict[str, int]) -> str:
    """Return the commonest code of the bin without the counts of one contig."""
    return top_code({code: n - own.get(code, 0) for code, n in total.items()})


def share_of(count: float, whole: float) -> float:
    """Return ``count / whole``, or NaN when ``whole`` is zero."""
    return count / whole if whole else float("nan")


def taxonomy_features(contigs: list[Contig]) -> pd.DataFrame:
    """Return per contig the agreement of its catalog hits with the rest of the bin.

    Family and genus codes of the NCLDV hits are compared with the commonest code of
    the other contigs; the class shares of the best hits (NCLDV, Mirus, other
    viruses) are compared with the share in the other contigs. A contig without
    hits, or a bin of one contig, has missing agreement.
    """
    taxa = [c.tax for c in contigs if c.tax is not None]
    totals = {
        level: sum((Counter(getattr(t, level)) for t in taxa), Counter())
        for level in TAXON_LEVELS
    }
    class_total = sum((Counter(t.classes) for t in taxa), Counter())
    hit_total = sum(t.n_hit for t in taxa)
    rows = []
    for tax in taxa:
        row = {
            "tx_hit_frac": share_of(tax.n_hit, tax.n_genes),
            "tx_ncldv_frac": share_of(tax.n_ncldv, tax.n_genes),
            "tx_pident": share_of(tax.pident_sum, tax.n_ncldv),
            "tx_bits": share_of(tax.bits_sum, tax.n_ncldv),
        }
        for level in TAXON_LEVELS:
            own = getattr(tax, level)
            top, other = top_code(own), rest_code(totals[level], own)
            row[f"tx_{level}_share"] = share_of(own.get(top, 0), tax.n_ncldv)
            row[f"tx_{level}_agree"] = (
                share_of(own.get(other, 0), tax.n_ncldv) if other else float("nan")
            )
            row[f"tx_{level}_differs"] = (
                float(top != other) if top and other else float("nan")
            )
        for name in TAXON_CLASSES:
            mine = share_of(tax.classes.get(name, 0), tax.n_hit)
            rest = share_of(
                class_total[name] - tax.classes.get(name, 0), hit_total - tax.n_hit
            )
            row[f"tx_cls_{name}"] = mine
            row[f"tx_cls_{name}_diff"] = mine - rest
        rows.append(row)
    return pd.DataFrame(rows)


def with_taxonomy(frame: pd.DataFrame, contigs: list[Contig]) -> pd.DataFrame:
    """Return the features with the hit-taxonomy columns when every contig has them."""
    if any(c.tax is None for c in contigs):
        return frame
    return pd.concat([frame, taxonomy_features(contigs)], axis=1)


class IdentityCache:
    """Amino-acid identity of two proteins (pyswrd, BLOSUM62), cached by hash pair."""

    def __init__(self) -> None:
        """Start empty."""
        self.table: dict[tuple[str, str], float] = {}

    def fill(self, items: list[tuple[str, str]]) -> None:
        """Compute identities for all pairs among (hash, protein) items."""
        todo = {}
        for (ha, _), (hb, _) in itertools.combinations(items, 2):
            key = (min(ha, hb), max(ha, hb))
            if ha != hb and key not in self.table:
                todo[key] = None
        if not todo:
            return
        sequences = dict(items)
        hashes = sorted({h for key in todo for h in key})
        proteins = [sequences[h].replace("*", "") for h in hashes]
        position = {h: i for i, h in enumerate(hashes)}
        found: dict[tuple[int, int], float] = {}
        for hit in pyswrd.search(
            proteins,
            proteins,
            scorer_name="BLOSUM62",
            score_threshold=0,
            max_alignments=len(proteins),
            threads=1,
        ):
            if hit.query_index != hit.target_index:
                found[(hit.query_index, hit.target_index)] = float(
                    hit.result.identity()
                )
        for key in todo:
            i, j = position[key[0]], position[key[1]]
            self.table[key] = max(found.get((i, j), 0.0), found.get((j, i), 0.0))

    def get(self, a: str, b: str) -> float:
        """Return the identity of two proteins (1.0 for the same hash)."""
        if a == b:
            return 1.0
        return self.table.get((min(a, b), max(a, b)), 0.0)


def pair_keys(models: list[str]) -> set[tuple[str, str]]:
    """Return unordered marker pairs at most ``SKIP`` markers apart."""
    keys = set()
    for i, first in enumerate(models):
        for second in models[i + 1 : i + 1 + SKIP]:
            keys.add((min(first, second), max(first, second)))
    return keys


def duplication_stats(
    index: int,
    contig: Contig,
    by_model: dict[str, list[tuple[int, str, str]]],
    priors: Priors,
    cache: IdentityCache,
) -> tuple[dict[str, float], dict[int, int]]:
    """Return copy-number features of one contig and its partner counts."""
    n_g = len(contig.models)
    identities: list[float] = []
    weighted = 0.0
    partners: dict[int, int] = {}
    n_sc = n_dup_sc = 0
    for model, digest in zip(contig.models, contig.hashes, strict=True):
        single = priors.single_copy(model)
        n_sc += single
        others = [(i, h) for i, h, _ in by_model[model] if i != index]
        if not others:
            continue
        n_dup_sc += single
        weighted += priors.weight(model)
        identities.append(max(cache.get(digest, h) for _, h in others))
        for other, _ in others:
            partners[other] = partners.get(other, 0) + 1
    arr = np.array(identities)

    def share(mask: np.ndarray) -> float:
        return float(mask.sum()) / n_g if n_g else np.nan

    return (
        {
            "dup_frac": len(identities) / n_g if n_g else np.nan,
            "dup_weighted": weighted / n_g if n_g else np.nan,
            "dup60_frac": share(arr >= 0.6),
            "dup80_frac": share(arr >= 0.8),
            "dup95_frac": share(arr >= 0.95),
            "dup_identity_mean": float(arr.mean()) if len(arr) else np.nan,
            "dup_identity_max": float(arr.max()) if len(arr) else np.nan,
            "n_sc": n_sc,
            "n_dup_sc": n_dup_sc,
            "dup_sc_frac": n_dup_sc / n_sc if n_sc else np.nan,
        },
        partners,
    )


def synteny_stats(
    index: int,
    keys: list[set[tuple[str, str]]],
    key_owner: dict[tuple[str, str], set[int]],
    length: np.ndarray,
) -> dict[str, float]:
    """Return the share of the contig's marker pairs found on other contigs."""
    shared = larger = 0
    for key in keys[index]:
        owners = key_owner[key] - {index}
        if owners:
            shared += 1
            larger += any(length[o] > length[index] for o in owners)
    denom = len(keys[index])
    return {
        "pair_shared": shared / denom if denom else np.nan,
        "pair_shared_larger": larger / denom if denom else np.nan,
        "n_pairs": denom,
    }


def bulk_distance(index: int, matrix: np.ndarray, length: np.ndarray) -> float:
    """Return the distance of one row to the length-weighted bulk of the others."""
    if len(length) < 2:
        return float("nan")
    weights = np.delete(length, index)
    bulk = (np.delete(matrix, index, axis=0) * weights[:, None]).sum(
        axis=0
    ) / weights.sum()
    return float(np.linalg.norm(matrix[index] - bulk))


def composition_stats(
    index: int,
    standard: np.ndarray,
    gc: np.ndarray,
    length: np.ndarray,
    partner: int,
) -> dict[str, float]:
    """Return composition distances of one contig to the rest of the bin."""
    if len(length) < 2:
        return {"d_bulk": np.nan, "d_partner": np.nan, "gc_diff": np.nan}
    weights = np.delete(length, index)
    bulk = (np.delete(standard, index, axis=0) * weights[:, None]).sum(
        axis=0
    ) / weights.sum()
    d_partner = (
        float(np.linalg.norm(standard[index] - standard[partner]))
        if partner >= 0
        else np.nan
    )
    return {
        "d_bulk": float(np.linalg.norm(standard[index] - bulk)),
        "d_partner": d_partner,
        "gc_diff": float(
            abs(gc[index] - np.average(np.delete(gc, index), weights=weights))
        ),
    }


def bin_features(
    contigs: list[Contig], priors: Priors, cache: IdentityCache
) -> pd.DataFrame:
    """Return one feature row per contig of a bin."""
    by_model: dict[str, list[tuple[int, str, str]]] = {}
    for index, contig in enumerate(contigs):
        for model, protein, digest in zip(
            contig.models, contig.proteins, contig.hashes, strict=True
        ):
            by_model.setdefault(model, []).append((index, digest, protein))
    for items in by_model.values():
        if len({i for i, _, _ in items}) > 1:
            cache.fill([(h, p) for _, h, p in items])
    keys = [pair_keys(c.models) for c in contigs]
    key_owner: dict[tuple[str, str], set[int]] = {}
    for index, owned in enumerate(keys):
        for key in owned:
            key_owner.setdefault(key, set()).add(index)
    length = np.array([c.length for c in contigs], dtype=float)
    comp = np.array([c.comp for c in contigs])
    standard = (comp - priors.mu) / priors.sd
    extra = {}
    for key, (mu, sd) in priors.extra.items():
        raw = np.array([getattr(c, key) for c in contigs])
        extra[key] = (raw - mu) / sd
    n_markers = sum(len(c.models) for c in contigs)
    rows = []
    partner_counts: list[dict[int, int]] = []
    for index, contig in enumerate(contigs):
        dup, partners = duplication_stats(index, contig, by_model, priors, cache)
        partner_counts.append(partners)
        partner = max(partners, key=partners.get) if partners else -1
        n_g = len(contig.models)
        row = {
            "contig": contig.name,
            "genome": contig.genome,
            "length_bp": contig.length,
            "log_len": math.log10(contig.length),
            "rel_len": contig.length / float(length.sum()),
            "len_rank": float((length > contig.length).sum())
            / max(len(contigs) - 1, 1),
            "n_markers": n_g,
            "marker_density": n_g / (contig.length / 1000.0),
            "partner_rel_len": float(length[partner] / contig.length)
            if partner >= 0
            else np.nan,
            "n_contigs": len(contigs),
            "bin_bp": float(length.sum()),
            "bin_markers": n_markers,
        }
        row.update(dup)
        row.update(synteny_stats(index, keys, key_owner, length))
        row.update(composition_stats(index, standard, comp[:, -1], length, partner))
        for key, matrix in extra.items():
            row[f"d_{key}"] = bulk_distance(index, matrix, length)
        row["d_partner_ratio"] = (
            row["d_partner"] / row["d_bulk"]
            if partner >= 0 and row["d_bulk"]
            else np.nan
        )
        rows.append(row)
    frame = pd.DataFrame(rows)
    partner_of = [max(p, key=p.get) if p else -1 for p in partner_counts]
    d_bulk = frame["d_bulk"].to_numpy()
    frame["d_bulk_partner"] = [d_bulk[q] if q >= 0 else np.nan for q in partner_of]
    frame["d_diff"] = frame["d_bulk"] - frame["d_bulk_partner"]
    duplicated = (frame["dup_frac"].fillna(0) * frame["n_markers"]).sum()
    frame["bin_n_sc"] = sum(frame["n_sc"])
    frame["bin_dup_sc_frac"] = frame["n_dup_sc"].sum() / max(frame["n_sc"].sum(), 1)
    frame["bin_dup_frac"] = duplicated / max(n_markers, 1)
    core = {m for c in contigs for m in c.models if m.startswith("GVOGm")}
    multi = sum(len(by_model[m]) >= 2 for m in core)
    frame["bin_core_multi"] = multi / max(len(core), 1)
    return with_taxonomy(frame, contigs)
