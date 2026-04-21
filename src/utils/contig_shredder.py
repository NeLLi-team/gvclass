"""Realistic MAG-like contig shredding.

Simulate metagenome-assembled-genome (MAG) fragmentation of reference isolate
genomes so trained contamination models see contig-length distributions
comparable to real binned assemblies.

Design
------
Real MAG contig length distributions (Tara Oceans GEM-MAGs, GORG, etc.)
approximate a log-normal shape with:

* assembler minimum-contig floor ~ 1.5-2 kb
* median ~ 5-30 kb
* long tail to 100-300 kb

We draw fragment counts and lengths from a log-normal(mu, sigma) model,
renormalise the lengths so they tile the source bp exactly (no deletion,
no overlap), and cut the source contig at those positions. Fragment order
may be shuffled to mimic binning randomness.

Determinism
-----------
All randomness flows through a seeded ``numpy.random.Generator``. Callers
are expected to derive the seed from a stable key such as
``hash(base_accession + scenario + variant_idx)`` so re-runs produce
byte-identical FASTAs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


DEFAULT_MIN_FRAGMENT_BP = 2_000
DEFAULT_MAX_FRAGMENT_BP = 200_000
DEFAULT_MEAN_FRAGMENT_BP = 20_000
DEFAULT_SIGMA = 0.7
DEFAULT_LONG_FRAGMENT_MIN_BP = 30_000


@dataclass(frozen=True)
class ShredConfig:
    """Parameters controlling log-normal shredding.

    Attributes
    ----------
    mean_fragment_bp:
        Target mean fragment length in bp. Drives the expected fragment count
        ``N ~ genome_bp / mean_fragment_bp``.
    sigma:
        Log-normal scale parameter. ``0.7`` gives a realistic MAG-like spread
        (coefficient of variation ~ 0.8).
    min_fragment_bp:
        Hard floor on fragment length, matching common assembler
        ``min_contig_len`` cutoffs.
    max_fragment_bp:
        Hard ceiling; caps very rare extreme draws.
    long_fragment_min_bp:
        Guarantee that at least one fragment meets this threshold. Real MAGs
        almost always contain a handful of long contigs on which HMM markers
        can land; without this floor the benchmark can become pathologically
        easy/hard depending on draw.
    shuffle:
        If True, randomise output fragment order (mimics binning). If False,
        fragments appear in genome order.
    """

    mean_fragment_bp: int = DEFAULT_MEAN_FRAGMENT_BP
    sigma: float = DEFAULT_SIGMA
    min_fragment_bp: int = DEFAULT_MIN_FRAGMENT_BP
    max_fragment_bp: int = DEFAULT_MAX_FRAGMENT_BP
    long_fragment_min_bp: int = DEFAULT_LONG_FRAGMENT_MIN_BP
    shuffle: bool = True

    def __post_init__(self) -> None:
        if self.mean_fragment_bp <= 0:
            raise ValueError("mean_fragment_bp must be positive")
        if self.sigma <= 0:
            raise ValueError("sigma must be positive")
        if self.min_fragment_bp < 1:
            raise ValueError("min_fragment_bp must be >= 1")
        if self.max_fragment_bp <= self.min_fragment_bp:
            raise ValueError("max_fragment_bp must exceed min_fragment_bp")
        if self.long_fragment_min_bp < self.min_fragment_bp:
            raise ValueError(
                "long_fragment_min_bp must be >= min_fragment_bp"
            )


def _draw_fragment_lengths(
    genome_bp: int,
    config: ShredConfig,
    rng: np.random.Generator,
) -> List[int]:
    """Draw and renormalise fragment lengths that tile ``genome_bp``.

    Steps:

    1. Target count ``N = max(1, round(genome_bp / mean_fragment_bp))``.
    2. Draw ``N`` log-normal samples, clip to [min, max].
    3. Scale lengths so they sum exactly to ``genome_bp``. Post-scaling,
       re-clip any fragment that fell below ``min_fragment_bp`` by merging
       it with a neighbour.
    4. Fix rounding residue by adjusting the last fragment.
    """
    if genome_bp < config.min_fragment_bp:
        return [genome_bp]

    target_count = max(1, round(genome_bp / config.mean_fragment_bp))

    mu = np.log(config.mean_fragment_bp) - 0.5 * config.sigma**2
    raw = rng.lognormal(mean=mu, sigma=config.sigma, size=target_count)
    raw = np.clip(raw, config.min_fragment_bp, config.max_fragment_bp)

    scale = genome_bp / raw.sum()
    scaled = raw * scale

    lengths = np.round(scaled).astype(int).tolist()

    # Fix rounding residue by adjusting the final fragment.
    residue = genome_bp - sum(lengths)
    lengths[-1] += residue

    # Merge any fragment below the min floor into its neighbour. This can
    # happen after scaling when many fragments were drawn for a small
    # genome. Iterative merge keeps floor invariant.
    merged: List[int] = []
    carry = 0
    for length in lengths:
        combined = length + carry
        if combined < config.min_fragment_bp:
            carry = combined
            continue
        merged.append(combined)
        carry = 0
    if carry > 0:
        if merged:
            merged[-1] += carry
        else:
            merged.append(carry)

    # Enforce the long-fragment guarantee. If no fragment meets the
    # threshold, merge the two largest until one does. This can only
    # fail if genome_bp < long_fragment_min_bp, in which case we return
    # the single-fragment result above.
    while merged and max(merged) < config.long_fragment_min_bp:
        if sum(merged) < config.long_fragment_min_bp:
            break
        merged.sort(reverse=True)
        big = merged.pop(0) + merged.pop(0)
        merged.append(big)

    # Sanity: preserve total bp.
    assert sum(merged) == genome_bp, (
        f"shredded bp {sum(merged)} != source {genome_bp}"
    )
    return merged


def shred_sequence(
    sequence: str,
    config: ShredConfig,
    rng: np.random.Generator,
) -> List[str]:
    """Cut a single DNA string into MAG-like fragments.

    Cuts are contiguous and non-overlapping; total bp is preserved.
    """
    if not sequence:
        return []
    lengths = _draw_fragment_lengths(len(sequence), config, rng)

    fragments: List[str] = []
    cursor = 0
    for length in lengths:
        fragments.append(sequence[cursor : cursor + length])
        cursor += length

    if config.shuffle and len(fragments) > 1:
        order = rng.permutation(len(fragments))
        fragments = [fragments[i] for i in order]

    return fragments


def shred_record(
    record: SeqRecord,
    config: ShredConfig,
    rng: np.random.Generator,
    id_prefix: str | None = None,
) -> List[SeqRecord]:
    """Shred one ``SeqRecord`` into MAG-like fragments.

    Fragments are renamed ``{id_prefix or record.id}_frag_{i:04d}``, preserving
    the source id as the prefix so provenance stays traceable.
    """
    sequence = str(record.seq).upper()
    base_id = id_prefix if id_prefix is not None else record.id
    fragments = shred_sequence(sequence, config, rng)

    out: List[SeqRecord] = []
    for idx, frag in enumerate(fragments):
        frag_record = SeqRecord(
            Seq(frag),
            id=f"{base_id}_frag_{idx:04d}",
            description="",
            name="",
        )
        out.append(frag_record)
    return out


def shred_records(
    records: Iterable[SeqRecord],
    config: ShredConfig,
    rng: np.random.Generator,
    id_prefix: str | None = None,
) -> List[SeqRecord]:
    """Shred an iterable of ``SeqRecord`` objects; fragments are
    concatenated in source-record order then (optionally) shuffled
    globally so fragments from different source contigs interleave."""
    fragments: List[SeqRecord] = []
    for record in records:
        # Disable per-record shuffle so we can control shuffle globally.
        per_record_cfg = ShredConfig(
            mean_fragment_bp=config.mean_fragment_bp,
            sigma=config.sigma,
            min_fragment_bp=config.min_fragment_bp,
            max_fragment_bp=config.max_fragment_bp,
            long_fragment_min_bp=config.long_fragment_min_bp,
            shuffle=False,
        )
        prefix = id_prefix if id_prefix is not None else record.id
        fragments.extend(shred_record(record, per_record_cfg, rng, id_prefix=prefix))

    if config.shuffle and len(fragments) > 1:
        order = rng.permutation(len(fragments))
        fragments = [fragments[i] for i in order]

    # Rename globally so contig ids stay unique even after shuffle.
    renamed: List[SeqRecord] = []
    for idx, frag in enumerate(fragments):
        renamed.append(
            SeqRecord(
                frag.seq,
                id=f"{frag.id}__g{idx:04d}",
                description="",
                name="",
            )
        )
    return renamed


def fragment_stats(records: Sequence[SeqRecord]) -> Tuple[int, int, int, int, float]:
    """Return ``(count, min_bp, max_bp, total_bp, mean_bp)`` for a fragment set."""
    if not records:
        return (0, 0, 0, 0, 0.0)
    lengths = [len(r.seq) for r in records]
    return (
        len(lengths),
        min(lengths),
        max(lengths),
        sum(lengths),
        sum(lengths) / len(lengths),
    )


def make_rng(seed_key: str) -> np.random.Generator:
    """Derive a reproducible ``Generator`` from a string key.

    Uses ``numpy``'s ``SeedSequence`` on a 64-bit hash so identical keys
    produce identical streams across runs and platforms. Python's built-in
    ``hash()`` is salted per-process and must not be used here.
    """
    import hashlib

    digest = hashlib.sha256(seed_key.encode("utf-8")).digest()
    seed_int = int.from_bytes(digest[:8], "little", signed=False)
    return np.random.default_rng(np.random.SeedSequence(seed_int))
