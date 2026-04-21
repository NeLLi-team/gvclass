"""Unit tests for realistic MAG-like contig shredding."""

from __future__ import annotations

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.utils.contig_shredder import (
    ShredConfig,
    fragment_stats,
    make_rng,
    shred_record,
    shred_records,
    shred_sequence,
)


def _random_seq(length: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=length))


def test_shred_sequence_preserves_total_bp():
    seq = _random_seq(500_000, seed=1)
    cfg = ShredConfig()
    rng = make_rng("test_total_bp")
    frags = shred_sequence(seq, cfg, rng)
    assert sum(len(f) for f in frags) == len(seq)


def test_shred_sequence_fragment_count_reasonable():
    seq = _random_seq(400_000, seed=2)
    cfg = ShredConfig(mean_fragment_bp=20_000)
    rng = make_rng("test_count")
    frags = shred_sequence(seq, cfg, rng)
    # target ~ 20 fragments (400k / 20k); allow wide tolerance for log-normal spread.
    assert 5 <= len(frags) <= 80


def test_shred_sequence_respects_min_floor():
    seq = _random_seq(300_000, seed=3)
    cfg = ShredConfig(min_fragment_bp=2_000)
    rng = make_rng("test_min")
    frags = shred_sequence(seq, cfg, rng)
    assert all(len(f) >= cfg.min_fragment_bp for f in frags)


def test_shred_sequence_guarantees_long_fragment():
    """Every shred of a reasonably sized genome must contain >=1 fragment
    at or above the long-fragment threshold."""
    seq = _random_seq(500_000, seed=4)
    cfg = ShredConfig(long_fragment_min_bp=30_000)
    for i in range(20):
        rng = make_rng(f"test_long_{i}")
        frags = shred_sequence(seq, cfg, rng)
        assert max(len(f) for f in frags) >= cfg.long_fragment_min_bp, (
            f"iteration {i}: longest fragment {max(len(f) for f in frags)} "
            f"below threshold {cfg.long_fragment_min_bp}"
        )


def test_shred_sequence_deterministic():
    """Same seed key yields identical fragment sequence."""
    seq = _random_seq(200_000, seed=5)
    cfg = ShredConfig()
    a = shred_sequence(seq, cfg, make_rng("dtm"))
    b = shred_sequence(seq, cfg, make_rng("dtm"))
    assert a == b


def test_shred_sequence_different_seeds_differ():
    seq = _random_seq(200_000, seed=6)
    cfg = ShredConfig()
    a = shred_sequence(seq, cfg, make_rng("seed_a"))
    b = shred_sequence(seq, cfg, make_rng("seed_b"))
    assert a != b


def test_shred_sequence_median_length_realistic():
    """Aggregate across many genomes: median fragment length should sit in
    the realistic MAG range (roughly 5-35 kb for the default config)."""
    cfg = ShredConfig(mean_fragment_bp=20_000)
    all_lengths: list[int] = []
    for i in range(10):
        seq = _random_seq(500_000, seed=100 + i)
        rng = make_rng(f"agg_{i}")
        all_lengths.extend(len(f) for f in shred_sequence(seq, cfg, rng))
    median = int(np.median(all_lengths))
    assert 5_000 <= median <= 35_000, f"median {median} outside realistic range"


def test_shred_sequence_short_genome_single_fragment():
    seq = _random_seq(1_500, seed=7)
    cfg = ShredConfig(min_fragment_bp=2_000)
    rng = make_rng("short")
    frags = shred_sequence(seq, cfg, rng)
    assert frags == [seq]


def test_shred_record_preserves_bp_and_ids():
    seq = _random_seq(400_000, seed=8)
    record = SeqRecord(Seq(seq), id="base_contig_A", description="")
    cfg = ShredConfig()
    rng = make_rng("record")
    frags = shred_record(record, cfg, rng)
    assert sum(len(f.seq) for f in frags) == len(seq)
    assert all(f.id.startswith("base_contig_A_frag_") for f in frags)
    assert len({f.id for f in frags}) == len(frags)


def test_shred_records_concatenates_and_renames_uniquely():
    rec_a = SeqRecord(Seq(_random_seq(200_000, seed=9)), id="A", description="")
    rec_b = SeqRecord(Seq(_random_seq(200_000, seed=10)), id="B", description="")
    cfg = ShredConfig()
    rng = make_rng("multi")
    frags = shred_records([rec_a, rec_b], cfg, rng)
    ids = [f.id for f in frags]
    assert len(ids) == len(set(ids)), "fragment ids must be globally unique"
    # Total bp is preserved across both sources.
    assert sum(len(f.seq) for f in frags) == 400_000


def test_shred_records_shuffles_across_sources_when_enabled():
    """With shuffle=True the global order should not be strictly grouped
    by source contig across at least one of several runs."""
    rec_a = SeqRecord(Seq(_random_seq(200_000, seed=11)), id="A", description="")
    rec_b = SeqRecord(Seq(_random_seq(200_000, seed=12)), id="B", description="")
    cfg = ShredConfig(shuffle=True)

    interleaved = False
    for i in range(5):
        rng = make_rng(f"shuffle_{i}")
        frags = shred_records([rec_a, rec_b], cfg, rng)
        prefixes = [f.id.split("_frag_")[0] for f in frags]
        # Find at least one A after a B or vice versa -> interleaved.
        if any(prefixes[j] != prefixes[j + 1] for j in range(len(prefixes) - 1)):
            interleaved = True
            break
    assert interleaved, "shuffle=True failed to interleave across sources"


def test_shred_config_validates():
    with pytest.raises(ValueError):
        ShredConfig(mean_fragment_bp=0)
    with pytest.raises(ValueError):
        ShredConfig(sigma=0)
    with pytest.raises(ValueError):
        ShredConfig(min_fragment_bp=0)
    with pytest.raises(ValueError):
        ShredConfig(min_fragment_bp=5_000, max_fragment_bp=1_000)
    with pytest.raises(ValueError):
        ShredConfig(min_fragment_bp=5_000, long_fragment_min_bp=1_000)


def test_fragment_stats():
    cfg = ShredConfig()
    rng = make_rng("stats")
    seq = _random_seq(200_000, seed=13)
    record = SeqRecord(Seq(seq), id="X")
    frags = shred_record(record, cfg, rng)
    count, mn, mx, total, mean = fragment_stats(frags)
    assert count == len(frags)
    assert total == len(seq)
    assert mn >= cfg.min_fragment_bp
    assert mx <= cfg.max_fragment_bp
    assert mean == pytest.approx(total / count)


def test_make_rng_deterministic_across_calls():
    a = make_rng("same_key").random(5)
    b = make_rng("same_key").random(5)
    np.testing.assert_array_equal(a, b)
    c = make_rng("other_key").random(5)
    assert not np.array_equal(a, c)
