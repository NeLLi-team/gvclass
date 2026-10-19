"""Find contigs of another giant virus in a giant-virus bin.

Usage: python -m src.contig_screen.giant_virus <bin.fna> --resources DIR --out out.tsv
                                               [--threads N] [--cellular cellular.tsv]

Calls genes on the bin, assigns the catalog markers (OG, GVOGm), and scores every
contig of 5 kb or more for the evidence of a second genome: markers that occur on
other contigs although they are single-copy in reference genomes, shared neighbour
pairs of markers, identity of the copies, the composition of the contig against the
rest of the bin, and the best catalog hits of its genes (DIAMOND against ``ref7x``)
with their agreement with the bin (``gv_order``, ``gv_taxa``). The output has one row
per contig with P(foreign) and a call at the threshold ``0.01`` of ``thresholds.json``
(1% of clean contigs called), unless ``--threshold`` is given. With ``--cellular``
(the output of ``cellular``) the cellular calls are joined to the giant-virus calls
and the union is reported as well. A call is a review flag, not a contamination
estimate.
"""

from __future__ import annotations

import argparse
import logging
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.contig_screen import layout, torch_stack
from src.contig_screen.gv_genes import call_genes, read_fasta
from src.contig_screen.gv_model import (
    OPERATING_POINT,
    Bundle,
    load_bundle,
    load_thresholds,
)
from src.contig_screen.gv_order import IdentityCache, bin_features
from src.contig_screen.gv_order_build import Genome
from src.contig_screen.gv_taxa import taxa_for_bin
from src.contig_screen.thresholds import sliding_calls

logger = logging.getLogger(__name__)

MIN_SCORED_BP = 5000
KEPT_COLUMNS = [
    "contig",
    "length_bp",
    "n_markers",
    "n_dup_sc",
    "dup_weighted",
    "pair_shared",
    "d_bulk",
    "p_foreign",
    "called_foreign",
]


@dataclass(frozen=True, slots=True)
class Settings:
    """How many threads the searches use and which DIAMOND binary runs."""

    threads: int = 4
    diamond: Path | None = None


def build_parser() -> argparse.ArgumentParser:
    """Return the command line parser."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fasta", type=Path, help="bin (FASTA, optionally gzipped)")
    parser.add_argument(
        "--resources",
        type=Path,
        required=True,
        help="resource bundle root (holds contig_screen/, hmm/)",
    )
    parser.add_argument("--out", type=Path, required=True, help="output table (TSV)")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help=f"P(foreign) threshold (default: thresholds.json key {OPERATING_POINT})",
    )
    parser.add_argument(
        "--cellular",
        type=Path,
        default=None,
        help="output table of the cellular screen, to report the union of the calls",
    )
    parser.add_argument(
        "--cellular-threshold",
        type=float,
        default=None,
        help="stage 2 threshold of the cellular calls (default: thresholds.json)",
    )
    parser.add_argument(
        "--diamond", type=Path, default=None, help="DIAMOND binary (default: on PATH)"
    )
    parser.add_argument(
        "--manifest-sha256",
        default=None,
        help="digest of manifest.json, for a bundle other than the release bundle",
    )
    parser.add_argument(
        "--verify-all",
        action="store_true",
        help="hash every file of the bundle, DIAMOND databases included",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Score one bin, write the table and print the summary."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    try:
        resources = layout.ScreenResources.open(
            args.resources,
            manifest_sha256=args.manifest_sha256,
            verify_all=args.verify_all,
        )
        bundle = load_bundle(resources)
        cuts = load_thresholds(resources)
    except (layout.IntegrityError, ValueError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    threshold = cuts[OPERATING_POINT] if args.threshold is None else args.threshold
    sequences = read_fasta(args.fasta)
    lengths = {name: len(sequence) for name, sequence in sequences.items()}
    total = float(sum(lengths.values()))
    foreign: set[str] = set()
    with torch_stack.limit_threads(args.threads):
        table = score_bin(
            sequences, resources, bundle, Settings(args.threads, args.diamond)
        )
    if table.empty:
        print("fewer than two contigs of 5 kb or more with genes: nothing to compare")
    else:
        table["called_foreign"] = table["p_foreign"] >= threshold
        args.out.parent.mkdir(parents=True, exist_ok=True)
        table[KEPT_COLUMNS].to_csv(args.out, sep="\t", index=False)
        foreign = summarize(table, threshold, total)
    if args.cellular is not None:
        cellular_threshold = args.cellular_threshold
        if cellular_threshold is None and resources.has(layout.THRESHOLDS_FILE):
            cellular_threshold = layout.load_thresholds(resources)["stage2"]
        report_union(args.cellular, cellular_threshold, lengths, foreign, total)
    return 0


def score_bin(
    sequences: Mapping[str, str],
    resources: layout.ScreenResources,
    bundle: Bundle,
    settings: Settings,
) -> pd.DataFrame:
    """Return the contig table of one bin with P(foreign).

    The table is empty when the bin has fewer than two contigs of 5 kb or more with
    genes.
    """
    calls = call_genes(sequences, resources.path(layout.HMM_FILE), settings.threads)
    taxa = None
    if any(column.startswith("tx_") for column in bundle.columns):
        with tempfile.TemporaryDirectory() as scratch:
            taxa = taxa_for_bin(
                calls.proteins,
                resources,
                Path(scratch),
                threads=settings.threads,
                diamond=settings.diamond,
            )
    contigs = Genome("bin", sequences, calls.table, taxa).contigs()
    if len(contigs) < 2:
        return pd.DataFrame()
    table = bin_features(contigs, bundle.priors, IdentityCache())
    for column in bundle.columns:
        if column not in table:
            table[column] = np.nan
    features = table[bundle.columns].to_numpy(dtype=float)
    table["p_foreign"] = bundle.model.predict_proba(features)[:, 1]
    return table


def summarize(table: pd.DataFrame, threshold: float, total: float) -> set[str]:
    """Print the called share of the input; return the called contigs."""
    scored = table[table["length_bp"] >= MIN_SCORED_BP]
    foreign_bp = float(scored.loc[scored["called_foreign"], "length_bp"].sum())
    print(f"contigs scored: {len(scored)} of {len(table)} (5 kb or more)")
    print(
        f"foreign giant-virus contigs: {int(scored['called_foreign'].sum())}, "
        f"{foreign_bp:.0f} bp = {100 * foreign_bp / total:.1f}% of the input "
        f"(threshold {threshold:.3f})"
    )
    return set(scored.loc[scored["called_foreign"], "contig"])


def report_union(
    cellular: Path,
    threshold: float | None,
    lengths: Mapping[str, int],
    foreign: set[str],
    total: float,
) -> None:
    """Print the cellular calls of ``--cellular`` and their union with contigs."""
    if threshold is None:
        print("--cellular needs --cellular-threshold; union not computed")
        return
    windows = pd.read_csv(cellular, sep="\t")
    score = "p_context" if "p_context" in windows else "p_cellular"
    calls = sliding_calls(windows.assign(p=windows[score]), dict(lengths), threshold)
    cellular_bp = union_bp = 0
    for contig, mask in calls.items():
        cellular_bp += int(mask.sum())
        union_bp += len(mask) if contig in foreign else int(mask.sum())
    for contig in foreign - set(calls):
        union_bp += lengths[contig]
    print(
        f"cellular calls: {cellular_bp} bp = {100 * cellular_bp / total:.1f}%; "
        f"union with giant-virus contigs: {union_bp} bp = {100 * union_bp / total:.1f}%"
    )


if __name__ == "__main__":
    raise SystemExit(main())
