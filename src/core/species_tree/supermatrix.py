"""Supermatrix construction for the species-tree feature (Section 4).

Aligns each GVOG8 group's representative proteins, trims alignment columns
(witchi chi-squared pruning, with pytrimal and then unpruned fallbacks),
namespaces leaves reversibly, and concatenates per taxon into a gap-padded
supermatrix plus a partition file.

Leaf namespacing keeps queries and references distinct and reversible even if a
query file stem itself begins with a reference domain prefix: queries become
``QUERY__<stem>`` while references keep their already-domain-prefixed genome id
(``NCLDV__<genome>``). :func:`denamespace` reverses it for output.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Set

import pytrimal
from Bio import SeqIO

from src.core.alignment import align_sequences_pyfamsa

logger = logging.getLogger(__name__)

QUERY_PREFIX = "QUERY__"


def namespace_query(stem: str) -> str:
    """Namespace a query file stem as a supermatrix leaf."""
    return f"{QUERY_PREFIX}{stem}"


def is_query_leaf(name: str) -> bool:
    return name.startswith(QUERY_PREFIX)


def denamespace(name: str) -> str:
    """Reverse leaf namespacing for output (queries lose the prefix; refs keep id)."""
    if name.startswith(QUERY_PREFIX):
        return name[len(QUERY_PREFIX):]
    return name


@dataclass
class SupermatrixResult:
    supermatrix_faa: Path
    partitions: Path
    group_widths: Dict[str, int]
    included_groups: List[str]
    leaves: List[str]
    dropped_groups: List[str] = field(default_factory=list)
    dropped_leaves: List[str] = field(default_factory=list)


def _read_fasta_dict(path: Path) -> Dict[str, str]:
    return {record.id: str(record.seq) for record in SeqIO.parse(str(path), "fasta")}


def _uniform_nonempty(block: Dict[str, str]) -> bool:
    if not block:
        return False
    widths = {len(seq) for seq in block.values()}
    return len(widths) == 1 and next(iter(widths)) > 0


def _witchi_prune(aln_path: Path, threads: int, timeout: int = 1200) -> Optional[Path]:
    """Run the witchi prune CLI on ``aln_path``; return the pruned fasta or None.

    witchi writes ``<stem>_<algorithm>_s<top_n>_pruned.fasta`` next to the input
    (Section 0), so the output is located by globbing rather than assuming a name.
    """
    witchi = shutil.which("witchi")
    if witchi is None:
        logger.warning("witchi not on PATH; falling back to pytrimal")
        return None
    try:
        proc = subprocess.run(
            [
                witchi, "prune",
                "--file", str(aln_path),
                "--format", "fasta",
                "--num_workers_chisq", str(max(1, threads)),
                "--num_workers_permute", str(max(1, threads)),
                # Speed: 20 permutations (vs default 100) for the stopping null and
                # 10 columns removed per iteration (vs 1) cut the recursive pruning
                # cost ~50x on the large supermatrix blocks with negligible effect
                # on the kept-column set. Output name becomes <stem>_wasserstein_s10
                # _pruned.fasta, still matched by the <stem>_*_pruned.fasta glob.
                "--permutations", "20",
                "--top_n", "10",
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except Exception as exc:  # subprocess/timeout guard
        logger.warning("witchi prune errored for %s: %s", aln_path.name, exc)
        return None
    if proc.returncode != 0:
        logger.warning(
            "witchi prune nonzero exit for %s: %s", aln_path.name, proc.stderr[-300:]
        )
        return None

    pruned = sorted(aln_path.parent.glob(f"{aln_path.stem}_*_pruned.fasta"))
    if not pruned:
        logger.warning("witchi produced no *_pruned.fasta for %s", aln_path.name)
        return None
    return pruned[0]


def _pytrimal_trim(aln_path: Path, out_path: Path) -> Path:
    alignment = pytrimal.Alignment.load(str(aln_path))
    trimmed = pytrimal.AutomaticTrimmer(method="automated1").trim(alignment)
    trimmed.dump(str(out_path))
    return out_path


def _valid_block(block: Dict[str, str], expected_leaves: Set[str]) -> bool:
    """A trim result is usable only if it covers EXACTLY the input taxa with a
    uniform, non-empty width. The taxon-set check guards against a stale or
    partial pruned file (uniform width alone is not enough)."""
    return _uniform_nonempty(block) and set(block) == expected_leaves


def _trim_group_alignment(
    aln_path: Path,
    threads: int,
    expected_leaves: Set[str],
    method: str = "witchi",
) -> Dict[str, str]:
    """Trim a group alignment by ``method``, with fallbacks. Returns {leaf: seq}.

    * ``"witchi"`` (default): witchi prune -> pytrimal -> unpruned.
    * ``"pytrimal"``: pytrimal -> unpruned (skip the slower witchi step).
    * ``"none"``: unpruned (raw alignment).

    Each candidate is accepted only if it covers exactly ``expected_leaves`` with a
    uniform, non-empty width, so a stale/partial pruned file or an over-aggressive
    0-column trim can never silently corrupt or empty a block. The unpruned raw
    alignment is the always-valid terminal fallback.
    """
    # Remove stale pruning artefacts first so this run's output is discovered
    # unambiguously even when work_dir is reused (or the method changed) across runs.
    for stale in aln_path.parent.glob(f"{aln_path.stem}_*_pruned.fasta"):
        stale.unlink()
    trimal_out = aln_path.parent / f"{aln_path.stem}.trimal.fasta"
    if trimal_out.exists():
        trimal_out.unlink()

    if method == "none":
        return _read_fasta_dict(aln_path)

    if method == "witchi":
        pruned = _witchi_prune(aln_path, threads)
        if pruned is not None:
            block = _read_fasta_dict(pruned)
            if _valid_block(block, expected_leaves):
                return block
            logger.warning(
                "witchi output for %s unusable; trying pytrimal", aln_path.name
            )

    try:
        _pytrimal_trim(aln_path, trimal_out)
        block = _read_fasta_dict(trimal_out)
        if _valid_block(block, expected_leaves):
            return block
        logger.warning("pytrimal output for %s unusable; using unpruned", aln_path.name)
    except Exception as exc:
        logger.warning("pytrimal fallback failed for %s: %s", aln_path.name, exc)

    # Unpruned raw alignment trivially covers exactly the expected taxa.
    return _read_fasta_dict(aln_path)


def build_supermatrix(
    taxa_group_seqs: Dict[str, Dict[str, str]],
    groups: Sequence[str],
    work_dir: Path,
    out_dir: Path,
    basename: str = "combined",
    threads: int = 4,
    min_taxa_per_group: int = 3,
    trim_method: str = "witchi",
) -> Optional[SupermatrixResult]:
    """Align + trim each group block and concatenate into a gap-padded supermatrix.

    Args:
        taxa_group_seqs: ``{leaf_name -> {group -> raw_sequence}}`` with leaves
            already namespaced (``QUERY__<stem>`` / ``<genome_id>``) and filtered
            to ``>= min_markers`` groups upstream (Section 3).
        groups: Ordered group names (defines partition/column order).
        work_dir: Scratch dir for per-group alignments/pruning artefacts.
        out_dir: Destination dir for the supermatrix + partitions.
        basename: Output basename (``<basename>.supermatrix.faa`` / ``.partitions.txt``).
        threads: Worker threads for alignment/pruning.
        min_taxa_per_group: Groups with fewer aligned taxa are dropped.

    Returns:
        :class:`SupermatrixResult`, or ``None`` if no group survives.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    group_blocks: Dict[str, Dict[str, str]] = {}
    group_widths: Dict[str, int] = {}
    dropped_groups: List[str] = []

    for group in groups:
        block_in = {
            leaf: seqs[group]
            for leaf, seqs in taxa_group_seqs.items()
            if seqs.get(group)
        }
        if len(block_in) < min_taxa_per_group:
            if block_in:
                dropped_groups.append(group)
                logger.info(
                    "Dropping group %s from supermatrix: %d taxa (< %d)",
                    group,
                    len(block_in),
                    min_taxa_per_group,
                )
            continue

        records = [
            SimpleNamespace(id=leaf, seq=seq) for leaf, seq in sorted(block_in.items())
        ]
        aligned = align_sequences_pyfamsa(records, threads)
        aln_path = work_dir / f"{group}.aln.fasta"
        with open(aln_path, "w") as handle:
            for leaf, seq in aligned:
                handle.write(f">{leaf}\n{seq.replace('*', '-')}\n")

        block = _trim_group_alignment(aln_path, threads, set(block_in), method=trim_method)
        if not _uniform_nonempty(block):
            dropped_groups.append(group)
            logger.warning("Group %s produced an empty/ragged block; dropped", group)
            continue
        group_blocks[group] = block
        group_widths[group] = len(next(iter(block.values())))

    included = [group for group in groups if group in group_blocks]
    if not included:
        logger.warning("No groups survived supermatrix construction")
        return None

    leaves = sorted(taxa_group_seqs.keys())
    rows: Dict[str, str] = {}
    dropped_leaves: List[str] = []
    for leaf in leaves:
        parts: List[str] = []
        present = 0
        for group in included:
            block = group_blocks[group]
            if leaf in block:
                parts.append(block[leaf])
                present += 1
            else:
                parts.append("-" * group_widths[group])
        if present == 0:
            dropped_leaves.append(leaf)
            continue
        rows[leaf] = "".join(parts)

    if not rows:
        logger.warning("No taxa survived supermatrix concatenation")
        return None

    kept_leaves = [leaf for leaf in leaves if leaf in rows]

    supermatrix_faa = out_dir / f"{basename}.supermatrix.faa"
    with open(supermatrix_faa, "w") as handle:
        for leaf in kept_leaves:
            handle.write(f">{leaf}\n{rows[leaf]}\n")

    partitions = out_dir / f"{basename}.partitions.txt"
    _write_partitions(partitions, included, group_widths)

    if dropped_leaves:
        logger.info(
            "Excluded %d taxa absent from every surviving group block", len(dropped_leaves)
        )

    return SupermatrixResult(
        supermatrix_faa=supermatrix_faa,
        partitions=partitions,
        group_widths=group_widths,
        included_groups=included,
        leaves=kept_leaves,
        dropped_groups=dropped_groups,
        dropped_leaves=dropped_leaves,
    )


def _write_partitions(
    path: Path, included_groups: Sequence[str], group_widths: Dict[str, int]
) -> None:
    lines: List[str] = []
    start = 1
    for group in included_groups:
        end = start + group_widths[group] - 1
        lines.append(f"AA, {group} = {start}-{end}")
        start = end + 1
    path.write_text("\n".join(lines) + "\n")
