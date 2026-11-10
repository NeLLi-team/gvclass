"""Reference label index shared by tree analysis, taxonomy votes and contamination.

One loader reads ``labels.tsv`` once per file, one resolver maps a protein or leaf
name to its genome ID, and one rule turns the nearest tree neighbours of a protein
into a lineage.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType

from src.core.summarize.consensus import TAX_LEVELS, Lineage

__all__ = [
    "LabelIndex",
    "Placement",
    "common_lineage",
    "domain_prefix",
    "nearest_placements",
]

logger = logging.getLogger(__name__)

_GENE_SUFFIX = re.compile(r"_\d+$")
_LABEL_FIELDS = 2
_CACHED_INDEXES = 8

# Database v2.1.0 holds one genome twice with different labels. Each key is a
# duplicate record that takes the lineage of its value, the canonical record.
# Temporary: delete this constant and _with_canonical_duplicates when a database
# release corrects the record.
DUPLICATE_LABEL_RECORDS = {"NCLDV__GCF_905367645.1": "NCLDV__PkV-RF01"}


@dataclass(frozen=True, slots=True)
class Placement:
    """One labelled reference neighbour of a query protein in one marker tree."""

    distance: float
    marker: str
    neighbor: str
    genome_id: str
    lineage: Lineage


@dataclass(frozen=True, slots=True, eq=False)
class LabelIndex:
    """Immutable seven-rank lineage of every reference genome ID.

    Attributes:
        lineages: Read-only map from genome ID to its seven ranks, from domain to
            species. The domain is the ID prefix before ``__``. Ranks the label
            does not give are empty.
    """

    lineages: Mapping[str, Lineage]

    @classmethod
    def from_taxonomy(cls, taxonomy: Mapping[str, Sequence[str]]) -> LabelIndex:
        """Build an index from ``{genome ID: ranks}``, padding short rank lists."""
        lineages = {
            genome_id: _lineage(genome_id, ranks)
            for genome_id, ranks in taxonomy.items()
        }
        return cls(MappingProxyType(lineages))

    @classmethod
    def load(cls, labels_file: Path) -> LabelIndex:
        """Return the index of ``labels_file``, read once while the file is unchanged.

        Raises:
            OSError: The file cannot be read.
            ValueError: A record is not ``genome ID<TAB>taxonomy``, or the file
                holds no record.
        """
        try:
            stat = labels_file.stat()
        except OSError as exc:
            raise OSError(
                f"Failed to read required labels {labels_file}: {exc}"
            ) from exc
        return _load_index((labels_file.resolve(), stat.st_mtime_ns, stat.st_size))

    def resolve_genome_id(self, protein_id: str) -> str:
        """Map a protein or tree-leaf name to its genome ID.

        Precedence: the name itself when it is a label key (genome IDs can end in
        ``_<digits>``); the name without a trailing ``_<digits>`` gene number; the
        part before the first ``|``; else the name without the gene number, which
        has no label.
        """
        if protein_id in self.lineages:
            return protein_id
        stripped = _GENE_SUFFIX.sub("", protein_id)
        if stripped in self.lineages:
            return stripped
        base_id = protein_id.split("|", 1)[0]
        if base_id != protein_id and base_id in self.lineages:
            return base_id
        return stripped

    def placements(
        self, marker: str, neighbors: Mapping[str, float]
    ) -> list[Placement]:
        """Return the neighbours of one protein in one tree that have a label."""
        placements: list[Placement] = []
        for neighbor, distance in neighbors.items():
            genome_id = self.resolve_genome_id(neighbor)
            lineage = self.lineages.get(genome_id)
            if lineage is not None:
                placements.append(
                    Placement(distance, marker, neighbor, genome_id, lineage)
                )
        return placements

    def nearest_lineage(
        self, marker: str, neighbors: Mapping[str, float]
    ) -> tuple[Lineage, float] | None:
        """Return the shared lineage and distance of a protein's nearest neighbours.

        The neighbours are those of one protein in one tree. The lineage is the one
        the labelled neighbours at the minimal distance have in common, so a tie
        across orders has no order and a tie across domains has no domain. None
        when no neighbour has a label.
        """
        nearest = nearest_placements(self.placements(marker, neighbors))
        if not nearest:
            return None
        return common_lineage([placement.lineage for placement in nearest]), (
            nearest[0].distance
        )


def domain_prefix(genome_id: str) -> str:
    """Return the ``DOMAIN`` of a ``DOMAIN__name`` genome ID, or ``""``."""
    return genome_id.split("__", 1)[0] if "__" in genome_id else ""


def nearest_placements(placements: Sequence[Placement]) -> list[Placement]:
    """Return the placements at the minimal distance, sorted by marker and neighbour.

    Distances compare for exact equality: tied neighbours sit at the same
    patristic distance in the tree.
    """
    if not placements:
        return []
    minimal = min(placement.distance for placement in placements)
    tied = [placement for placement in placements if placement.distance == minimal]
    return sorted(tied, key=lambda placement: (placement.marker, placement.neighbor))


def common_lineage(lineages: Sequence[Lineage]) -> Lineage:
    """Return the longest common prefix of the lineages, padded with empty ranks."""
    shared: list[str] = []
    for ranks in zip(*lineages, strict=True):
        if len(set(ranks)) > 1:
            break
        shared.append(ranks[0])
    padding = [""] * (len(TAX_LEVELS) - len(shared))
    return (*shared, *padding)


def _lineage(genome_id: str, ranks: Sequence[str]) -> Lineage:
    """Return seven stripped ranks, with the domain taken from the genome ID."""
    kept = [rank.strip() for rank in ranks[: len(TAX_LEVELS)]]
    kept += [""] * (len(TAX_LEVELS) - len(kept))
    kept[0] = domain_prefix(genome_id) or kept[0]
    return tuple(kept)


@lru_cache(maxsize=_CACHED_INDEXES)
def _load_index(file_key: tuple[Path, int, int]) -> LabelIndex:
    """Read one labels file; the modification time and size in the key expire it."""
    labels_file, _mtime_ns, _size = file_key
    index = LabelIndex.from_taxonomy(
        _with_canonical_duplicates(_read_taxonomy(labels_file))
    )
    logger.info("Loaded %d labels from %s", len(index.lineages), labels_file)
    return index


def _with_canonical_duplicates(
    taxonomy: Mapping[str, Sequence[str]],
) -> dict[str, Sequence[str]]:
    """Return ``taxonomy`` with each duplicate record holding its canonical ranks.

    A duplicate keeps its own ranks when its canonical record is absent, so a
    corrected labels file needs no code change.
    """
    resolved = dict(taxonomy)
    for duplicate, canonical in DUPLICATE_LABEL_RECORDS.items():
        if duplicate in resolved and canonical in resolved:
            resolved[duplicate] = resolved[canonical]
    return resolved


def _read_taxonomy(labels_file: Path) -> dict[str, list[str]]:
    """Read ``genome ID<TAB>domain|phylum|...`` records; reject any other record."""
    taxonomy: dict[str, list[str]] = {}
    try:
        with labels_file.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                record = line.rstrip("\r\n")
                if not record or record.startswith("#"):
                    continue
                fields = record.split("\t")
                if len(fields) != _LABEL_FIELDS or not all(fields):
                    raise ValueError(
                        f"Malformed required label record in {labels_file} at line "
                        f"{line_number}: expected nonempty genome ID and taxonomy "
                        "separated by one tab"
                    )
                taxonomy[fields[0]] = fields[1].split("|")
    except OSError as exc:
        raise OSError(f"Failed to read required labels {labels_file}: {exc}") from exc
    if not taxonomy:
        raise ValueError(f"Required labels file {labels_file} has no records")
    return taxonomy
