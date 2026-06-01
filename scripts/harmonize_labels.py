#!/usr/bin/env python3
"""Harmonize GVClass taxonomy labels into an active label bundle.

The script produces a candidate resource-label bundle without modifying the
input resources:

- ``labels.tsv``: active, FAA-backed labels only
- ``aliases.tsv``: old/source/inactive IDs mapped to active canonical IDs
- ``label_context.tsv``: host/context/original-taxonomy annotations
- ``inactive_labels.tsv``: labels with no FAA evidence and no safe alias target
- ``harmonization_summary.txt``: counts and validation summary

It is intentionally conservative: labels are not collapsed just because their
taxonomy strings match. Duplicate taxonomy labels are dereplicated only when
there is a unique FAA-backed canonical label or an explicit virophage PHAGE/PLV
to VP mapping.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
DOMAIN_PLACEHOLDERS = {
    "ARC": "ARC_unclassified",
    "BAC": "BAC_unclassified",
    "EUK": "EUK_unclassified",
    "MIRUS": "MIRUS_unclassified",
    "MITO": "MITO_unclassified",
    "NCLDV": "NCLDV_unclassified",
    "PHAGE": "PHAGE_unclassified",
    "PLASTID": "PLASTID_unclassified",
    "PLV": "PLV_unclassified",
    "UNCLASSIFIED": "UNCLASSIFIED_unclassified",
    "VP": "VP_unclassified",
}
KNOWN_SCAFFOLDS = {
    "MIRUS": ("MIRUS", "Mirusviricota"),
    "NCLDV": ("NCLDV", "Nucleocytoviricota"),
    "PLV": ("PLV", "Preplasmiviricota", "Polintoviricetes"),
    "VP": ("VP", "Preplasmiviricota", "Maveriviricetes"),
}
PLACEHOLDER_VALUES = {"", "na", "unknown", "unclassified"}
EUK_FORMAL_RANKS = ("phylum", "class", "order", "family", "genus", "species")


@dataclass(frozen=True)
class LabelRecord:
    label_id: str
    ranks: tuple[str, ...]
    line_no: int

    @property
    def taxonomy(self) -> str:
        return "|".join(self.ranks)


@dataclass
class Evidence:
    record_count: int = 0
    marker_files: set[str] | None = None

    def add(self, marker_file: str) -> None:
        self.record_count += 1
        if self.marker_files is None:
            self.marker_files = set()
        self.marker_files.add(marker_file)


def prefix_of(label_id: str) -> str:
    return label_id.split("__", 1)[0] if "__" in label_id else "NO_PREFIX"


def suffix_of(label_id: str) -> str:
    return label_id.split("__", 1)[1] if "__" in label_id else label_id


def accession_from_label(label_id: str) -> str:
    suffix = suffix_of(label_id)
    match = re.match(r"^(GC[AF])[-_](\d+)[-_\.](\d+)$", suffix)
    if match:
        return f"{match.group(1)}_{match.group(2)}.{match.group(3)}"
    return suffix


def accession_key(label_id: str) -> str:
    accession = accession_from_label(label_id)
    match = re.match(r"^(GC[AF])_(\d+)\.(\d+)$", accession)
    if match:
        return f"{match.group(2)}.{match.group(3)}"
    return accession


def preferred_accession_id(label_ids: list[str]) -> str:
    def score(label_id: str) -> tuple[int, int, str]:
        accession = accession_from_label(label_id)
        return (
            1 if accession.startswith("GCF_") else 0,
            1 if accession.startswith("GCA_") else 0,
            accession,
        )

    return sorted(label_ids, key=score, reverse=True)[0]


def placeholder_for(domain: str) -> str:
    return DOMAIN_PLACEHOLDERS.get(domain, f"{domain}_unclassified")


def read_labels(path: Path) -> dict[str, LabelRecord]:
    labels: dict[str, LabelRecord] = {}
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            try:
                label_id, taxonomy = line.split("\t", 1)
            except ValueError as exc:
                raise ValueError(f"Malformed labels.tsv line {line_no}: {line!r}") from exc
            ranks = tuple(taxonomy.split("|"))
            if len(ranks) != 7:
                raise ValueError(
                    f"Malformed taxonomy at line {line_no}: expected 7 ranks, got {len(ranks)}"
                )
            if label_id in labels:
                raise ValueError(f"Duplicate label id in labels.tsv: {label_id}")
            labels[label_id] = LabelRecord(label_id=label_id, ranks=ranks, line_no=line_no)
    return labels


def candidate_label_ids(protein_id: str) -> list[str]:
    candidates = [protein_id]

    stripped = re.sub(r"_\d+$", "", protein_id)
    if stripped != protein_id:
        candidates.append(stripped)

    if "|" in protein_id:
        base, payload = protein_id.split("|", 1)
        candidates.append(base)
        contig = payload.rsplit("_", 1)[0] if "_" in payload else payload
        if contig:
            candidates.append(f"{base}_{contig}")

    return list(dict.fromkeys(candidates))


def iter_faa_headers(faa_dir: Path) -> tuple[str, str]:
    for faa in sorted(faa_dir.glob("*.faa")):
        rel = faa.name
        with faa.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    yield line[1:].split()[0], rel


def collect_evidence(
    labels: dict[str, LabelRecord], faa_dir: Path
) -> tuple[dict[str, Evidence], list[str], Counter]:
    label_ids = set(labels)
    evidence: dict[str, Evidence] = defaultdict(Evidence)
    unresolved_headers: list[str] = []
    strategy_counts: Counter = Counter()

    for header_id, marker_file in iter_faa_headers(faa_dir):
        resolved = None
        for candidate in candidate_label_ids(header_id):
            if candidate in label_ids:
                resolved = candidate
                break
        if resolved is None:
            unresolved_headers.append(header_id)
        else:
            evidence[resolved].add(marker_file)
            strategy_counts[prefix_of(header_id)] += 1

    return evidence, unresolved_headers, strategy_counts


def read_taxonkit_euk_lineages(
    species_names: set[str], taxonkit_data_dir: Path | None
) -> dict[str, dict[str, str]]:
    if not species_names or taxonkit_data_dir is None:
        return {}
    if not (taxonkit_data_dir / "names.dmp").exists():
        return {}

    name2taxid = subprocess.run(
        [
            "taxonkit",
            "name2taxid",
            "--data-dir",
            str(taxonkit_data_dir),
            "--sci-name",
        ],
        input="\n".join(sorted(species_names)) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    taxid_rows = []
    for line in name2taxid.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1]:
            taxid_rows.append((fields[0], fields[1]))
    if not taxid_rows:
        return {}

    lineage_input = "\n".join(f"{name}\t{taxid}" for name, taxid in taxid_rows) + "\n"
    lineage = subprocess.run(
        [
            "taxonkit",
            "lineage",
            "--data-dir",
            str(taxonkit_data_dir),
            "-i",
            "2",
            "-r",
            "-R",
            "-n",
        ],
        input=lineage_input,
        text=True,
        capture_output=True,
        check=True,
    )

    candidates: dict[str, list[dict[str, str]]] = defaultdict(list)
    for line in lineage.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        query_name, taxid, lineage_names, sci_name, terminal_rank, lineage_ranks = fields[:6]
        names = lineage_names.split(";")
        ranks = lineage_ranks.split(";")
        if "Eukaryota" not in names:
            continue
        formal = {rank: "" for rank in EUK_FORMAL_RANKS}
        for name, rank in zip(names, ranks):
            if rank in formal and not formal[rank]:
                formal[rank] = name
        if not formal["species"] and terminal_rank == "species":
            formal["species"] = sci_name
        if not formal["species"]:
            continue
        formal["taxid"] = taxid
        formal["terminal_rank"] = terminal_rank
        formal["scientific_name"] = sci_name
        score = sum(1 for rank in EUK_FORMAL_RANKS if formal[rank])
        formal["_score"] = str(score)
        candidates[query_name].append(formal)

    selected: dict[str, dict[str, str]] = {}
    for query_name, rows in candidates.items():
        selected[query_name] = sorted(
            rows,
            key=lambda row: (
                int(row["_score"]),
                row["terminal_rank"] == "species",
                row["terminal_rank"] in {"strain", "isolate", "subspecies"},
            ),
            reverse=True,
        )[0]
    return selected


def read_taxonkit_lineages(
    species_names: set[str], taxonkit_data_dir: Path | None
) -> dict[str, dict[str, str]]:
    if not species_names or taxonkit_data_dir is None:
        return {}
    if not (taxonkit_data_dir / "names.dmp").exists():
        return {}

    name2taxid = subprocess.run(
        [
            "taxonkit",
            "name2taxid",
            "--data-dir",
            str(taxonkit_data_dir),
            "--sci-name",
        ],
        input="\n".join(sorted(species_names)) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    name_to_taxid: dict[str, str] = {}
    for line in name2taxid.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1]:
            name_to_taxid[fields[0]] = fields[1]
    if not name_to_taxid:
        return {}

    lineage = subprocess.run(
        [
            "taxonkit",
            "lineage",
            "--data-dir",
            str(taxonkit_data_dir),
            "-r",
            "-n",
        ],
        input="\n".join(name_to_taxid.values()) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    taxid_rows: dict[str, dict[str, str]] = {}
    for line in lineage.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        taxid, lineage_names, scientific_name, terminal_rank = fields[:4]
        taxid_rows[taxid] = {
            "taxid": taxid,
            "lineage": lineage_names,
            "scientific_name": scientific_name,
            "terminal_rank": terminal_rank,
            "domain": lineage_names.split(";", 1)[0],
        }
    return {
        name: taxid_rows[taxid]
        for name, taxid in name_to_taxid.items()
        if taxid in taxid_rows
    }


def normalize_unclassified(value: str, domain: str) -> tuple[str, str | None]:
    stripped = value.strip()
    low = stripped.lower()
    placeholder = placeholder_for(domain)
    if low in PLACEHOLDER_VALUES:
        return placeholder, None
    if low.startswith("unclassified "):
        return placeholder, stripped
    if stripped == domain:
        return placeholder, None
    if stripped.startswith(f"{domain}__"):
        return placeholder, stripped
    if stripped.startswith("MCP in "):
        return placeholder, stripped
    return stripped, None


def scaffold_known_ranks(domain: str, ranks: list[str]) -> None:
    if domain == "MIRUS" and (
        ranks[1] in {"MIRUS", "unclassified", "MIRUS_unclassified"}
        or ranks[1].startswith("MCP in ")
    ):
        ranks[1] = "Mirusviricota"
    elif domain == "NCLDV" and ranks[1] in {"NCLDV", "unclassified", "NCLDV_unclassified"}:
        ranks[1] = "Nucleocytoviricota"
    elif domain in {"PLV", "VP"}:
        scaffold = KNOWN_SCAFFOLDS[domain]
        if ranks[1] in {domain, "unclassified", f"{domain}_unclassified"} or ranks[1].startswith("MCP in "):
            ranks[1] = scaffold[1]
        if ranks[2] in {domain, "unclassified", f"{domain}_unclassified"} or ranks[2].startswith("MCP in "):
            ranks[2] = scaffold[2]


def euk_formal_taxonomy(
    ranks: tuple[str, ...],
    euk_lineages: dict[str, dict[str, str]],
) -> tuple[list[str], dict[str, str] | None]:
    lineage = euk_lineages.get(ranks[6])
    if not lineage:
        return list(ranks), None
    normalized = ["EUK"]
    for rank in EUK_FORMAL_RANKS:
        normalized.append(lineage.get(rank) or "EUK_unclassified")
    return normalized, lineage


def harmonize_taxonomy(
    label_id: str,
    ranks: tuple[str, ...],
    euk_lineages: dict[str, dict[str, str]],
) -> tuple[tuple[str, ...], list[tuple[str, str]]]:
    domain = prefix_of(label_id)
    contexts: list[tuple[str, str]] = []

    if domain == "EUK":
        normalized, lineage = euk_formal_taxonomy(ranks, euk_lineages)
        if lineage and tuple(normalized) != ranks:
            contexts.append(("original_taxonomy", "|".join(ranks)))
            contexts.append(("ncbi_taxid", lineage["taxid"]))
            contexts.append(("ncbi_terminal_rank", lineage["terminal_rank"]))
            ranks = tuple(normalized)

    if domain == "UNCLASSIFIED":
        ranks = ("UNCLASSIFIED", *ranks[1:])

    normalized = list(ranks)
    normalized[0] = domain
    pre_scaffold = normalized.copy()
    scaffold_known_ranks(domain, normalized)
    for index in range(1, 3):
        if pre_scaffold[index] != normalized[index] and pre_scaffold[index].startswith("MCP in "):
            contexts.append((RANKS[index], pre_scaffold[index]))

    # Preserve species names/IDs unless they are an exact bare domain placeholder.
    for index in range(1, 6):
        new_value, context = normalize_unclassified(normalized[index], domain)
        if context:
            contexts.append((RANKS[index], context))
        normalized[index] = new_value

    if normalized[6].startswith("MCP in "):
        contexts.append(("species", normalized[6]))
        normalized[6] = suffix_of(label_id)
    elif normalized[6] == domain or normalized[6].lower() in PLACEHOLDER_VALUES:
        normalized[6] = placeholder_for(domain)

    return tuple(normalized), contexts


def canonical_virophage_taxonomy(
    phage_id: str, labels: dict[str, LabelRecord]
) -> tuple[str, ...]:
    suffix = suffix_of(phage_id)
    plv = labels.get(f"PLV__{suffix}")
    species = suffix
    ranks = None
    if plv:
        ranks = list(plv.ranks)
        species = ranks[6]
    if ranks and ranks[1] == "Preplasmiviricota" and ranks[2] == "Maveriviricetes":
        return ("VP", *ranks[1:])
    phage = labels[phage_id]
    phage_species = phage.ranks[6]
    if phage_species.startswith("PHAGE__"):
        phage_species = phage_species.split("__", 1)[1]
    species = phage_species if phage_species else species
    return (
        "VP",
        "Preplasmiviricota",
        "Maveriviricetes",
        "VP_unclassified",
        "VP_unclassified",
        "VP_unclassified",
        species,
    )


def ncldv_taxonomy_from_lineage(lineage: str) -> tuple[str, ...] | None:
    names = lineage.split(";") if lineage else []
    if "Nucleocytoviricota" not in names:
        return None

    phylum = "Nucleocytoviricota"
    cls = "Megaviricetes" if "Megaviricetes" in names else "NCLDV_unclassified"
    order = "NCLDV_unclassified"
    for candidate in ("Imitervirales", "Pimascovirales"):
        if candidate in names:
            order = candidate
            break
    family = "NCLDV_unclassified"
    for candidate in ("Mimiviridae", "Pithoviridae"):
        if candidate in names:
            family = candidate
            break
    genus = "NCLDV_unclassified"
    if family in names:
        family_index = names.index(family)
        for name in names[family_index + 1 : -1]:
            if (
                name.startswith("unclassified ")
                or name.endswith("virinae")
                or " " in name
            ):
                continue
            genus = name
            break
    species = names[-1] if names else "NCLDV_unclassified"
    return ("NCLDV", phylum, cls, order, family, genus, species)


def plv_taxonomy_from_lineage(
    lineage: str, fallback_ranks: tuple[str, ...]
) -> tuple[str, ...] | None:
    taxonomy_text = "|".join(fallback_ranks)
    names = lineage.split(";") if lineage else []
    if "Polintoviricetes" not in names and "Polintoviricetes" not in taxonomy_text:
        return None

    phylum = "Preplasmiviricota"
    cls = "Polintoviricetes"
    order = "PLV_unclassified"
    for name in names:
        if name.endswith("virales"):
            order = name
            break
    if order == "PLV_unclassified" and len(fallback_ranks) > 3:
        order = fallback_ranks[3]

    family = "PLV_unclassified"
    for name in names:
        if name.endswith("viridae") and not name.startswith("unclassified "):
            family = name
            break
    if family == "PLV_unclassified" and len(fallback_ranks) > 4:
        family = fallback_ranks[4]

    genus = "PLV_unclassified"
    if family in names:
        family_index = names.index(family)
        for name in names[family_index + 1 : -1]:
            if name.startswith("unclassified ") or " " in name:
                continue
            genus = name
            break
    if genus == "PLV_unclassified" and len(fallback_ranks) > 5:
        fallback_genus = fallback_ranks[5]
        if not fallback_genus.startswith("unclassified "):
            genus = fallback_genus

    species = names[-1] if names else fallback_ranks[6]
    return ("PLV", phylum, cls, order, family, genus, species)


def build_cross_domain_viral_relabels(
    labels: dict[str, LabelRecord],
    evidence: dict[str, Evidence],
    taxonkit_lineages: dict[str, dict[str, str]],
) -> dict[str, tuple[str, tuple[str, ...], str]]:
    relabels: dict[str, tuple[str, tuple[str, ...], str]] = {}
    ncldv_groups: dict[str, list[str]] = defaultdict(list)

    for label_id, record in labels.items():
        if label_id not in evidence:
            continue
        domain = prefix_of(label_id)
        lineage = taxonkit_lineages.get(record.ranks[6], {}).get("lineage", "")
        if domain != "NCLDV" and "Nucleocytoviricota" in lineage:
            ncldv_groups[accession_key(label_id)].append(label_id)
        elif domain == "PHAGE":
            plv_taxonomy = plv_taxonomy_from_lineage(lineage, record.ranks)
            if plv_taxonomy is not None:
                canonical_id = f"PLV__{suffix_of(label_id)}"
                relabels[label_id] = (
                    canonical_id,
                    plv_taxonomy,
                    "phage_plv_relabel",
                )

    for label_ids in ncldv_groups.values():
        preferred = preferred_accession_id(label_ids)
        preferred_accession = accession_from_label(preferred)
        canonical_id = f"NCLDV__{preferred_accession}"
        preferred_record = labels[preferred]
        lineage = taxonkit_lineages.get(preferred_record.ranks[6], {}).get("lineage", "")
        ncldv_taxonomy = ncldv_taxonomy_from_lineage(lineage)
        if ncldv_taxonomy is None:
            continue
        for label_id in label_ids:
            reason = f"{prefix_of(label_id).lower()}_ncldv_relabel"
            relabels[label_id] = (canonical_id, ncldv_taxonomy, reason)

    return relabels


def build_excluded_phage_labels(
    labels: dict[str, LabelRecord],
    evidence: dict[str, Evidence],
    taxonkit_lineages: dict[str, dict[str, str]],
    relabels: dict[str, tuple[str, tuple[str, ...], str]],
) -> set[str]:
    excluded: set[str] = set()
    for label_id, record in labels.items():
        if label_id not in evidence or prefix_of(label_id) != "PHAGE":
            continue
        if label_id in relabels:
            continue
        taxonomy_text = record.taxonomy
        lineage = taxonkit_lineages.get(record.ranks[6], {}).get("lineage", "")
        joined = f"{taxonomy_text}|{lineage}"
        if "Preplasmiviricota" not in joined:
            continue
        if any(
            term in joined
            for term in (
                "Maveriviricetes",
                "Lavidaviridae",
                "virophage",
                "Polintoviricetes",
                "Adintoviridae",
                "Eupolintoviridae",
            )
        ):
            continue
        excluded.add(label_id)
    return excluded


def build_taxonomy_aliases(
    labels: dict[str, LabelRecord],
    evidence: dict[str, Evidence],
) -> dict[str, str]:
    taxonomy_groups: dict[str, list[str]] = defaultdict(list)
    for label_id, record in labels.items():
        taxonomy_groups[record.taxonomy].append(label_id)

    aliases: dict[str, str] = {}
    for ids in taxonomy_groups.values():
        active = [label_id for label_id in ids if label_id in evidence]
        inactive = [label_id for label_id in ids if label_id not in evidence]
        if len(active) == 1:
            for label_id in inactive:
                aliases[label_id] = active[0]
    return aliases


def harmonize_bundle(
    labels: dict[str, LabelRecord],
    evidence: dict[str, Evidence],
    euk_lineages: dict[str, dict[str, str]],
    taxonkit_lineages: dict[str, dict[str, str]] | None = None,
) -> tuple[
    dict[str, tuple[str, ...]],
    dict[str, tuple[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    set[str],
]:
    active_labels: dict[str, tuple[str, ...]] = {}
    aliases: dict[str, tuple[str, str]] = {}
    context_rows: list[dict[str, str]] = []
    inactive_rows: list[dict[str, str]] = []
    if taxonkit_lineages is None:
        taxonkit_lineages = {}

    taxonomy_aliases = build_taxonomy_aliases(labels, evidence)

    # Explicit PHAGE/PLV virophage bridge to canonical VP IDs.
    for label_id in labels:
        if not label_id.startswith("PHAGE__GCA-") or label_id not in evidence:
            continue
        suffix = suffix_of(label_id)
        plv_id = f"PLV__{suffix}"
        if plv_id not in labels:
            continue
        canonical_id = f"VP__{suffix}"
        active_labels[canonical_id] = canonical_virophage_taxonomy(label_id, labels)
        aliases[label_id] = (canonical_id, "phage_virophage_relabel")
        aliases[plv_id] = (canonical_id, "plv_virophage_alias")
        context_rows.append(
            {
                "label_id": canonical_id,
                "context_type": "original_phage_taxonomy",
                "context_value": labels[label_id].taxonomy,
                "source": label_id,
            }
        )
        context_rows.append(
            {
                "label_id": canonical_id,
                "context_type": "original_plv_taxonomy",
                "context_value": labels[plv_id].taxonomy,
                "source": plv_id,
            }
        )

    viral_relabels = build_cross_domain_viral_relabels(
        labels, evidence, taxonkit_lineages
    )
    excluded_labels = build_excluded_phage_labels(
        labels, evidence, taxonkit_lineages, viral_relabels
    )
    for label_id, (canonical_id, ranks, reason) in sorted(viral_relabels.items()):
        if label_id in aliases:
            continue
        active_labels[canonical_id] = ranks
        aliases[label_id] = (canonical_id, reason)
        context_rows.append(
            {
                "label_id": canonical_id,
                "context_type": "original_taxonomy",
                "context_value": labels[label_id].taxonomy,
                "source": label_id,
            }
        )

    for label_id, record in labels.items():
        if label_id in aliases:
            continue
        if label_id in excluded_labels:
            inactive_rows.append(
                {
                    "label_id": label_id,
                    "taxonomy": record.taxonomy,
                    "reason": "excluded_non_phage_preplasmiviricota_from_phage",
                }
            )
            continue
        if label_id not in evidence:
            target = taxonomy_aliases.get(label_id)
            if target and target not in aliases:
                aliases[label_id] = (target, "duplicate_taxonomy_without_faa")
            else:
                inactive_rows.append(
                    {
                        "label_id": label_id,
                        "taxonomy": record.taxonomy,
                        "reason": "no_faa_evidence",
                    }
                )
            continue

        if label_id.startswith("PHAGE__GCA-") and label_id in aliases:
            continue
        ranks, contexts = harmonize_taxonomy(label_id, record.ranks, euk_lineages)
        active_labels[label_id] = ranks
        for context_type, context_value in contexts:
            context_rows.append(
                {
                    "label_id": label_id,
                    "context_type": context_type,
                    "context_value": context_value,
                    "source": "harmonize_labels.py",
                }
            )

    # Alias source duplicates only after all canonical active IDs are known.
    for alias_id, target in list(taxonomy_aliases.items()):
        if alias_id in aliases or alias_id in active_labels:
            continue
        if target in active_labels:
            aliases[alias_id] = (target, "duplicate_taxonomy_without_faa")

    return active_labels, aliases, context_rows, inactive_rows, excluded_labels


def write_labels(path: Path, labels: dict[str, tuple[str, ...]]) -> None:
    with path.open("w") as handle:
        for label_id in sorted(labels):
            handle.write(f"{label_id}\t{'|'.join(labels[label_id])}\n")


def write_aliases(
    path: Path,
    aliases: dict[str, tuple[str, str]],
    labels: dict[str, LabelRecord],
    active_labels: dict[str, tuple[str, ...]],
) -> None:
    with path.open("w", newline="") as handle:
        fieldnames = [
            "alias_id",
            "canonical_id",
            "reason",
            "alias_taxonomy",
            "canonical_taxonomy",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for alias_id in sorted(aliases):
            canonical_id, reason = aliases[alias_id]
            writer.writerow(
                {
                    "alias_id": alias_id,
                    "canonical_id": canonical_id,
                    "reason": reason,
                    "alias_taxonomy": labels[alias_id].taxonomy if alias_id in labels else "",
                    "canonical_taxonomy": "|".join(active_labels[canonical_id]),
                }
            )


def write_context(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        fieldnames = ["label_id", "context_type", "context_value", "source"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["label_id"], row["context_type"])))


def write_inactive(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        fieldnames = ["label_id", "taxonomy", "reason"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["label_id"]))


def validate_output(
    active_labels: dict[str, tuple[str, ...]],
    aliases: dict[str, tuple[str, str]],
    labels: dict[str, LabelRecord],
    evidence: dict[str, Evidence],
    excluded_labels: set[str] | None = None,
) -> tuple[list[str], list[str]]:
    alias_targets = {alias_id: target for alias_id, (target, _reason) in aliases.items()}
    active_ids = set(active_labels)
    excluded_labels = excluded_labels or set()
    unresolved = []
    active_with_evidence: Counter = Counter()

    for source_id in evidence:
        if source_id in excluded_labels:
            continue
        target = alias_targets.get(source_id, source_id)
        if target not in active_ids:
            unresolved.append(source_id)
        else:
            active_with_evidence[target] += evidence[source_id].record_count

    no_evidence = sorted(active_ids - set(active_with_evidence))
    missing_alias_targets = [
        f"{alias_id}->{target}"
        for alias_id, target in alias_targets.items()
        if alias_id in labels and target not in active_ids
    ]
    return sorted(unresolved) + sorted(missing_alias_targets), no_evidence


def rewrite_header_token(
    header_id: str, alias_targets: dict[str, str]
) -> tuple[str, str | None, str | None]:
    for candidate in candidate_label_ids(header_id):
        canonical_id = alias_targets.get(candidate)
        if canonical_id is None:
            continue
        new_header_id = None
        if header_id.startswith(candidate):
            new_header_id = canonical_id + header_id[len(candidate) :]
        if "|" in header_id:
            _base, payload = header_id.split("|", 1)
            new_header_id = f"{canonical_id}|{payload}"
        if new_header_id is not None and new_header_id != header_id:
            return new_header_id, candidate, canonical_id
    return header_id, None, None


def rewrite_ambiguous_header_token(
    header_id: str, active_label_ids: set[str]
) -> tuple[str, str | None, str | None]:
    _ = active_label_ids
    return header_id, None, None


def excluded_header_source(
    header_id: str, excluded_label_ids: set[str]
) -> str | None:
    for candidate in candidate_label_ids(header_id):
        if candidate in excluded_label_ids:
            return candidate
    return None


def write_rewritten_faa(
    source_faa_dir: Path,
    output_faa_dir: Path,
    aliases: dict[str, tuple[str, str]],
    active_labels: dict[str, tuple[str, ...]],
    excluded_label_ids: set[str] | None = None,
) -> dict[str, int]:
    alias_targets = {
        alias_id: canonical_id
        for alias_id, (canonical_id, _reason) in aliases.items()
    }
    _ = active_labels
    excluded_label_ids = excluded_label_ids or set()
    output_faa_dir.mkdir(parents=True, exist_ok=True)

    stats = {
        "files_seen": 0,
        "files_symlinked": 0,
        "files_rewritten": 0,
        "headers_seen": 0,
        "headers_rewritten": 0,
        "headers_dropped": 0,
    }
    rewrite_rows: list[dict[str, str]] = []

    for source in sorted(source_faa_dir.glob("*.faa")):
        stats["files_seen"] += 1
        target = output_faa_dir / source.name
        temp_target = output_faa_dir / f".{source.name}.tmp"
        changed_lines: list[str] = []
        rewrites_by_header_id: dict[str, str] = {}
        dropped_by_header_id: dict[str, str] = {}
        if temp_target.exists() or temp_target.is_symlink():
            temp_target.unlink()

        with source.open() as in_handle:
            for line in in_handle:
                if not line.startswith(">"):
                    continue
                stats["headers_seen"] += 1
                body = line[1:].rstrip("\n")
                header_id = body.split()[0]
                excluded_source = excluded_header_source(header_id, excluded_label_ids)
                if excluded_source is not None:
                    stats["headers_dropped"] += 1
                    dropped_by_header_id[header_id] = excluded_source
                    changed_lines.append(
                        "\t".join(
                            [
                                source.name,
                                header_id,
                                "",
                                "dropped_excluded_label",
                                excluded_source,
                                "",
                            ]
                        )
                    )
                    continue
                new_header_id, alias_id, canonical_id = rewrite_header_token(
                    header_id, alias_targets
                )
                reason = "alias_rewrite"

                if alias_id is not None and canonical_id is not None:
                    stats["headers_rewritten"] += 1
                    rewrites_by_header_id[header_id] = new_header_id
                    changed_lines.append(
                        "\t".join(
                            [
                                source.name,
                                header_id,
                                new_header_id,
                                reason,
                                alias_id,
                                canonical_id,
                            ]
                        )
                    )

        if rewrites_by_header_id or dropped_by_header_id:
            with source.open() as in_handle, temp_target.open("w") as out_handle:
                write_record = True
                for line in in_handle:
                    if not line.startswith(">"):
                        if write_record:
                            out_handle.write(line)
                        continue
                    body = line[1:].rstrip("\n")
                    header_id, sep, rest = body.partition(" ")
                    if header_id in dropped_by_header_id:
                        write_record = False
                        continue
                    write_record = True
                    new_header_id = rewrites_by_header_id.get(header_id, header_id)
                    out_handle.write(f">{new_header_id}{sep}{rest}\n")
            if target.exists() or target.is_symlink():
                target.unlink()
            temp_target.replace(target)
            stats["files_rewritten"] += 1
            for row in changed_lines:
                (
                    faa,
                    old_header,
                    new_header,
                    reason,
                    source_id,
                    target_id,
                ) = row.split("\t")
                rewrite_rows.append(
                    {
                        "faa": faa,
                        "old_header": old_header,
                        "new_header": new_header,
                        "rewrite_reason": reason,
                        "source_id": source_id,
                        "target_id": target_id,
                    }
                )
        else:
            if target.exists() or target.is_symlink():
                target.unlink()
            os.symlink(source.resolve(), target)
            stats["files_symlinked"] += 1

    rewrite_map = output_faa_dir.parents[1] / "faa_rewrite_map.tsv"
    with rewrite_map.open("w", newline="") as handle:
        fieldnames = [
            "faa",
            "old_header",
            "new_header",
            "rewrite_reason",
            "source_id",
            "target_id",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rewrite_rows)
    return stats


def write_summary(
    path: Path,
    labels: dict[str, LabelRecord],
    evidence: dict[str, Evidence],
    active_labels: dict[str, tuple[str, ...]],
    aliases: dict[str, tuple[str, str]],
    inactive_rows: list[dict[str, str]],
    unresolved: list[str],
    no_evidence: list[str],
    faa_rewrite_stats: dict[str, int] | None = None,
) -> None:
    active_by_prefix = Counter(prefix_of(label_id) for label_id in active_labels)
    inactive_by_prefix = Counter(prefix_of(row["label_id"]) for row in inactive_rows)
    inactive_reasons = Counter(row["reason"] for row in inactive_rows)
    alias_reasons = Counter(reason for _target, reason in aliases.values())

    with path.open("w") as handle:
        handle.write("GVClass labels harmonization summary\n\n")
        handle.write(f"input_labels: {len(labels)}\n")
        handle.write(f"input_labels_with_faa: {len(evidence)}\n")
        handle.write(f"active_output_labels: {len(active_labels)}\n")
        handle.write(f"aliases: {len(aliases)}\n")
        handle.write(f"inactive_labels: {len(inactive_rows)}\n")
        handle.write(f"validation_unresolved_evidence_sources: {len(unresolved)}\n")
        handle.write(f"validation_active_labels_without_evidence: {len(no_evidence)}\n")
        handle.write("\nactive_output_labels_by_prefix:\n")
        for prefix, count in sorted(active_by_prefix.items()):
            handle.write(f"  {prefix}: {count}\n")
        handle.write("\ninactive_labels_by_prefix:\n")
        for prefix, count in sorted(inactive_by_prefix.items()):
            handle.write(f"  {prefix}: {count}\n")
        handle.write("\ninactive_reasons:\n")
        for reason, count in sorted(inactive_reasons.items()):
            handle.write(f"  {reason}: {count}\n")
        handle.write("\nalias_reasons:\n")
        for reason, count in sorted(alias_reasons.items()):
            handle.write(f"  {reason}: {count}\n")
        if unresolved:
            handle.write("\nexamples_unresolved:\n")
            for item in unresolved[:20]:
                handle.write(f"  {item}\n")
        if no_evidence:
            handle.write("\nexamples_active_without_evidence:\n")
            for item in no_evidence[:20]:
                handle.write(f"  {item}\n")
        if faa_rewrite_stats is not None:
            handle.write("\nfaa_rewrite:\n")
            for key, value in faa_rewrite_stats.items():
                handle.write(f"  {key}: {value}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=Path("resources/labels.tsv"))
    parser.add_argument("--faa-dir", type=Path, default=Path("resources/database/faa"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tasks/issue16_harmonized_resources"),
    )
    parser.add_argument(
        "--taxonkit-data-dir",
        type=Path,
        default=Path("/home/fschulz/dev/nelli-genomes-db/resources/taxonomy"),
    )
    parser.add_argument(
        "--no-euk-taxonkit",
        action="store_true",
        help="skip TaxonKit lineage normalization for EUK species",
    )
    parser.add_argument(
        "--rewrite-faa",
        action="store_true",
        help="write output-dir/database/faa with canonicalized headers",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    labels = read_labels(args.labels)
    evidence, unresolved_headers, _strategy_counts = collect_evidence(labels, args.faa_dir)
    if unresolved_headers:
        raise SystemExit(
            f"{len(unresolved_headers)} FAA headers do not resolve to input labels; "
            f"first example: {unresolved_headers[0]}"
        )

    euk_species = {
        record.ranks[6]
        for record in labels.values()
        if record.label_id.startswith("EUK__") and record.label_id in evidence
    }
    euk_lineages = {}
    taxonkit_lineages = {}
    if not args.no_euk_taxonkit:
        euk_lineages = read_taxonkit_euk_lineages(euk_species, args.taxonkit_data_dir)
        active_species = {
            record.ranks[6]
            for record in labels.values()
            if record.label_id in evidence
        }
        taxonkit_lineages = read_taxonkit_lineages(
            active_species, args.taxonkit_data_dir
        )

    active_labels, aliases, context_rows, inactive_rows, excluded_labels = harmonize_bundle(
        labels, evidence, euk_lineages, taxonkit_lineages
    )
    unresolved, no_evidence = validate_output(
        active_labels, aliases, labels, evidence, excluded_labels
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_labels(args.output_dir / "labels.tsv", active_labels)
    write_aliases(args.output_dir / "aliases.tsv", aliases, labels, active_labels)
    write_context(args.output_dir / "label_context.tsv", context_rows)
    write_inactive(args.output_dir / "inactive_labels.tsv", inactive_rows)
    faa_rewrite_stats = None
    if args.rewrite_faa:
        faa_rewrite_stats = write_rewritten_faa(
            args.faa_dir,
            args.output_dir / "database" / "faa",
            aliases,
            active_labels,
            excluded_labels,
        )

    write_summary(
        args.output_dir / "harmonization_summary.txt",
        labels,
        evidence,
        active_labels,
        aliases,
        inactive_rows,
        unresolved,
        no_evidence,
        faa_rewrite_stats,
    )

    if unresolved or no_evidence:
        raise SystemExit(
            "harmonized bundle failed validation; see harmonization_summary.txt"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
