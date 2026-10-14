#!/usr/bin/env python
"""Prepare, run, and summarize nucleotide-space quality benchmarks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

from Bio import SeqIO

from scripts._contamination_features import (
    COMPLETENESS_MODE,
    compute_query_features,
    contamination_feature_vector,
)
from src.core.contamination_scoring import CONTAMINATION_MODEL_FEATURES
from src.utils.common import sha256_file

RETENTION_LEVELS = (0.10, 0.25, 0.50, 0.75, 1.00)
CELLULAR_LEVELS = (0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
VIRAL_LEVELS = (0.10, 0.25, 0.75, 0.90)
CHIMERA_LEVEL = 0.25
MANIFEST_VERSION = 1
AUDIT_FILENAME = "quality_feature_audit.mode_v2.json"
ORIGINAL_AUDIT_FILENAME = "quality_feature_audit.json"


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _validate_identifier(value: str, label: str) -> None:
    if (
        not value
        or value in {".", ".."}
        or any(character.isspace() or character in "/\\" for character in value)
    ):
        raise ValueError(f"Unsafe {label}: {value!r}")


def _resolve_path(value: str, manifest_path: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    if path.exists():
        return path.resolve()
    return (manifest_path.parent / path).resolve()


def _load_json_rows(path: Path) -> List[Dict[str, Any]]:
    data = json.loads(path.read_text())
    if not isinstance(data, list) or not all(isinstance(row, dict) for row in data):
        raise ValueError(f"Expected a JSON list of objects in {path}")
    return data


def _load_donor_rows(path: Path) -> List[Dict[str, Any]]:
    data = json.loads(path.read_text())
    rows = data.get("donors", []) if isinstance(data, dict) else []
    if not rows or not all(isinstance(row, dict) for row in rows):
        raise ValueError(f"Expected a donors list in {path}")
    if not data.get("donor_summary", {}).get("all_checks_pass"):
        raise ValueError(f"Donor QC did not pass in {path}")
    for row in rows:
        if not all(row.get("checks", {}).values()):
            raise ValueError(f"Donor QC failed for {row.get('accession', 'unknown')}")
    return rows


def _load_genome(row: Dict[str, Any], manifest_path: Path) -> Dict[str, Any]:
    identifier = str(
        row.get("source_genome_id") or row.get("donor_id") or row.get("accession") or ""
    )
    if not identifier:
        raise ValueError(f"Genome row has no source or donor ID: {row}")
    _validate_identifier(identifier, "genome ID")
    row_paths = row.get("paths", {})
    fna_path = _resolve_path(
        str(row.get("fna_path") or row_paths["fasta"]), manifest_path
    )
    genbank_path = _resolve_path(
        str(row.get("genbank_path") or row_paths["genbank"]), manifest_path
    )
    row_hashes = row.get("sha256", {})
    expected_file_hashes = {
        fna_path: str(row.get("fna_sha256") or row_hashes.get("fasta_file") or ""),
        genbank_path: str(
            row.get("genbank_sha256") or row_hashes.get("genbank_file") or ""
        ),
    }
    for path, expected in expected_file_hashes.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        if expected and sha256_file(path) != expected:
            raise ValueError(f"Checksum mismatch for {path}")

    fna_records = list(SeqIO.parse(str(fna_path), "fasta"))
    genbank_records = list(SeqIO.parse(str(genbank_path), "genbank"))
    if len(fna_records) != 1 or len(genbank_records) != 1:
        raise ValueError(
            f"Benchmark genomes must contain one sequence record: {identifier}"
        )
    sequence = str(fna_records[0].seq).upper()
    if sequence != str(genbank_records[0].seq).upper():
        raise ValueError(f"FASTA and GenBank sequences differ for {identifier}")
    expected_sequence = str(
        row.get("sequence_sha256") or row_hashes.get("normalized_sequence") or ""
    )
    if expected_sequence and _sha256_bytes(sequence.encode()) != expected_sequence:
        raise ValueError(f"Sequence checksum mismatch for {identifier}")

    cds_intervals: List[Tuple[int, int]] = []
    cds_sequences: List[Dict[str, Any]] = []
    record = genbank_records[0]
    for feature_index, feature in enumerate(record.features):
        if feature.type != "CDS":
            continue
        parts = getattr(feature.location, "parts", [feature.location])
        for part in parts:
            start, end = int(part.start), int(part.end)
            if 0 <= start < end <= len(sequence):
                cds_intervals.append((start, end))
        extracted = str(feature.extract(record.seq)).upper()
        if extracted:
            cds_sequences.append(
                {
                    "feature_index": feature_index,
                    "start": min(int(part.start) for part in parts),
                    "end": max(int(part.end) for part in parts),
                    "sequence": extracted,
                }
            )
    return {
        "id": identifier,
        "record_id": fna_records[0].id,
        "sequence": sequence,
        "cds_intervals": sorted(cds_intervals),
        "cds_sequences": cds_sequences,
        "fna_path": str(fna_path),
        "genbank_path": str(genbank_path),
        "fna_sha256": sha256_file(fna_path),
        "genbank_sha256": sha256_file(genbank_path),
        "sequence_sha256": _sha256_bytes(sequence.encode()),
        "row": row,
    }


def _select_span(genome: Dict[str, Any], target_bp: int, key: str) -> Tuple[int, int]:
    sequence_length = len(genome["sequence"])
    target_bp = min(max(1, int(target_bp)), sequence_length)
    if target_bp == sequence_length:
        return 0, sequence_length
    max_start = sequence_length - target_bp
    start = int(hashlib.sha256(key.encode()).hexdigest()[:16], 16) % (max_start + 1)
    return start, start + target_bp


def _cut_through_cds(genome: Dict[str, Any], position: int) -> bool:
    return any(start < position < end for start, end in genome["cds_intervals"])


def _span_contigs(
    genome: Dict[str, Any],
    start: int,
    end: int,
    fragment_bp: int,
    role: str,
    id_prefix: str,
) -> List[Dict[str, Any]]:
    cuts = list(range(start, end, fragment_bp))
    if cuts[-1] != end:
        cuts.append(end)
    contigs = []
    for index, (left, right) in enumerate(zip(cuts, cuts[1:]), start=1):
        sequence = genome["sequence"][left:right]
        contigs.append(
            {
                "id": f"{id_prefix}_{index:04d}",
                "sequence": sequence,
                "parts": [
                    {
                        "role": role,
                        "parent_id": genome["id"],
                        "parent_record_id": genome["record_id"],
                        "start": left,
                        "end": right,
                        "length_bp": right - left,
                        "sequence_sha256": _sha256_bytes(sequence.encode()),
                        "left_cut_through_cds": _cut_through_cds(genome, left),
                        "right_cut_through_cds": _cut_through_cds(genome, right),
                    }
                ],
            }
        )
    return contigs


def _write_fasta(contigs: Sequence[Dict[str, Any]], path: Path) -> None:
    lines: List[str] = []
    for contig in contigs:
        lines.append(f">{contig['id']}")
        sequence = str(contig["sequence"])
        lines.extend(
            sequence[index : index + 80] for index in range(0, len(sequence), 80)
        )
    data = ("\n".join(lines) + "\n").encode()
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(data)
    temporary.replace(path)


def _taxonomy_overlap(left: Dict[str, Any], right: Dict[str, Any]) -> int:
    left_taxonomy = set(str(item) for item in left.get("taxonomy", []))
    return len(
        left_taxonomy.intersection(str(item) for item in right.get("taxonomy", []))
    )


def _viral_donors(
    source: Dict[str, Any], rows: Sequence[Dict[str, Any]]
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    candidates = [
        row
        for row in rows
        if row["split"] == source["split"]
        and row["source_genome_id"] != source["source_genome_id"]
    ]
    if len(candidates) < 2:
        raise ValueError(
            f"Need at least two split-local viral donors for {source['split']}"
        )
    same_panel = [row for row in candidates if row.get("panel") == source.get("panel")]
    other_panels = [
        row for row in candidates if row.get("panel") != source.get("panel")
    ]
    if not same_panel or not other_panels:
        raise ValueError(
            f"Need same-panel and other-panel donors for {source['source_genome_id']}"
        )
    related = sorted(
        same_panel,
        key=lambda row: (
            -_taxonomy_overlap(source, row),
            str(row["source_genome_id"]),
        ),
    )[0]
    divergent = sorted(
        other_panels,
        key=lambda row: (
            _taxonomy_overlap(source, row),
            str(row.get("panel", "")),
            str(row["source_genome_id"]),
        ),
    )[0]
    return related, divergent


def _sample_row(
    source: Dict[str, Any],
    sample_id: str,
    scenario: str,
    contigs: Sequence[Dict[str, Any]],
    source_total_bp: int,
    nominal_contamination: float,
    matched_control_sample_id: str,
    interpretation: str,
) -> Dict[str, Any]:
    output_bp = sum(len(str(contig["sequence"])) for contig in contigs)
    role_bp: Dict[str, int] = {}
    for contig in contigs:
        for part in contig["parts"]:
            role = str(part["role"])
            role_bp[role] = role_bp.get(role, 0) + int(part["length_bp"])
    base_retained_bp = role_bp.get("intended_target", 0)
    contamination_truth_bp = role_bp.get("contaminant", 0)
    foreign_sequence_bp = output_bp - base_retained_bp
    parent_ids = sorted(
        {str(part["parent_id"]) for contig in contigs for part in contig["parts"]}
    )
    truncating_boundaries = sorted(
        {
            (str(part["parent_id"]), int(position), str(part["role"]))
            for contig in contigs
            for part in contig["parts"]
            for position, crosses in (
                (part["start"], part.get("left_cut_through_cds", False)),
                (part["end"], part.get("right_cut_through_cds", False)),
            )
            if crosses
        }
    )
    return {
        "sample_id": sample_id,
        "source_genome_id": source["source_genome_id"],
        "source_split": source["split"],
        "source_family": source.get("panel", ""),
        "scenario": scenario,
        "truth_interpretation": interpretation,
        "source_total_bp": source_total_bp,
        "base_retained_bp": base_retained_bp,
        "foreign_sequence_bp": foreign_sequence_bp,
        "contamination_truth_bp": contamination_truth_bp,
        "output_bp": output_bp,
        "nominal_contamination_pct": round(nominal_contamination * 100.0, 6),
        "true_source_retention_pct": round(
            (base_retained_bp / source_total_bp) * 100.0, 10
        ),
        "true_contamination_pct": round(
            (contamination_truth_bp / output_bp) * 100.0, 10
        ),
        "matched_control_sample_id": matched_control_sample_id,
        "all_parent_ids": parent_ids,
        "cds_truncating_boundaries": [
            {"parent_id": parent, "position": position, "role": role}
            for parent, position, role in truncating_boundaries
        ],
        "contigs": [
            {
                "contig_id": contig["id"],
                "length_bp": len(str(contig["sequence"])),
                "sequence_sha256": _sha256_bytes(str(contig["sequence"]).encode()),
                "parts": contig["parts"],
            }
            for contig in contigs
        ],
    }


def _add_sample(
    samples: List[Dict[str, Any]],
    source: Dict[str, Any],
    contigs: Sequence[Dict[str, Any]],
    inputs_dir: Path,
    **fields: Any,
) -> Dict[str, Any]:
    row = _sample_row(source=source, contigs=contigs, **fields)
    _validate_identifier(str(row["sample_id"]), "sample ID")
    path = inputs_dir / f"{row['sample_id']}.fna"
    _write_fasta(contigs, path)
    row["input_path"] = str(path)
    row["input_sha256"] = sha256_file(path)
    samples.append(row)
    return row


def _base_contigs(
    source_genome: Dict[str, Any], target_bp: int, fragment_bp: int, key: str
) -> Tuple[List[Dict[str, Any]], int, Tuple[int, int]]:
    start, end = _select_span(source_genome, target_bp, key)
    contigs = _span_contigs(
        source_genome,
        start,
        end,
        fragment_bp,
        "intended_target",
        f"{source_genome['id']}__target",
    )
    return contigs, end - start, (start, end)


def _mixed_sample_parts(
    source_genome: Dict[str, Any],
    donor_genome: Dict[str, Any],
    scenario: str,
    fraction: float,
    fragment_bp: int,
) -> Tuple[List[Dict[str, Any]], int, int]:
    donor_available = len(donor_genome["sequence"])
    source_available = len(source_genome["sequence"])
    base_target = min(
        source_available,
        max(1, math.floor(donor_available * (1.0 - fraction) / fraction)),
    )
    base_contigs, base_bp, _ = _base_contigs(
        source_genome, base_target, fragment_bp, f"{scenario}:base"
    )
    donor_target = min(
        donor_available,
        max(1, round(base_bp * fraction / (1.0 - fraction))),
    )
    donor_start, donor_end = _select_span(
        donor_genome, donor_target, f"{scenario}:donor"
    )
    donor_contigs = _span_contigs(
        donor_genome,
        donor_start,
        donor_end,
        fragment_bp,
        "contaminant",
        f"{donor_genome['id']}__donor",
    )
    return [*base_contigs, *donor_contigs], base_bp, donor_end - donor_start


def _copy_contigs(contigs: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return json.loads(json.dumps(contigs))


def _target_only(contigs: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [
        contig
        for contig in _copy_contigs(contigs)
        if all(part["role"] == "intended_target" for part in contig["parts"])
    ]


def _add_mixture_series(
    samples: List[Dict[str, Any]],
    source: Dict[str, Any],
    source_genome: Dict[str, Any],
    donor_genome: Dict[str, Any],
    inputs_dir: Path,
    scenario_prefix: str,
    fractions: Sequence[float],
    fragment_bp: int,
    clean_full_id: str,
    interpretation: str,
) -> None:
    source_total_bp = len(source_genome["sequence"])
    for fraction in fractions:
        scenario = f"{scenario_prefix}_{round(fraction * 100):03d}"
        sample_id = f"{source_genome['id']}__{scenario}"
        contigs, base_bp, _ = _mixed_sample_parts(
            source_genome, donor_genome, sample_id, fraction, fragment_bp
        )
        matched_id = clean_full_id
        if base_bp != source_total_bp:
            matched_id = f"{sample_id}__matched_clean"
            _add_sample(
                samples,
                source,
                _target_only(contigs),
                inputs_dir,
                sample_id=matched_id,
                scenario=f"clean_matched_{scenario}",
                source_total_bp=source_total_bp,
                nominal_contamination=0.0,
                matched_control_sample_id=matched_id,
                interpretation="matched_source_retention_control",
            )
        _add_sample(
            samples,
            source,
            contigs,
            inputs_dir,
            sample_id=sample_id,
            scenario=scenario,
            source_total_bp=source_total_bp,
            nominal_contamination=fraction,
            matched_control_sample_id=matched_id,
            interpretation=interpretation,
        )


def _validate_sample(row: Dict[str, Any], expected_split: str) -> None:
    _validate_identifier(str(row.get("sample_id", "")), "sample ID")
    if row["source_split"] != expected_split:
        raise ValueError(f"Sample {row['sample_id']} has the wrong split")
    input_path = Path(row["input_path"])
    if not input_path.is_file() or sha256_file(input_path) != row["input_sha256"]:
        raise ValueError(f"Invalid prepared input for {row['sample_id']}")
    contig_bp = sum(int(contig["length_bp"]) for contig in row["contigs"])
    if contig_bp != int(row["output_bp"]):
        raise ValueError(f"Contig lengths do not sum for {row['sample_id']}")
    part_bp = sum(
        int(part["length_bp"]) for contig in row["contigs"] for part in contig["parts"]
    )
    if part_bp != contig_bp:
        raise ValueError(f"Parent spans do not sum for {row['sample_id']}")
    seen: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for contig in row["contigs"]:
        for part in contig["parts"]:
            if part["role"] == "hypothetical_hgt_intended":
                continue
            key = (str(part["role"]), str(part["parent_id"]))
            interval = (int(part["start"]), int(part["end"]))
            prior = seen.setdefault(key, [])
            if any(interval[0] < end and start < interval[1] for start, end in prior):
                raise ValueError(f"Repeated parent span in {row['sample_id']}: {key}")
            prior.append(interval)
    expected_contamination = int(row["contamination_truth_bp"]) / int(row["output_bp"])
    if not math.isclose(
        float(row["true_contamination_pct"]),
        expected_contamination * 100.0,
        abs_tol=1e-8,
    ):
        raise ValueError(f"Contamination truth ratio differs for {row['sample_id']}")


def prepare_dataset(
    source_manifest: Path,
    donor_manifest: Path,
    work_dir: Path,
    split: str,
    fragment_bp: int = 20_000,
    force: bool = False,
    allow_final_test: bool = False,
) -> Path:
    if split == "final_test" and not allow_final_test:
        raise ValueError("Final-test preparation requires --allow-final-test")
    if fragment_bp <= 0:
        raise ValueError("Fragment length must be positive")
    source_manifest = source_manifest.resolve()
    donor_manifest = donor_manifest.resolve()
    work_dir = work_dir.resolve()
    profile = {
        "manifest_version": MANIFEST_VERSION,
        "source_manifest": str(source_manifest),
        "source_manifest_sha256": sha256_file(source_manifest),
        "donor_manifest": str(donor_manifest),
        "donor_manifest_sha256": sha256_file(donor_manifest),
        "split": split,
        "fragment_bp": fragment_bp,
        "retention_levels": list(RETENTION_LEVELS),
        "cellular_levels": list(CELLULAR_LEVELS),
        "viral_levels": list(VIRAL_LEVELS),
    }
    profile["profile_sha256"] = _sha256_bytes(_canonical_json(profile).encode())
    profile_path = work_dir / "prepare_profile.json"
    manifest_path = work_dir / "benchmark_manifest.json"
    if profile_path.exists():
        existing = json.loads(profile_path.read_text())
        if existing != profile:
            raise RuntimeError(
                "Existing benchmark uses different inputs or split; use a new work directory"
            )
        if manifest_path.exists() and not force:
            manifest = json.loads(manifest_path.read_text())
            for row in manifest["samples"]:
                _validate_sample(row, split)
            return manifest_path

    source_rows = _load_json_rows(source_manifest)
    for row in source_rows:
        _validate_identifier(str(row.get("source_genome_id", "")), "source genome ID")
    selected_sources = sorted(
        (row for row in source_rows if row.get("split") == split),
        key=lambda row: str(row["source_genome_id"]),
    )
    if not selected_sources:
        raise ValueError(f"No sources found for split {split}")
    source_ids_by_split = {
        name: {
            str(row["source_genome_id"])
            for row in source_rows
            if row.get("split") == name
        }
        for name in ("development", "final_test")
    }
    if source_ids_by_split["development"].intersection(
        source_ids_by_split["final_test"]
    ):
        raise ValueError("Source genome IDs overlap development and final-test splits")

    donor_rows = _load_donor_rows(donor_manifest)
    for row in donor_rows:
        _validate_identifier(str(row.get("accession", "")), "cellular donor ID")
    donor_ids_by_split = {
        name: {str(row["accession"]) for row in donor_rows if row.get("split") == name}
        for name in ("development", "final_test")
    }
    if donor_ids_by_split["development"].intersection(donor_ids_by_split["final_test"]):
        raise ValueError("Cellular donor IDs overlap benchmark splits")
    if (
        set()
        .union(*source_ids_by_split.values())
        .intersection(set().union(*donor_ids_by_split.values()))
    ):
        raise ValueError("A cellular donor is also a viral source genome")
    cellular_rows = sorted(
        (
            row
            for row in donor_rows
            if row.get("split") == split
            and str(row.get("role", "")).endswith("_contamination_donor")
        ),
        key=lambda row: str(row["accession"]),
    )
    roles = {str(row["role"]) for row in cellular_rows}
    required_roles = {
        "bacterial_contamination_donor",
        "eukaryotic_contamination_donor",
    }
    if roles != required_roles:
        raise ValueError(f"Split {split} needs one bacterial and one eukaryotic donor")
    cellular_genomes = [_load_genome(row, donor_manifest) for row in cellular_rows]
    hgt_genome = next(
        genome
        for genome in cellular_genomes
        if genome["row"]["role"] == "eukaryotic_contamination_donor"
    )
    if len(hgt_genome["cds_sequences"]) < 5:
        raise ValueError(
            "Eukaryotic donor must have at least five annotated CDS features"
        )

    work_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = work_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    samples: List[Dict[str, Any]] = []
    genome_cache: Dict[str, Dict[str, Any]] = {
        genome["id"]: genome for genome in cellular_genomes
    }

    for source in selected_sources:
        source_id = str(source["source_genome_id"])
        source_genome = _load_genome(source, source_manifest)
        genome_cache[source_id] = source_genome
        source_total_bp = len(source_genome["sequence"])
        clean_full_id = f"{source_id}__clean_retention_100"

        for retention in RETENTION_LEVELS:
            sample_id = f"{source_id}__clean_retention_{round(retention * 100):03d}"
            contigs, _, _ = _base_contigs(
                source_genome,
                round(source_total_bp * retention),
                fragment_bp,
                sample_id,
            )
            _add_sample(
                samples,
                source,
                contigs,
                inputs_dir,
                sample_id=sample_id,
                scenario=f"clean_retention_{round(retention * 100):03d}",
                source_total_bp=source_total_bp,
                nominal_contamination=0.0,
                matched_control_sample_id=sample_id,
                interpretation="source_retention_truth_not_absolute_genome_completeness",
            )

        full_contigs, _, _ = _base_contigs(
            source_genome, source_total_bp, fragment_bp, clean_full_id
        )
        for gene_count in (1, 5):
            sample_id = f"{source_id}__hgt_like_{gene_count}cds"
            contigs = _copy_contigs(full_contigs)
            selected_cds = hgt_genome["cds_sequences"][:gene_count]
            insert = "".join(str(cds["sequence"]) for cds in selected_cds)
            target_contig = max(contigs, key=lambda contig: len(contig["sequence"]))
            target_contig["sequence"] += insert
            for cds in selected_cds:
                target_contig["parts"].append(
                    {
                        "role": "hypothetical_hgt_intended",
                        "parent_id": hgt_genome["id"],
                        "parent_record_id": hgt_genome["record_id"],
                        "start": cds["start"],
                        "end": cds["end"],
                        "length_bp": len(cds["sequence"]),
                        "sequence_sha256": _sha256_bytes(cds["sequence"].encode()),
                    }
                )
            _add_sample(
                samples,
                source,
                contigs,
                inputs_dir,
                sample_id=sample_id,
                scenario=f"hgt_like_{gene_count}cds",
                source_total_bp=source_total_bp,
                nominal_contamination=0.0,
                matched_control_sample_id=clean_full_id,
                interpretation="hypothetical_hgt_control_not_a_biological_claim",
            )

        related_row, divergent_row = _viral_donors(source, source_rows)
        relation_rows = {
            "viral_related": related_row,
            "viral_divergent": divergent_row,
        }
        for relation, donor_row in relation_rows.items():
            donor_id = str(donor_row["source_genome_id"])
            donor_genome = genome_cache.get(donor_id)
            if donor_genome is None:
                donor_genome = _load_genome(donor_row, source_manifest)
                genome_cache[donor_id] = donor_genome
            _add_mixture_series(
                samples,
                source,
                source_genome,
                donor_genome,
                inputs_dir,
                relation,
                VIRAL_LEVELS,
                fragment_bp,
                clean_full_id,
                "known_viral_donor_bp_is_contamination",
            )

        for cellular_genome in cellular_genomes:
            donor_type = str(cellular_genome["row"]["role"]).removesuffix(
                "_contamination_donor"
            )
            _add_mixture_series(
                samples,
                source,
                source_genome,
                cellular_genome,
                inputs_dir,
                f"cellular_{donor_type}_contigs",
                CELLULAR_LEVELS,
                fragment_bp,
                clean_full_id,
                "known_cellular_donor_bp_is_contamination",
            )

        chimera_genome = hgt_genome
        scenario = "cellular_eukaryotic_chimera_025"
        sample_id = f"{source_id}__{scenario}"
        contigs, base_bp, _ = _mixed_sample_parts(
            source_genome,
            chimera_genome,
            sample_id,
            CHIMERA_LEVEL,
            fragment_bp,
        )
        target_contigs = [
            contig
            for contig in contigs
            if all(part["role"] == "intended_target" for part in contig["parts"])
        ]
        donor_contigs = [
            contig
            for contig in contigs
            if all(part["role"] == "contaminant" for part in contig["parts"])
        ]
        matched_contigs = _copy_contigs(target_contigs)
        longest = max(target_contigs, key=lambda contig: len(contig["sequence"]))
        for donor_contig in donor_contigs:
            longest["sequence"] += donor_contig["sequence"]
            longest["parts"].extend(donor_contig["parts"])
        matched_id = clean_full_id
        if base_bp != source_total_bp:
            matched_id = f"{sample_id}__matched_clean"
            _add_sample(
                samples,
                source,
                matched_contigs,
                inputs_dir,
                sample_id=matched_id,
                scenario=f"clean_matched_{scenario}",
                source_total_bp=source_total_bp,
                nominal_contamination=0.0,
                matched_control_sample_id=matched_id,
                interpretation="matched_source_retention_control",
            )
        _add_sample(
            samples,
            source,
            target_contigs,
            inputs_dir,
            sample_id=sample_id,
            scenario=scenario,
            source_total_bp=source_total_bp,
            nominal_contamination=CHIMERA_LEVEL,
            matched_control_sample_id=matched_id,
            interpretation="known_cellular_donor_bp_on_chimeric_contig_is_contamination",
        )

    samples.sort(key=lambda row: row["sample_id"])
    for row in samples:
        _validate_sample(row, split)
    used_parent_ids = {
        parent_id for row in samples for parent_id in row["all_parent_ids"]
    }
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "profile_sha256": profile["profile_sha256"],
        "split": split,
        "truth_scope": "source_retention_not_absolute_genome_completeness",
        "parents": [
            {
                "parent_id": genome_cache[parent_id]["id"],
                "record_id": genome_cache[parent_id]["record_id"],
                "split": genome_cache[parent_id]["row"]["split"],
                "role": genome_cache[parent_id]["row"].get("role", "viral_source"),
                "fna_path": genome_cache[parent_id]["fna_path"],
                "genbank_path": genome_cache[parent_id]["genbank_path"],
                "fna_sha256": genome_cache[parent_id]["fna_sha256"],
                "genbank_sha256": genome_cache[parent_id]["genbank_sha256"],
                "sequence_sha256": genome_cache[parent_id]["sequence_sha256"],
            }
            for parent_id in sorted(used_parent_ids)
        ],
        "samples": samples,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    profile_path.write_text(json.dumps(profile, indent=2, sort_keys=True) + "\n")
    return manifest_path


def _read_summary(path: Path, sample_id: str) -> Dict[str, str]:
    if not path.is_file():
        raise ValueError(f"Missing GVClass summary: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or rows[0].get("query") != sample_id:
        raise ValueError(
            f"GVClass summary is not a valid one-query result for {sample_id}"
        )
    return rows[0]


def _run_settings(
    manifest_path: Path,
    sample: Dict[str, Any],
    database: Path,
    threads: int,
    mode: str,
) -> Dict[str, Any]:
    if threads <= 0:
        raise ValueError("Thread count must be positive")
    repo_root = Path(__file__).resolve().parents[1]
    database_files = [
        database / "DB_VERSION",
        database / "RESOURCE_METADATA.json",
        database / "contamination" / "model.joblib",
        database / "completeness" / "config.json",
        database / "completeness" / "model.joblib",
    ]
    missing_database_files = [path for path in database_files if not path.is_file()]
    if missing_database_files:
        raise ValueError(
            "Database identity files are missing: "
            + ", ".join(str(path) for path in missing_database_files)
        )
    database_hashes = {
        str(path.relative_to(database)): sha256_file(path) for path in database_files
    }
    code_files = [
        repo_root / "gvclass",
        repo_root / "scripts" / "_contamination_features.py",
        Path(__file__).resolve(),
        *sorted((repo_root / "src").rglob("*.py")),
    ]
    code_hashes = {
        str(path.relative_to(repo_root)): sha256_file(path) for path in code_files
    }
    settings = {
        "manifest_sha256": sha256_file(manifest_path),
        "sample_id": sample["sample_id"],
        "sample_input_sha256": sample["input_sha256"],
        "source_split": sample["source_split"],
        "database": str(database.resolve()),
        "database_hashes": database_hashes,
        "threads": threads,
        "mode": mode,
        "sensitive": True,
        "completeness_mode": COMPLETENESS_MODE,
        "code_hashes": code_hashes,
        "config_sha256": sha256_file(repo_root / "config" / "gvclass_config.yaml"),
        "pixi_lock_sha256": sha256_file(repo_root / "pixi.lock"),
        "model_features": CONTAMINATION_MODEL_FEATURES,
    }
    settings["settings_sha256"] = _sha256_bytes(_canonical_json(settings).encode())
    return settings


def _equal_number(left: Any, right: Any, tolerance: float = 1e-8) -> bool:
    try:
        left_number = float(left)
        right_number = float(right)
    except (TypeError, ValueError):
        return str(left) == str(right)
    if math.isnan(left_number) and math.isnan(right_number):
        return True
    return math.isclose(left_number, right_number, rel_tol=0.0, abs_tol=tolerance)


def _extract_artifacts(archive: Path, sample_id: str, destination: Path) -> Path:
    wanted = ("/hmmout/", "/stats/", "/blastp_out/", "/query_faa/")
    with tarfile.open(archive, "r:gz") as handle:
        members = []
        for member in handle.getmembers():
            member_path = Path(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise ValueError(f"Unsafe archive member: {member.name}")
            if any(token in member.name for token in wanted) or member.name.endswith(
                f"/{sample_id}.summary.tab"
            ):
                members.append(member)
        handle.extractall(destination, members=members)
    query_root = destination / sample_id
    if not query_root.is_dir():
        raise ValueError(f"Archive has no query root for {sample_id}")
    return query_root


def _audit_metadata(
    result_dir: Path,
    sample_id: str,
    mode_fast: bool,
    runtime_settings: Dict[str, Any],
) -> Dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[1]
    original_audit = result_dir / ORIGINAL_AUDIT_FILENAME
    return {
        "feature_audit_version": "mode_v2",
        "mode_fast": mode_fast,
        "runtime_settings_sha256": runtime_settings["settings_sha256"],
        "runtime_profile_sha256": sha256_file(result_dir / "quality_run_profile.json"),
        "runtime_artifact_sha256": {
            "gvclass_summary.tsv": sha256_file(result_dir / "gvclass_summary.tsv"),
            f"{sample_id}.tar.gz": sha256_file(result_dir / f"{sample_id}.tar.gz"),
        },
        "audit_code_sha256": {
            "scripts/_contamination_features.py": sha256_file(
                repo_root / "scripts" / "_contamination_features.py"
            ),
            "scripts/benchmark_quality_nucleotide.py": sha256_file(
                Path(__file__).resolve()
            ),
        },
        "original_audit_sha256": (
            sha256_file(original_audit) if original_audit.is_file() else None
        ),
    }


def _audit_matches(
    audit_path: Path,
    sample_id: str,
    expected_metadata: Dict[str, Any],
) -> bool:
    try:
        audit = json.loads(audit_path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return False
    return (
        audit.get("sample_id") == sample_id
        and audit.get("observed_subset_parity") == "passed"
        and audit.get("runtime_prediction_parity") == "passed"
        and all(audit.get(key) == value for key, value in expected_metadata.items())
    )


def _persist_feature_audit(
    result_dir: Path,
    sample: Dict[str, Any],
    database: Path,
    final_runtime_summary: Dict[str, str],
    mode_fast: bool,
    runtime_settings: Dict[str, Any],
) -> Dict[str, Any]:
    sample_id = str(sample["sample_id"])
    archive = result_dir / f"{sample_id}.tar.gz"
    if not archive.is_file():
        raise ValueError(f"Missing query archive for {sample_id}")
    with tempfile.TemporaryDirectory(prefix="gvclass_quality_") as temporary:
        query_root = _extract_artifacts(archive, sample_id, Path(temporary))
        query_fna_dir = query_root / "query_fna"
        query_fna_dir.mkdir(exist_ok=True)
        shutil.copyfile(sample["input_path"], query_fna_dir / f"{sample_id}.fna")
        persisted_runtime_summary = _read_summary(
            query_root / f"{sample_id}.summary.tab", sample_id
        )
        recomputed = compute_query_features(
            sample_id,
            query_root,
            database,
            sensitive_mode=True,
            completeness_mode=COMPLETENESS_MODE,
            mode_fast=mode_fast,
        )
    vector = contamination_feature_vector(recomputed)
    observed_subset = {
        name: persisted_runtime_summary[name]
        for name in CONTAMINATION_MODEL_FEATURES
        if name in persisted_runtime_summary and persisted_runtime_summary[name] != ""
    }
    mismatches = [
        name
        for name, value in observed_subset.items()
        if not _equal_number(value, vector[name], tolerance=0.011)
    ]
    if mismatches:
        raise ValueError(
            "Recomputed feature values differ from persisted runtime fields: "
            + ", ".join(mismatches)
        )
    if not _equal_number(
        final_runtime_summary.get("estimated_contamination"),
        recomputed.get("estimated_contamination"),
        tolerance=0.011,
    ):
        raise ValueError("Recomputed model estimate differs from the runtime estimate")
    audit = {
        **_audit_metadata(result_dir, sample_id, mode_fast, runtime_settings),
        "sample_id": sample_id,
        "runtime_persisted_summary": persisted_runtime_summary,
        "runtime_observed_feature_subset": observed_subset,
        "reconstructed_model_feature_vector": vector,
        "recomputed_cellular_coherent_bp_fraction": recomputed[
            "cellular_coherent_bp_fraction"
        ],
        "model_feature_order": CONTAMINATION_MODEL_FEATURES,
        "observed_subset_parity": "passed",
        "runtime_prediction_parity": "passed",
        "limitation": (
            "The query archive does not store all 21 raw model inputs. The complete "
            "vector is reconstructed through scripts._contamination_features; parity "
            "checks for persisted inputs and the saved runtime prediction use the "
            "0.011 runtime-rounding tolerance."
        ),
    }
    (result_dir / AUDIT_FILENAME).write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    return audit


def run_one(
    manifest_path: Path,
    index: int,
    database: Path,
    threads: int,
    mode: str,
    force: bool = False,
    allow_final_test: bool = False,
) -> Path:
    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest["split"] == "final_test" and not allow_final_test:
        raise ValueError("Final-test execution requires --allow-final-test")
    samples = manifest["samples"]
    if index < 0 or index >= len(samples):
        raise IndexError(f"Index {index} outside 0..{len(samples) - 1}")
    sample = samples[index]
    _validate_sample(sample, manifest["split"])
    database = database.resolve()
    settings = _run_settings(manifest_path, sample, database, threads, mode)
    results_root = manifest_path.parent / "results"
    result_dir = results_root / sample["sample_id"]
    profile_path = result_dir / "quality_run_profile.json"
    combined_summary = result_dir / "gvclass_summary.tsv"
    if result_dir.exists():
        cached_profile = (
            json.loads(profile_path.read_text()) if profile_path.exists() else {}
        )
        cache_valid = False
        try:
            _read_summary(combined_summary, sample["sample_id"])
            cache_valid = _audit_matches(
                result_dir / AUDIT_FILENAME,
                sample["sample_id"],
                _audit_metadata(
                    result_dir,
                    sample["sample_id"],
                    mode == "fast",
                    cached_profile,
                ),
            )
        except (FileNotFoundError, KeyError, ValueError, json.JSONDecodeError):
            cache_valid = False
        if (
            cache_valid
            and cached_profile.get("settings_sha256") == settings["settings_sha256"]
        ):
            return result_dir
        if not force:
            raise RuntimeError(
                f"Cached result for {sample['sample_id']} is invalid or has different settings"
            )
        shutil.rmtree(result_dir)

    result_dir.parent.mkdir(parents=True, exist_ok=True)
    repo_root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(
        prefix=f"quality_{index:05d}_", dir=manifest_path.parent
    ) as temporary:
        query_dir = Path(temporary) / "query"
        query_dir.mkdir()
        query_path = query_dir / f"{sample['sample_id']}.fna"
        query_path.symlink_to(Path(sample["input_path"]).resolve())
        command = [
            str(repo_root / "gvclass"),
            str(query_dir),
            "--output-dir",
            str(result_dir),
            "--database",
            str(database),
            "--threads",
            str(threads),
            "--max-workers",
            "1",
            "--threads-per-worker",
            str(threads),
            "--completeness-mode",
            COMPLETENESS_MODE,
            "--sensitive",
            "--min-length",
            "0",
            "--plain-output",
            "--mode-fast" if mode == "fast" else "--extended",
        ]
        subprocess.run(command, cwd=repo_root, check=True)

    runtime_summary = _read_summary(combined_summary, sample["sample_id"])
    profile_path.write_text(json.dumps(settings, indent=2, sort_keys=True) + "\n")
    _persist_feature_audit(
        result_dir,
        sample,
        database,
        runtime_summary,
        mode_fast=mode == "fast",
        runtime_settings=settings,
    )
    return result_dir


def audit_one(
    manifest_path: Path,
    index: int,
    database: Path,
    threads: int,
    mode: str,
    force: bool = False,
    allow_final_test: bool = False,
) -> Path:
    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest["split"] == "final_test" and not allow_final_test:
        raise ValueError("Final-test audit requires --allow-final-test")
    if index < 0 or index >= len(manifest["samples"]):
        raise IndexError(f"Index {index} outside 0..{len(manifest['samples']) - 1}")
    sample = manifest["samples"][index]
    _validate_sample(sample, manifest["split"])
    result_dir = manifest_path.parent / "results" / sample["sample_id"]
    runtime_summary = _read_summary(
        result_dir / "gvclass_summary.tsv", sample["sample_id"]
    )
    profile_path = result_dir / "quality_run_profile.json"
    if not profile_path.is_file():
        raise ValueError(f"Missing original run profile for {sample['sample_id']}")
    runtime_settings = json.loads(profile_path.read_text())
    stored_settings_hash = runtime_settings.get("settings_sha256")
    settings_body = {
        key: value
        for key, value in runtime_settings.items()
        if key != "settings_sha256"
    }
    if stored_settings_hash != _sha256_bytes(_canonical_json(settings_body).encode()):
        raise ValueError(f"Invalid original run profile for {sample['sample_id']}")

    current_settings = _run_settings(
        manifest_path, sample, database.resolve(), threads, mode
    )
    common_keys = set(current_settings).difference({"code_hashes", "settings_sha256"})
    changed_settings = sorted(
        key
        for key in common_keys
        if runtime_settings.get(key) != current_settings.get(key)
    )
    original_code = runtime_settings.get("code_hashes", {})
    current_code = current_settings["code_hashes"]
    original_runtime_code = {
        key: value
        for key, value in original_code.items()
        if key == "gvclass" or key.startswith("src/")
    }
    current_runtime_code = {
        key: value
        for key, value in current_code.items()
        if key == "gvclass" or key.startswith("src/")
    }
    if changed_settings or original_runtime_code != current_runtime_code:
        details = changed_settings or ["src/**/*.py or gvclass hash"]
        raise ValueError(
            f"Runtime settings changed for {sample['sample_id']}: " + ", ".join(details)
        )

    audit_path = result_dir / AUDIT_FILENAME
    expected_metadata = _audit_metadata(
        result_dir, sample["sample_id"], mode == "fast", runtime_settings
    )
    if audit_path.exists() and not force:
        if _audit_matches(audit_path, sample["sample_id"], expected_metadata):
            return audit_path
        raise RuntimeError(
            f"Existing revised audit is invalid for {sample['sample_id']}"
        )
    _persist_feature_audit(
        result_dir,
        sample,
        database.resolve(),
        runtime_summary,
        mode_fast=mode == "fast",
        runtime_settings=runtime_settings,
    )
    return audit_path


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _q90(values: Sequence[float]) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    position = 0.9 * (len(ordered) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def metric_row(rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    n = len(rows)
    output: Dict[str, Any] = {"n": n}
    for prefix, estimate_key, truth_key in (
        ("completeness", "estimated_completeness", "true_source_retention_pct"),
        ("contamination", "estimated_contamination", "true_contamination_pct"),
    ):
        finite_rows = [row for row in rows if _finite(row.get(estimate_key))]
        errors = [
            float(row[estimate_key]) - float(row[truth_key]) for row in finite_rows
        ]
        output[f"{prefix}_finite_coverage"] = (
            len(finite_rows) / n if n else float("nan")
        )
        output[f"{prefix}_nan_rate"] = (
            1.0 - (len(finite_rows) / n) if n else float("nan")
        )
        output[f"{prefix}_mae"] = (
            sum(abs(error) for error in errors) / len(errors)
            if errors
            else float("nan")
        )
        output[f"{prefix}_signed_bias"] = (
            sum(errors) / len(errors) if errors else float("nan")
        )
        output[f"{prefix}_q90_error"] = _q90([abs(error) for error in errors])

    clean = [row for row in rows if float(row["true_contamination_pct"]) == 0.0]
    contaminated = [row for row in rows if float(row["true_contamination_pct"]) >= 20.0]
    output["contamination_fpr_at_10"] = (
        sum(
            _finite(row.get("estimated_contamination"))
            and float(row["estimated_contamination"]) >= 10.0
            for row in clean
        )
        / len(clean)
        if clean
        else float("nan")
    )
    output["contamination_sensitivity_at_10_truth_ge_20"] = (
        sum(
            _finite(row.get("estimated_contamination"))
            and float(row["estimated_contamination"]) >= 10.0
            for row in contaminated
        )
        / len(contaminated)
        if contaminated
        else float("nan")
    )
    output["unsupported_contaminated_count"] = sum(
        not _finite(row.get("estimated_contamination")) for row in contaminated
    )
    return output


def _write_tsv(rows: Sequence[Dict[str, Any]], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_results(
    manifest_path: Path, allow_final_test: bool = False
) -> Tuple[Path, Path]:
    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest["split"] == "final_test" and not allow_final_test:
        raise ValueError("Final-test summary requires --allow-final-test")
    predictions: List[Dict[str, Any]] = []
    results_root = manifest_path.parent / "results"
    manifest_hash = sha256_file(manifest_path)
    for sample in manifest["samples"]:
        _validate_sample(sample, manifest["split"])
        result_dir = results_root / sample["sample_id"]
        summary = _read_summary(result_dir / "gvclass_summary.tsv", sample["sample_id"])
        audit_path = result_dir / AUDIT_FILENAME
        profile_path = result_dir / "quality_run_profile.json"
        archive_path = result_dir / f"{sample['sample_id']}.tar.gz"
        if (
            not audit_path.is_file()
            or not profile_path.is_file()
            or not archive_path.is_file()
        ):
            raise ValueError(f"Incomplete benchmark result for {sample['sample_id']}")
        runtime_settings = json.loads(profile_path.read_text())
        settings_body = {
            key: value
            for key, value in runtime_settings.items()
            if key != "settings_sha256"
        }
        profile_is_bound = (
            runtime_settings.get("settings_sha256")
            == _sha256_bytes(_canonical_json(settings_body).encode())
            and runtime_settings.get("manifest_sha256") == manifest_hash
            and runtime_settings.get("sample_id") == sample["sample_id"]
            and runtime_settings.get("sample_input_sha256") == sample["input_sha256"]
            and runtime_settings.get("source_split") == sample["source_split"]
            and runtime_settings.get("mode") in {"fast", "extended"}
        )
        if not profile_is_bound or not _audit_matches(
            audit_path,
            sample["sample_id"],
            _audit_metadata(
                result_dir,
                sample["sample_id"],
                runtime_settings["mode"] == "fast",
                runtime_settings,
            ),
        ):
            raise ValueError(f"Audit is not bound to result for {sample['sample_id']}")
        audit = json.loads(audit_path.read_text())
        legacy = audit["runtime_persisted_summary"]
        row = {
            "sample_id": sample["sample_id"],
            "source_genome_id": sample["source_genome_id"],
            "source_family": sample["source_family"],
            "source_split": sample["source_split"],
            "scenario": sample["scenario"],
            "all_parent_ids": ",".join(sample["all_parent_ids"]),
            "true_source_retention_pct": sample["true_source_retention_pct"],
            "true_contamination_pct": sample["true_contamination_pct"],
            "estimated_completeness": summary.get("estimated_completeness", "NaN"),
            "completeness_query_support": summary.get("completeness_query_support", ""),
            "estimated_contamination": summary.get("estimated_contamination", "NaN"),
            "contamination_query_support": summary.get(
                "contamination_query_support", ""
            ),
            "run_settings_sha256": runtime_settings["settings_sha256"],
            "cellular_coherent_bp_fraction": audit[
                "recomputed_cellular_coherent_bp_fraction"
            ],
            "feature_observed_subset_parity": audit["observed_subset_parity"],
            "feature_runtime_prediction_parity": audit["runtime_prediction_parity"],
        }
        row.update(audit["reconstructed_model_feature_vector"])
        row.update(
            {f"runtime_{key}": value for key, value in legacy.items() if key != "query"}
        )
        predictions.append(row)

    metrics: List[Dict[str, Any]] = []
    for group_by in ("scenario", "source_genome_id", "source_family"):
        values = sorted({str(row[group_by]) for row in predictions})
        for value in values:
            group = [row for row in predictions if str(row[group_by]) == value]
            metrics.append(
                {"group_by": group_by, "group_value": value, **metric_row(group)}
            )
    predictions_path = manifest_path.parent / "quality_predictions.tsv"
    metrics_path = manifest_path.parent / "quality_metrics.tsv"
    _write_tsv(predictions, predictions_path)
    _write_tsv(metrics, metrics_path)
    return predictions_path, metrics_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--source-manifest", type=Path, required=True)
    prepare.add_argument("--donor-manifest", type=Path, required=True)
    prepare.add_argument("--work-dir", type=Path, required=True)
    prepare.add_argument(
        "--split", choices=("development", "final_test"), default="development"
    )
    prepare.add_argument("--fragment-bp", type=int, default=20_000)
    prepare.add_argument("--force", action="store_true")
    prepare.add_argument("--allow-final-test", action="store_true")

    run = subparsers.add_parser("run")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--index", type=int, required=True)
    run.add_argument("--database", type=Path, required=True)
    run.add_argument("--threads", type=int, required=True)
    run.add_argument("--mode", choices=("fast", "extended"), required=True)
    run.add_argument("--audit-only", action="store_true")
    run.add_argument("--force", action="store_true")
    run.add_argument("--allow-final-test", action="store_true")

    summarize = subparsers.add_parser("summarize")
    summarize.add_argument("--manifest", type=Path, required=True)
    summarize.add_argument("--allow-final-test", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "prepare":
        manifest = prepare_dataset(
            args.source_manifest,
            args.donor_manifest,
            args.work_dir,
            args.split,
            fragment_bp=args.fragment_bp,
            force=args.force,
            allow_final_test=args.allow_final_test,
        )
        count = len(json.loads(manifest.read_text())["samples"])
        print(f"Prepared {count} samples: {manifest}")
    elif args.command == "run":
        operation = audit_one if args.audit_only else run_one
        result = operation(
            args.manifest,
            args.index,
            args.database,
            args.threads,
            args.mode,
            force=args.force,
            allow_final_test=args.allow_final_test,
        )
        print(result)
    else:
        predictions, metrics = summarize_results(
            args.manifest, allow_final_test=args.allow_final_test
        )
        print(predictions)
        print(metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
