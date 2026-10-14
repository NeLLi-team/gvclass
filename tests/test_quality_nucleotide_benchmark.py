import json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from scripts import benchmark_quality_nucleotide as benchmark
from scripts._contamination_features import compute_query_features
from src.utils.common import sha256_file


def _genome_files(root: Path, genome_id: str, sequence: str) -> dict:
    fna = root / f"{genome_id}.fna"
    genbank = root / f"{genome_id}.gb"
    record = SeqRecord(Seq(sequence), id=genome_id, name=genome_id, description="test")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(FeatureLocation(start, start + 90), type="CDS")
        for start in (150, 300, 450, 600, 750)
        if start + 90 <= len(sequence)
    ]
    SeqIO.write(record, fna, "fasta")
    SeqIO.write(record, genbank, "genbank")
    return {
        "fna": str(fna),
        "genbank": str(genbank),
        "fna_sha256": sha256_file(fna),
        "genbank_sha256": sha256_file(genbank),
        "sequence_sha256": benchmark._sha256_bytes(sequence.encode()),
    }


@pytest.fixture
def prepared(tmp_path):
    source_rows = []
    for index, genome_id in enumerate(("V1", "V2", "V3", "V4")):
        files = _genome_files(tmp_path, genome_id, "ACGT" * 250)
        source_rows.append(
            {
                "source_genome_id": genome_id,
                "panel": f"family{index // 2}",
                "taxonomy": ["Viruses", "NCLDV", f"family{index // 2}"],
                "split": "development",
                "length_bp": 1000,
                "fna_path": files["fna"],
                "genbank_path": files["genbank"],
                "fna_sha256": files["fna_sha256"],
                "genbank_sha256": files["genbank_sha256"],
                "sequence_sha256": files["sequence_sha256"],
            }
        )
    final_files = _genome_files(tmp_path, "VF", "TGCA" * 250)
    source_rows.append(
        {
            "source_genome_id": "VF",
            "panel": "final_family",
            "taxonomy": ["Viruses", "NCLDV", "final_family"],
            "split": "final_test",
            "length_bp": 1000,
            "fna_path": final_files["fna"],
            "genbank_path": final_files["genbank"],
            "fna_sha256": final_files["fna_sha256"],
            "genbank_sha256": final_files["genbank_sha256"],
            "sequence_sha256": final_files["sequence_sha256"],
        }
    )
    source_manifest = tmp_path / "sources.json"
    source_manifest.write_text(json.dumps(source_rows))

    donors = []
    for genome_id, role, base in (
        ("BDEV", "bacterial_contamination_donor", "A"),
        ("EDEV", "eukaryotic_contamination_donor", "C"),
        ("BFINAL", "bacterial_contamination_donor", "G"),
        ("EFINAL", "eukaryotic_contamination_donor", "T"),
    ):
        files = _genome_files(tmp_path, genome_id, base * 900)
        split = "development" if genome_id.endswith("DEV") else "final_test"
        donors.append(
            {
                "accession": genome_id,
                "role": role,
                "split": split,
                "paths": {"fasta": files["fna"], "genbank": files["genbank"]},
                "sha256": {
                    "fasta_file": files["fna_sha256"],
                    "genbank_file": files["genbank_sha256"],
                    "normalized_sequence": files["sequence_sha256"],
                },
                "checks": {"sequence_valid": True},
            }
        )
    donor_manifest = tmp_path / "donors.json"
    donor_manifest.write_text(
        json.dumps({"donor_summary": {"all_checks_pass": True}, "donors": donors})
    )
    work_dir = tmp_path / "prepared"
    manifest_path = benchmark.prepare_dataset(
        source_manifest, donor_manifest, work_dir, "development", fragment_bp=200
    )
    assert (
        benchmark.prepare_dataset(
            source_manifest, donor_manifest, work_dir, "development", fragment_bp=200
        )
        == manifest_path
    )
    return manifest_path, source_manifest, donor_manifest


def test_prepare_has_exact_truth_and_no_repeated_parent_span(prepared):
    manifest_path, source_manifest, _ = prepared
    manifest = json.loads(manifest_path.read_text())
    samples = {row["sample_id"]: row for row in manifest["samples"]}

    for retention in benchmark.RETENTION_LEVELS:
        sample = samples[f"V1__clean_retention_{round(retention * 100):03d}"]
        assert sample["base_retained_bp"] == round(1000 * retention)
        assert sample["contamination_truth_bp"] == 0
        assert max(contig["length_bp"] for contig in sample["contigs"]) <= 200
    assert samples["V1__clean_retention_100"]["cds_truncating_boundaries"]

    mixed = samples["V1__cellular_bacterial_contigs_095"]
    assert mixed["true_contamination_pct"] == pytest.approx(
        100 * mixed["contamination_truth_bp"] / mixed["output_bp"]
    )
    assert mixed["contamination_truth_bp"] <= 900
    intervals = []
    for contig in mixed["contigs"]:
        intervals.extend(
            (part["parent_id"], part["start"], part["end"]) for part in contig["parts"]
        )
    for index, (parent, start, end) in enumerate(intervals):
        assert not any(
            parent == other_parent and start < other_end and other_start < end
            for other_parent, other_start, other_end in intervals[:index]
        )

    hgt = samples["V1__hgt_like_5cds"]
    assert hgt["foreign_sequence_bp"] > 0
    assert hgt["contamination_truth_bp"] == 0
    chimera = samples["V1__cellular_eukaryotic_chimera_025"]
    assert chimera["contamination_truth_bp"] > 0
    assert any(len(contig["parts"]) > 1 for contig in chimera["contigs"])

    panels = {
        row["source_genome_id"]: row["panel"]
        for row in json.loads(source_manifest.read_text())
    }
    for sample in samples.values():
        if not sample["scenario"].startswith(("viral_related_", "viral_divergent_")):
            continue
        donor_id = next(
            parent
            for parent in sample["all_parent_ids"]
            if parent != sample["source_genome_id"]
        )
        same_panel = panels[donor_id] == sample["source_family"]
        assert same_panel == sample["scenario"].startswith("viral_related_")


def test_prepare_rejects_source_or_donor_split_leakage(prepared, tmp_path):
    _, source_manifest, donor_manifest = prepared
    sources = json.loads(source_manifest.read_text())
    unsafe_sources = [{**sources[0], "source_genome_id": "../escape"}, *sources[1:]]
    unsafe_source_path = tmp_path / "unsafe_sources.json"
    unsafe_source_path.write_text(json.dumps(unsafe_sources))
    unsafe_work_dir = tmp_path / "unsafe"
    with pytest.raises(ValueError, match="Unsafe source genome ID"):
        benchmark.prepare_dataset(
            unsafe_source_path, donor_manifest, unsafe_work_dir, "development"
        )
    assert not unsafe_work_dir.exists()

    leaked_sources = [*sources, {**sources[0], "split": "final_test"}]
    leaked_source_path = tmp_path / "leaked_sources.json"
    leaked_source_path.write_text(json.dumps(leaked_sources))
    with pytest.raises(ValueError, match="Source genome IDs overlap"):
        benchmark.prepare_dataset(
            leaked_source_path, donor_manifest, tmp_path / "bad_source", "development"
        )

    donor_data = json.loads(donor_manifest.read_text())
    donor_data["donors"].append({**donor_data["donors"][0], "split": "final_test"})
    leaked_donor_path = tmp_path / "leaked_donors.json"
    leaked_donor_path.write_text(json.dumps(donor_data))
    with pytest.raises(ValueError, match="Cellular donor IDs overlap"):
        benchmark.prepare_dataset(
            source_manifest, leaked_donor_path, tmp_path / "bad_donor", "development"
        )


def test_metrics_keep_nan_in_denominators_and_count_missed():
    rows = [
        {
            "true_source_retention_pct": 50,
            "estimated_completeness": "NaN",
            "true_contamination_pct": 0,
            "estimated_contamination": "NaN",
        },
        {
            "true_source_retention_pct": 50,
            "estimated_completeness": 40,
            "true_contamination_pct": 25,
            "estimated_contamination": 15,
        },
        {
            "true_source_retention_pct": 50,
            "estimated_completeness": "NaN",
            "true_contamination_pct": 25,
            "estimated_contamination": "NaN",
        },
    ]
    metrics = benchmark.metric_row(rows)
    assert metrics["contamination_finite_coverage"] == pytest.approx(1 / 3)
    assert metrics["contamination_nan_rate"] == pytest.approx(2 / 3)
    assert metrics["contamination_fpr_at_10"] == 0
    assert metrics["contamination_sensitivity_at_10_truth_ge_20"] == 0.5
    assert metrics["unsupported_contaminated_count"] == 1


def test_feature_reconstruction_passes_explicit_pipeline_mode(tmp_path):
    class FakeSummarizer:
        completeness_mode = "novelty-aware"
        sensitive_mode = True

        def summarize_query_full(
            self, query_id, query_output_dir, tree_nn_results, mode_fast=True
        ):
            self.mode_fast = mode_fast
            return {
                **{name: 0 for name in benchmark.CONTAMINATION_MODEL_FEATURES},
                "contamination_feature_version": "feature_v2_unvalidated",
            }

    summarizer = FakeSummarizer()
    compute_query_features(
        "query",
        tmp_path,
        tmp_path,
        summarizer=summarizer,
        mode_fast=False,
    )
    assert summarizer.mode_fast is False


def test_resume_rejects_changed_run_settings(prepared, tmp_path, monkeypatch):
    manifest_path, _, _ = prepared
    manifest = json.loads(manifest_path.read_text())
    sample = manifest["samples"][0]
    database = tmp_path / "database"
    for relative in (
        "DB_VERSION",
        "RESOURCE_METADATA.json",
        "contamination/model.joblib",
        "completeness/config.json",
        "completeness/model.joblib",
    ):
        path = database / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative)
    result_dir = manifest_path.parent / "results" / sample["sample_id"]
    result_dir.mkdir(parents=True)
    (result_dir / "gvclass_summary.tsv").write_text(f"query\n{sample['sample_id']}\n")
    archive_path = result_dir / f"{sample['sample_id']}.tar.gz"
    archive_path.write_bytes(b"archive")
    settings = benchmark._run_settings(manifest_path, sample, database, 2, "fast")
    (result_dir / "quality_run_profile.json").write_text(json.dumps(settings))
    (result_dir / benchmark.AUDIT_FILENAME).write_text(
        json.dumps(
            {
                **benchmark._audit_metadata(
                    result_dir, sample["sample_id"], True, settings
                ),
                "sample_id": sample["sample_id"],
                "observed_subset_parity": "passed",
                "runtime_prediction_parity": "passed",
            }
        )
    )

    assert benchmark.run_one(manifest_path, 0, database, 2, "fast") == result_dir
    legacy_settings = {
        key: value
        for key, value in settings.items()
        if key not in {"config_sha256", "settings_sha256"}
    }
    legacy_settings["settings_sha256"] = benchmark._sha256_bytes(
        benchmark._canonical_json(legacy_settings).encode()
    )
    (result_dir / "quality_run_profile.json").write_text(json.dumps(legacy_settings))
    with pytest.raises(ValueError, match="config_sha256"):
        benchmark.audit_one(manifest_path, 0, database, 2, "fast")
    (result_dir / "quality_run_profile.json").write_text(json.dumps(settings))

    archive_path.write_bytes(b"changed")
    with pytest.raises(RuntimeError, match="invalid or has different settings"):
        benchmark.run_one(manifest_path, 0, database, 2, "fast")
    archive_path.write_bytes(b"archive")
    with pytest.raises(RuntimeError, match="different settings"):
        benchmark.run_one(manifest_path, 0, database, 3, "fast")

    original_sha256_file = benchmark.sha256_file
    config_path = (
        Path(benchmark.__file__).resolve().parents[1] / "config" / "gvclass_config.yaml"
    )

    def changed_config_hash(path):
        if Path(path).resolve() == config_path.resolve():
            return "0" * 64
        return original_sha256_file(path)

    monkeypatch.setattr(benchmark, "sha256_file", changed_config_hash)
    with pytest.raises(RuntimeError, match="different settings"):
        benchmark.run_one(manifest_path, 0, database, 2, "fast")
