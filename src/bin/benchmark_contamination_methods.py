"""Benchmark contamination detection methods on synthetic isolate-derived bins."""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import subprocess
import tarfile
import tempfile
from collections import Counter
from pathlib import Path
from statistics import mean
from typing import Any, Dict, List, Tuple

import pandas as pd
from Bio import SeqIO
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import average_precision_score, mean_absolute_error, roc_auc_score
from sklearn.model_selection import GroupKFold

from src.core.contamination_scoring import create_contamination_scorer

RANDOM_SEED = 42
BASE_FAMILIES = [
    "Iridoviridae",
    "Phycodnaviridae",
    "Mimiviridae",
    "Marseilleviridae",
    "Ascoviridae",
    "Poxviridae",
    "Asfarviridae",
]
SCENARIOS = [
    ("clean", 0.0),
    ("hgt_cellular_1gene", 0.0),
    ("hgt_cellular_5genes", 0.0),
    ("cellular_contig_10", 0.10),
    ("cellular_contig_25", 0.25),
    ("viral_divergent_10", 0.10),
    ("viral_divergent_25", 0.25),
    ("viral_related_10", 0.10),
    ("viral_related_25", 0.25),
]
METHOD3_FEATURES = [
    "contamination_score_v1",
    "contamination_cellular_signal_v1",
    "contamination_phage_signal_v1",
    "contamination_duplication_signal_v1",
    "contamination_viral_mixture_signal_v1",
    "contamination_nonviral_hit_fraction_v1",
    "order_dup",
    "gvog8_dup",
    "cellular_unique",
    "cellular_total",
    "phage_unique",
    "phage_total",
    "contigs",
    "suspicious_bp_fraction_v2",
    "suspicious_contig_count_v2",
    "estimated_completeness",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="benchmarking/contamination/contamination_benchmark")
    parser.add_argument("--input-fna-dir", default="benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10/fna")
    parser.add_argument("--metadata", default="benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10/metadata_table.tsv")
    parser.add_argument("--database", default="resources")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--force-regenerate", action="store_true")
    parser.add_argument("--force-rerun", action="store_true")
    return parser.parse_args()


def family_from_taxonomy(tax: str) -> str:
    parts = [p.strip() for p in tax.split(";")]
    return parts[6] if len(parts) > 6 else "NA"


def order_from_taxonomy(tax: str) -> str:
    parts = [p.strip() for p in tax.split(";")]
    return parts[5] if len(parts) > 5 else "NA"


def load_metadata(path: Path) -> List[Dict[str, Any]]:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        row["family"] = family_from_taxonomy(row["taxonomy"])
        row["order_name"] = order_from_taxonomy(row["taxonomy"])
        row["genome_size_bp"] = int(row["genome_size_bp"])
        row["num_contigs"] = int(row["num_contigs"])
    return rows


def select_base_genomes(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    selected = []
    for family in BASE_FAMILIES:
        family_rows = [
            row for row in rows
            if row["family"] == family
            and row["num_contigs"] == 1
            and 100_000 <= row["genome_size_bp"] <= 700_000
        ]
        family_rows.sort(key=lambda row: row["genome_size_bp"])
        selected.extend(family_rows[:2])
    if len(selected) != len(BASE_FAMILIES) * 2:
        raise RuntimeError(f"Expected {len(BASE_FAMILIES) * 2} base genomes, found {len(selected)}")
    return selected


def load_fna_records(path: Path) -> List[Any]:
    return list(SeqIO.parse(str(path), "fasta"))


def slice_sequence(seq: str, target_len: int, rng: random.Random) -> str:
    if target_len >= len(seq):
        return seq
    start = rng.randint(0, max(0, len(seq) - target_len))
    return seq[start:start + target_len]


def contamination_bp(base_bp: int, target_fraction: float) -> int:
    return max(2000, round((base_bp * target_fraction) / max(1e-9, 1.0 - target_fraction)))


def pick_related_donor(base: Dict[str, Any], all_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    family_pool = [row for row in all_rows if row["family"] == base["family"] and row["accession_number"] != base["accession_number"] and row["num_contigs"] == 1 and row["genome_size_bp"] >= 100_000]
    if family_pool:
        return sorted(family_pool, key=lambda row: row["genome_size_bp"], reverse=True)[0]
    order_pool = [row for row in all_rows if row["order_name"] == base["order_name"] and row["accession_number"] != base["accession_number"] and row["num_contigs"] == 1 and row["genome_size_bp"] >= 100_000]
    return sorted(order_pool, key=lambda row: row["genome_size_bp"], reverse=True)[0]


def pick_divergent_donor(base: Dict[str, Any], all_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    pool = [row for row in all_rows if row["order_name"] != base["order_name"] and row["num_contigs"] == 1 and row["genome_size_bp"] >= 100_000]
    return sorted(pool, key=lambda row: row["genome_size_bp"], reverse=True)[0]


def build_cellular_gene_pool() -> List[str]:
    gene_path = Path('.pixi/envs/default/lib/python3.11/site-packages/pyhmmer/tests/data/seqs/CP040672.1.genes_100.fna')
    return [str(record.seq) for record in SeqIO.parse(str(gene_path), 'fasta') if len(record.seq) >= 300]


def build_cellular_contig_source() -> str:
    contig_path = Path('.pixi/envs/default/lib/python3.11/site-packages/pyhmmer/tests/data/seqs/1390.SAMEA104415756.OFHT01000022.fna')
    return str(next(SeqIO.parse(str(contig_path), 'fasta')).seq)


def write_synthetic_bin(
    output_path: Path,
    base_row: Dict[str, Any],
    all_rows: List[Dict[str, Any]],
    base_records: List[Any],
    cellular_genes: List[str],
    cellular_contig: str,
    rng: random.Random,
    scenario: str,
    target_fraction: float,
    input_fna_dir: Path,
) -> Dict[str, Any]:
    records = [record[:] for record in base_records]
    base_bp = sum(len(record.seq) for record in records)
    truth_fraction = 0.0
    contamination_type = "clean"
    added_bp = 0

    if scenario.startswith("hgt_cellular"):
        contamination_type = "hgt_like"
        n_genes = 1 if scenario.endswith("1gene") else 5
        largest = max(records, key=lambda rec: len(rec.seq))
        pieces = [slice_sequence(gene, min(len(gene), 2000), rng) for gene in rng.sample(cellular_genes, n_genes)]
        hgt_block = ("N" * 200).join(pieces)
        largest.seq = largest.seq + ("N" * 200) + hgt_block
        added_bp = len(hgt_block)
    elif scenario.startswith("cellular_contig"):
        contamination_type = "cellular_contig"
        added_bp = contamination_bp(base_bp, target_fraction)
        seq = slice_sequence(cellular_contig, added_bp, rng)
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        records.append(SeqRecord(Seq(seq), id=f"cellular_contam_{scenario}", description=""))
        truth_fraction = (len(seq) / (base_bp + len(seq))) * 100.0
    elif scenario.startswith("viral_divergent") or scenario.startswith("viral_related"):
        contamination_type = "viral_contig_divergent" if "divergent" in scenario else "viral_contig_related"
        donor_row = pick_divergent_donor(base_row, all_rows) if "divergent" in scenario else pick_related_donor(base_row, all_rows)
        donor_records = load_fna_records(input_fna_dir / f"{donor_row['accession_number']}.fna")
        donor_record = max(donor_records, key=lambda rec: len(rec.seq))
        added_bp = contamination_bp(base_bp, target_fraction)
        seq = slice_sequence(str(donor_record.seq), added_bp, rng)
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        records.append(SeqRecord(Seq(seq), id=f"viral_contam_{scenario}", description=""))
        truth_fraction = (len(seq) / (base_bp + len(seq))) * 100.0
    else:
        contamination_type = "clean"

    SeqIO.write(records, str(output_path), 'fasta')
    return {
        "sample_id": output_path.stem,
        "base_accession": base_row["accession_number"],
        "base_family": base_row["family"],
        "base_order": base_row["order_name"],
        "scenario": scenario,
        "split": base_row["split"],
        "contamination_type": contamination_type,
        "true_contamination_pct": round(truth_fraction, 2),
        "added_bp": added_bp,
        "base_bp": base_bp,
    }


def generate_dataset(run_dir: Path, input_fna_dir: Path, metadata_rows: List[Dict[str, Any]], force: bool) -> Path:
    inputs_dir = run_dir / 'inputs'
    manifest_path = run_dir / 'synthetic_manifest.tsv'
    if manifest_path.exists() and not force:
        return manifest_path
    if inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    inputs_dir.mkdir(parents=True, exist_ok=True)

    cellular_genes = build_cellular_gene_pool()
    cellular_contig = build_cellular_contig_source()
    base_rows = select_base_genomes(metadata_rows)
    base_rows_by_family = {}
    for family in BASE_FAMILIES:
        pair = [row for row in base_rows if row['family'] == family]
        pair.sort(key=lambda row: row['genome_size_bp'], reverse=True)
        pair[0]['split'] = 'train'
        pair[1]['split'] = 'test'
        base_rows_by_family[family] = pair

    rng = random.Random(RANDOM_SEED)
    manifest_rows: List[Dict[str, Any]] = []
    for family in BASE_FAMILIES:
        for base_row in base_rows_by_family[family]:
            base_records = load_fna_records(input_fna_dir / f"{base_row['accession_number']}.fna")
            for scenario, target_fraction in SCENARIOS:
                sample_id = f"{base_row['accession_number']}__{scenario}"
                output_path = inputs_dir / f"{sample_id}.fna"
                manifest_rows.append(
                    write_synthetic_bin(
                        output_path=output_path,
                        base_row=base_row,
                        all_rows=metadata_rows,
                        base_records=base_records,
                        cellular_genes=cellular_genes,
                        cellular_contig=cellular_contig,
                        rng=rng,
                        scenario=scenario,
                        target_fraction=target_fraction,
                        input_fna_dir=input_fna_dir,
                    )
                )

    with manifest_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()), delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest_rows)
    return manifest_path


def run_gvclass(inputs_dir: Path, results_dir: Path, workers: int, force: bool) -> None:
    if results_dir.exists() and (results_dir / 'gvclass_summary.tsv').exists() and not force:
        return
    if results_dir.exists() and force:
        shutil.rmtree(results_dir)
    cmd = [
        './gvclass',
        str(inputs_dir),
        '-o', str(results_dir),
        '--threads', str(workers),
        '-j', str(workers),
        '--threads-per-worker', '1',
        '--completeness-mode', 'novelty-aware',
        '--plain-output',
    ]
    subprocess.run(cmd, check=True)


def parse_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def extract_primary_token(value: str) -> str:
    if not value:
        return ''
    token = value.split(',')[0].split(':')[0]
    return token.split('__', 1)[1] if '__' in token else token


def compute_method2_features(results_dir: Path, summary_rows: Dict[str, Dict[str, Any]], database: Path) -> Dict[str, Dict[str, Any]]:
    scorer = create_contamination_scorer(database)
    outputs: Dict[str, Dict[str, Any]] = {}
    for sample_id, row in summary_rows.items():
        tar_path = results_dir / f"{sample_id}.tar.gz"
        if not tar_path.exists():
            outputs[sample_id] = {"suspicious_bp_fraction_v2": 0.0, "suspicious_contig_count_v2": 0}
            continue
        primary_order = extract_primary_token(row.get('order', ''))
        primary_family = extract_primary_token(row.get('family', ''))
        with tempfile.TemporaryDirectory(prefix='gvclass_contam_') as tmpdir:
            tmpdir = Path(tmpdir)
            with tarfile.open(tar_path, 'r:gz') as tf:
                targets = [m for m in tf.getmembers() if '/query_fna/' in m.name or '/query_faa/' in m.name or '/blastp_out/' in m.name]
                tf.extractall(path=tmpdir, members=targets)
            sample_root = next(tmpdir.iterdir())
            query_fna_dir = sample_root / 'query_fna'
            query_faa_dir = sample_root / 'query_faa'
            blast_dir = sample_root / 'blastp_out'
            query_fna = next(query_fna_dir.glob('*.fna')) if query_fna_dir.exists() else Path()
            query_faa = next(query_faa_dir.glob('*.faa')) if query_faa_dir.exists() else Path()
            features = scorer.collect_contig_features(query_fna, query_faa, blast_dir, primary_order, primary_family)
            outputs[sample_id] = {
                'suspicious_bp_fraction_v2': features['suspicious_bp_fraction'],
                'suspicious_contig_count_v2': features['suspicious_contig_count'],
            }
    return outputs


def choose_best_model(X: pd.DataFrame, y: pd.Series, groups: pd.Series) -> Any:
    candidates = [
        ('random_forest', RandomForestRegressor(n_estimators=400, random_state=RANDOM_SEED, min_samples_leaf=2)),
        ('extra_trees', ExtraTreesRegressor(n_estimators=500, random_state=RANDOM_SEED, min_samples_leaf=2)),
        ('hist_gbm', HistGradientBoostingRegressor(random_state=RANDOM_SEED, max_depth=6, min_samples_leaf=5)),
    ]
    group_kfold = GroupKFold(n_splits=min(5, len(set(groups))))
    best_name = ''
    best_model = None
    best_mae = None
    for name, model in candidates:
        maes = []
        for train_idx, valid_idx in group_kfold.split(X, y, groups):
            model.fit(X.iloc[train_idx], y.iloc[train_idx])
            pred = model.predict(X.iloc[valid_idx])
            maes.append(mean_absolute_error(y.iloc[valid_idx], pred))
        score = mean(maes)
        if best_mae is None or score < best_mae:
            best_mae = score
            best_name = name
            best_model = model
    best_model.fit(X, y)
    return best_name, best_model


def optimize_threshold(y_true: List[int], scores: List[float]) -> float:
    best_threshold = 0.0
    best_f1 = -1.0
    for threshold in [x / 2 for x in range(2, 101)]:
        preds = [1 if score >= threshold else 0 for score in scores]
        tp = sum(1 for yt, yp in zip(y_true, preds) if yt == 1 and yp == 1)
        fp = sum(1 for yt, yp in zip(y_true, preds) if yt == 0 and yp == 1)
        fn = sum(1 for yt, yp in zip(y_true, preds) if yt == 1 and yp == 0)
        if tp == 0:
            f1 = 0.0
        else:
            precision = tp / (tp + fp) if (tp + fp) else 0.0
            recall = tp / (tp + fn) if (tp + fn) else 0.0
            f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
        if f1 > best_f1:
            best_f1 = f1
            best_threshold = threshold
    return best_threshold


def evaluate_predictions(df: pd.DataFrame, score_col: str, threshold: float) -> Dict[str, Any]:
    y_true = df['contamination_present'].astype(int).tolist()
    scores = df[score_col].astype(float).tolist()
    preds = [1 if score >= threshold else 0 for score in scores]
    tp = sum(1 for yt, yp in zip(y_true, preds) if yt == 1 and yp == 1)
    fp = sum(1 for yt, yp in zip(y_true, preds) if yt == 0 and yp == 1)
    fn = sum(1 for yt, yp in zip(y_true, preds) if yt == 1 and yp == 0)
    tn = sum(1 for yt, yp in zip(y_true, preds) if yt == 0 and yp == 0)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
    auroc = roc_auc_score(y_true, scores) if len(set(y_true)) > 1 else 0.0
    aupr = average_precision_score(y_true, scores) if len(set(y_true)) > 1 else 0.0
    mae_all = mean_absolute_error(df['true_contamination_pct'], df[score_col])
    contaminated = df[df['contamination_present'] == 1]
    mae_contaminated = mean_absolute_error(contaminated['true_contamination_pct'], contaminated[score_col]) if not contaminated.empty else 0.0
    hgt = df[df['contamination_type'] == 'hgt_like']
    hgt_false_positive_rate = float((hgt[score_col] >= threshold).mean()) if not hgt.empty else 0.0
    return {
        'threshold': round(threshold, 2),
        'auroc': round(float(auroc), 4),
        'aupr': round(float(aupr), 4),
        'precision': round(float(precision), 4),
        'recall': round(float(recall), 4),
        'f1': round(float(f1), 4),
        'mae_all': round(float(mae_all), 4),
        'mae_contaminated': round(float(mae_contaminated), 4),
        'hgt_false_positive_rate': round(float(hgt_false_positive_rate), 4),
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'tn': tn,
    }


def build_feature_table(manifest_path: Path, results_dir: Path, database: Path) -> pd.DataFrame:
    manifest_df = pd.read_csv(manifest_path, sep='\t')
    summary_df = pd.read_csv(results_dir / 'gvclass_summary.tsv', sep='\t')
    summary_rows = {row['query']: row.to_dict() for _, row in summary_df.iterrows()}
    method2 = compute_method2_features(results_dir, summary_rows, database)

    rows = []
    for _, manifest_row in manifest_df.iterrows():
        sample_id = manifest_row['sample_id']
        summary = summary_rows.get(sample_id)
        if summary is None:
            continue
        row = dict(manifest_row)
        row.update(method2.get(sample_id, {'suspicious_bp_fraction_v2': 0.0, 'suspicious_contig_count_v2': 0}))
        for key in ['contamination_score_v1','contamination_cellular_signal_v1','contamination_phage_signal_v1','contamination_duplication_signal_v1','contamination_viral_mixture_signal_v1','contamination_nonviral_hit_fraction_v1','order_dup','gvog8_dup','cellular_unique','cellular_total','phage_unique','phage_total','contigs','estimated_completeness']:
            row[key] = parse_float(summary.get(key, 0.0))
        row['contamination_present'] = 1 if row['true_contamination_pct'] > 0 else 0
        rows.append(row)
    return pd.DataFrame(rows)


def train_and_benchmark(feature_df: pd.DataFrame, run_dir: Path) -> Dict[str, Any]:
    train_df = feature_df[feature_df['split'] == 'train'].copy()
    test_df = feature_df[feature_df['split'] == 'test'].copy()

    method1_threshold = optimize_threshold(train_df['contamination_present'].astype(int).tolist(), train_df['contamination_score_v1'].astype(float).tolist())
    method2_threshold = optimize_threshold(train_df['contamination_present'].astype(int).tolist(), train_df['suspicious_bp_fraction_v2'].astype(float).tolist())

    X_train = train_df[METHOD3_FEATURES].fillna(0.0)
    y_train = train_df['true_contamination_pct'].astype(float)
    groups = train_df['base_accession']
    model_name, model = choose_best_model(X_train, y_train, groups)

    train_df['method3_score'] = model.predict(X_train)
    test_df['method3_score'] = model.predict(test_df[METHOD3_FEATURES].fillna(0.0))
    method3_threshold = optimize_threshold(train_df['contamination_present'].astype(int).tolist(), train_df['method3_score'].astype(float).tolist())

    metrics_rows = []
    for method_name, score_col, threshold in [
        ('method1_rule_based', 'contamination_score_v1', method1_threshold),
        ('method2_contig_aware', 'suspicious_bp_fraction_v2', method2_threshold),
        (f'method3_{model_name}', 'method3_score', method3_threshold),
    ]:
        metrics = evaluate_predictions(test_df, score_col, threshold)
        metrics['method'] = method_name
        metrics_rows.append(metrics)

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(run_dir / 'benchmark_metrics.tsv', sep='\t', index=False)
    feature_df_out = pd.concat([train_df, test_df], ignore_index=True)
    feature_df_out.to_csv(run_dir / 'benchmark_feature_table.tsv', sep='\t', index=False)

    report_lines = [
        '# Contamination Benchmark',
        '',
        '## Design',
        '- Base genomes: one-contig isolate references from 7 families (train/test split by genome within family).',
        '- Scenarios: clean, HGT-like cellular gene insertions, cellular contig contamination, related viral contig contamination, divergent viral contig contamination.',
        '- Truth label: added contaminant contig fraction (%); HGT-like scenarios labeled as 0 contamination.',
        '',
        '## Test Metrics',
        '',
        metrics_df.to_markdown(index=False),
        '',
        '## Scenario Means (test split)',
        '',
        test_df.groupby('scenario')[['contamination_score_v1','suspicious_bp_fraction_v2','method3_score','true_contamination_pct']].mean().round(2).to_markdown(),
        '',
        '## HGT False Positives',
        '',
        test_df[test_df['contamination_type'] == 'hgt_like'][['sample_id','contamination_score_v1','suspicious_bp_fraction_v2','method3_score']].round(2).to_markdown(index=False),
    ]
    (run_dir / 'benchmark_report.md').write_text('\n'.join(report_lines) + '\n')

    model_meta = {
        'model_name': model_name,
        'method1_threshold': method1_threshold,
        'method2_threshold': method2_threshold,
        'method3_threshold': method3_threshold,
    }
    (run_dir / 'method3_model_meta.json').write_text(json.dumps(model_meta, indent=2) + '\n')
    return {
        'metrics_df': metrics_df,
        'test_df': test_df,
        'model_name': model_name,
        'thresholds': model_meta,
    }


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    input_fna_dir = Path(args.input_fna_dir).resolve()
    metadata_rows = load_metadata(Path(args.metadata).resolve())
    manifest_path = generate_dataset(run_dir, input_fna_dir, metadata_rows, args.force_regenerate)
    run_gvclass(run_dir / 'inputs', run_dir / 'results', args.workers, args.force_rerun)
    feature_df = build_feature_table(manifest_path, run_dir / 'results', Path(args.database).resolve())
    outcome = train_and_benchmark(feature_df, run_dir)
    print(f"benchmark_rows\t{len(feature_df)}")
    print(f"method3_model\t{outcome['model_name']}")
    for _, row in outcome['metrics_df'].iterrows():
        print(f"{row['method']}\tAUROC={row['auroc']}\tAUPR={row['aupr']}\tMAE={row['mae_all']}\tF1={row['f1']}")


if __name__ == '__main__':
    main()
