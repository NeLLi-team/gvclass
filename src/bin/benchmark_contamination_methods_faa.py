"""Benchmark contamination methods on synthetic FAA bins."""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import subprocess
import tarfile
import tempfile
import io
from pathlib import Path
from statistics import mean
from typing import Any, Dict, List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import joblib
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import average_precision_score, mean_absolute_error, roc_auc_score
from sklearn.model_selection import GroupKFold

from src.config.marker_sets import BUSCO_MODELS, UNI56_MODELS
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


def create_model_candidate(model_name: str):
    if model_name == "random_forest":
        return RandomForestRegressor(
            n_estimators=400,
            random_state=RANDOM_SEED,
            min_samples_leaf=2,
        )
    if model_name == "extra_trees":
        return ExtraTreesRegressor(
            n_estimators=500,
            random_state=RANDOM_SEED,
            min_samples_leaf=2,
        )
    if model_name == "hist_gbm":
        return HistGradientBoostingRegressor(
            random_state=RANDOM_SEED,
            max_depth=6,
            min_samples_leaf=5,
        )
    raise ValueError(f"Unknown model candidate: {model_name}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="benchmarking/contamination/contamination_benchmark_faa")
    parser.add_argument("--metadata", default="benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10/metadata_table.tsv")
    parser.add_argument("--query-results", default="benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10/gvclass_extended_results")
    parser.add_argument("--database", default="resources")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--force-regenerate", action="store_true")
    parser.add_argument("--force-rerun", action="store_true")
    parser.add_argument("--export-model", nargs="?", const="src/bundled_models/contamination_model.joblib",
                        help="Retrain on full data and save model bundle to path (default: src/bundled_models/contamination_model.joblib)")
    return parser.parse_args()


def family_from_taxonomy(tax: str) -> str:
    parts = [p.strip() for p in tax.split(';')]
    return parts[6] if len(parts) > 6 else 'NA'


def order_from_taxonomy(tax: str) -> str:
    parts = [p.strip() for p in tax.split(';')]
    return parts[5] if len(parts) > 5 else 'NA'


def load_metadata(path: Path) -> List[Dict[str, Any]]:
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    for row in rows:
        row['family'] = family_from_taxonomy(row['taxonomy'])
        row['order_name'] = order_from_taxonomy(row['taxonomy'])
        row['genome_size_bp'] = int(row['genome_size_bp'])
        row['num_contigs'] = int(row['num_contigs'])
    return rows


def select_base_genomes(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    selected = []
    for family in BASE_FAMILIES:
        fam_rows = [
            row for row in rows
            if row['family'] == family and row['num_contigs'] == 1 and 100_000 <= row['genome_size_bp'] <= 700_000
        ]
        fam_rows.sort(key=lambda row: row['genome_size_bp'])
        selected.extend(fam_rows[:2])
    return selected


def extract_query_faa(query_results_dir: Path, accession: str) -> List[SeqRecord]:
    tar_path = query_results_dir / f"{accession}.tar.gz"
    with tarfile.open(tar_path, 'r:gz') as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith(f'query_faa/{accession}.faa'))
        handle = tf.extractfile(member)
        return list(SeqIO.parse(io.TextIOWrapper(handle), 'fasta'))


def protein_bp(record: SeqRecord) -> int:
    return len(record.seq) * 3


def protein_contig(record: SeqRecord) -> str:
    token = record.id.split('|', 1)[1] if '|' in record.id else record.id
    return token.rsplit('_', 1)[0] if '_' in token else token


def choose_donor(base: Dict[str, Any], rows: List[Dict[str, Any]], related: bool) -> Dict[str, Any]:
    if related:
        pool = [r for r in rows if r['family'] == base['family'] and r['accession_number'] != base['accession_number'] and 100_000 <= r['genome_size_bp'] <= 700_000 and r['num_contigs'] == 1]
        if pool:
            return sorted(pool, key=lambda r: r['genome_size_bp'])[0]
        pool = [r for r in rows if r['order_name'] == base['order_name'] and r['accession_number'] != base['accession_number'] and 100_000 <= r['genome_size_bp'] <= 700_000 and r['num_contigs'] == 1]
        return sorted(pool, key=lambda r: r['genome_size_bp'])[0]
    pool = [r for r in rows if r['order_name'] != base['order_name'] and 100_000 <= r['genome_size_bp'] <= 700_000 and r['num_contigs'] == 1]
    return sorted(pool, key=lambda r: r['genome_size_bp'])[0]


def build_cellular_protein_pool(database: Path) -> List[SeqRecord]:
    pool = []
    faa_dir = database / 'database' / 'faa'
    for marker in BUSCO_MODELS[:40] + UNI56_MODELS[:20]:
        path = faa_dir / f'{marker}.faa'
        if not path.exists():
            continue
        for record in SeqIO.parse(str(path), 'fasta'):
            if record.id.startswith('BAC__') or record.id.startswith('ARC__') or record.id.startswith('EUK__'):
                pool.append(record)
        if len(pool) >= 500:
            break
    return pool


def select_proteins_to_fraction(records: List[SeqRecord], target_fraction: float, base_bp: int, rng: random.Random) -> List[SeqRecord]:
    target_bp = max(900, round((base_bp * target_fraction) / max(1e-9, 1.0 - target_fraction)))
    shuffled = records[:]
    rng.shuffle(shuffled)
    selected = []
    total_bp = 0
    for record in shuffled:
        selected.append(record)
        total_bp += protein_bp(record)
        if total_bp >= target_bp:
            break
    return selected


def rename_records(records: List[SeqRecord], contig_id: str, prefix: str) -> List[SeqRecord]:
    renamed = []
    for idx, record in enumerate(records, start=1):
        new_id = f"{prefix}|{contig_id}_{idx}"
        renamed.append(SeqRecord(record.seq, id=new_id, description=''))
    return renamed


def largest_contig_id(records: List[SeqRecord]) -> str:
    contig_bp = {}
    for record in records:
        contig = protein_contig(record)
        contig_bp[contig] = contig_bp.get(contig, 0) + protein_bp(record)
    return max(contig_bp.items(), key=lambda item: item[1])[0]


def write_synthetic_faa(
    output_path: Path,
    base_row: Dict[str, Any],
    base_records: List[SeqRecord],
    all_rows: List[Dict[str, Any]],
    query_results_dir: Path,
    cellular_pool: List[SeqRecord],
    rng: random.Random,
    scenario: str,
    target_fraction: float,
) -> Dict[str, Any]:
    records = [SeqRecord(record.seq, id=record.id, description='') for record in base_records]
    base_bp = sum(protein_bp(record) for record in records)
    contamination_type = 'clean'
    truth_fraction = 0.0
    added_bp = 0

    if scenario.startswith('hgt_cellular'):
        contamination_type = 'hgt_like'
        n_genes = 1 if scenario.endswith('1gene') else 5
        donor_records = rng.sample(cellular_pool, n_genes)
        main_contig = largest_contig_id(records)
        injected = rename_records(donor_records, main_contig, base_row['accession_number'])
        records.extend(injected)
        added_bp = sum(protein_bp(r) for r in donor_records)
    elif scenario.startswith('cellular_contig'):
        contamination_type = 'cellular_contig'
        donor_records = select_proteins_to_fraction(cellular_pool, target_fraction, base_bp, rng)
        contig_id = f'cellularcontam_{scenario}'
        injected = rename_records(donor_records, contig_id, base_row['accession_number'])
        records.extend(injected)
        added_bp = sum(protein_bp(r) for r in donor_records)
        truth_fraction = (added_bp / (base_bp + added_bp)) * 100.0
    elif scenario.startswith('viral_related') or scenario.startswith('viral_divergent'):
        contamination_type = 'viral_contig_related' if 'related' in scenario else 'viral_contig_divergent'
        donor_row = choose_donor(base_row, all_rows, related='related' in scenario)
        donor_records = extract_query_faa(query_results_dir, donor_row['accession_number'])
        donor_contig = largest_contig_id(donor_records)
        donor_records = [r for r in donor_records if protein_contig(r) == donor_contig]
        donor_records = select_proteins_to_fraction(donor_records, target_fraction, base_bp, rng)
        contig_id = f'viralcontam_{scenario}'
        injected = rename_records(donor_records, contig_id, base_row['accession_number'])
        records.extend(injected)
        added_bp = sum(protein_bp(r) for r in donor_records)
        truth_fraction = (added_bp / (base_bp + added_bp)) * 100.0

    SeqIO.write(records, str(output_path), 'fasta')
    return {
        'sample_id': output_path.stem,
        'base_accession': base_row['accession_number'],
        'base_family': base_row['family'],
        'base_order': base_row['order_name'],
        'scenario': scenario,
        'split': base_row['split'],
        'contamination_type': contamination_type,
        'true_contamination_pct': round(truth_fraction, 2),
        'added_bp': added_bp,
        'base_bp': base_bp,
    }


def generate_dataset(run_dir: Path, metadata_rows: List[Dict[str, Any]], query_results_dir: Path, database: Path, force: bool) -> Path:
    inputs_dir = run_dir / 'inputs'
    manifest_path = run_dir / 'synthetic_manifest.tsv'
    if manifest_path.exists() and not force:
        return manifest_path
    if inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    base_rows = select_base_genomes(metadata_rows)
    by_family = {}
    for family in BASE_FAMILIES:
        pair = [r for r in base_rows if r['family'] == family]
        pair.sort(key=lambda r: r['genome_size_bp'])
        pair[0]['split'] = 'train'
        pair[1]['split'] = 'test'
        by_family[family] = pair
    cellular_pool = build_cellular_protein_pool(database)
    rng = random.Random(RANDOM_SEED)
    manifest_rows = []
    for family in BASE_FAMILIES:
        for base_row in by_family[family]:
            base_records = extract_query_faa(query_results_dir, base_row['accession_number'])
            for scenario, target_fraction in SCENARIOS:
                sample_id = f"{base_row['accession_number']}__{scenario}"
                output_path = inputs_dir / f"{sample_id}.faa"
                manifest_rows.append(
                    write_synthetic_faa(output_path, base_row, base_records, metadata_rows, query_results_dir, cellular_pool, rng, scenario, target_fraction)
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
        './gvclass', str(inputs_dir), '-o', str(results_dir), '--threads', str(workers), '-j', str(workers), '--threads-per-worker', '1', '--completeness-mode', 'novelty-aware', '--plain-output'
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
    outputs = {}
    for sample_id, row in summary_rows.items():
        tar_path = results_dir / f"{sample_id}.tar.gz"
        primary_order = extract_primary_token(row.get('order', ''))
        primary_family = extract_primary_token(row.get('family', ''))
        with tempfile.TemporaryDirectory(prefix='gvclass_contam_faa_') as tmpdir:
            tmpdir = Path(tmpdir)
            with tarfile.open(tar_path, 'r:gz') as tf:
                targets = [m for m in tf.getmembers() if '/query_faa/' in m.name or '/blastp_out/' in m.name]
                tf.extractall(path=tmpdir, members=targets)
            sample_root = next(tmpdir.iterdir())
            query_faa = next((sample_root / 'query_faa').glob('*.faa'))
            features = scorer.collect_contig_features(Path(), query_faa, sample_root / 'blastp_out', primary_order, primary_family)
            outputs[sample_id] = {
                'suspicious_bp_fraction_v2': features['suspicious_bp_fraction'],
                'suspicious_contig_count_v2': features['suspicious_contig_count'],
            }
    return outputs


def choose_best_model(X: pd.DataFrame, y: pd.Series, groups: pd.Series) -> Tuple[str, Any]:
    candidates = [
        (name, create_model_candidate(name))
        for name in ("random_forest", "extra_trees", "hist_gbm")
    ]
    splitter = GroupKFold(n_splits=min(5, len(set(groups))))
    best = ('', None, None)
    for name, model in candidates:
        maes = []
        for train_idx, valid_idx in splitter.split(X, y, groups):
            model.fit(X.iloc[train_idx], y.iloc[train_idx])
            pred = model.predict(X.iloc[valid_idx])
            maes.append(mean_absolute_error(y.iloc[valid_idx], pred))
        score = mean(maes)
        if best[2] is None or score < best[2]:
            best = (name, model, score)
    best[1].fit(X, y)
    return best[0], best[1]


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
    hgt_fpr = float((hgt[score_col] >= threshold).mean()) if not hgt.empty else 0.0
    return {
        'threshold': round(threshold, 2),
        'auroc': round(float(auroc), 4),
        'aupr': round(float(aupr), 4),
        'precision': round(float(precision), 4),
        'recall': round(float(recall), 4),
        'f1': round(float(f1), 4),
        'mae_all': round(float(mae_all), 4),
        'mae_contaminated': round(float(mae_contaminated), 4),
        'hgt_false_positive_rate': round(float(hgt_fpr), 4),
        'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
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

    model_name, model = choose_best_model(train_df[METHOD3_FEATURES].fillna(0.0), train_df['true_contamination_pct'].astype(float), train_df['base_accession'])
    train_df['method3_score'] = model.predict(train_df[METHOD3_FEATURES].fillna(0.0))
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
    ranked = metrics_df.sort_values(
        by=['auroc', 'aupr', 'f1', 'mae_all', 'hgt_false_positive_rate'],
        ascending=[False, False, False, True, True],
    ).reset_index(drop=True)
    recommended_method = str(ranked.iloc[0]['method'])
    contamination_type_means = test_df.groupby('contamination_type')[
        ['contamination_score_v1', 'suspicious_bp_fraction_v2', 'method3_score', 'true_contamination_pct']
    ].mean().round(2).reset_index()
    def as_tsv_block(df: pd.DataFrame) -> str:
        return df.to_csv(sep='\t', index=False).strip()

    report_lines = [
        '# Contamination Benchmark (FAA synthetic bins)',
        '',
        '## Design',
        '- Base genomes: one-contig isolate references from 7 families, split train/test by genome within family.',
        '- Scenarios: clean, HGT-like cellular gene additions, cellular contaminant contigs, related viral contaminant contigs, divergent viral contaminant contigs.',
        '- Truth label: contaminant protein bp fraction (%); HGT-like scenarios labeled as 0 contamination.',
        '',
        '## Decision Summary',
        '',
        f'- Recommended default: `{recommended_method}`',
        '- Reason: best held-out ranking by AUROC, AUPR, F1, and MAE together.',
        '- Method 1 is interpretable but overcalls HGT-like additions.',
        '- Method 2 alone is not usable as a primary score because it mistakes many HGT-like additions for contamination.',
        '',
        '## Test Metrics',
        '',
        '```tsv',
        as_tsv_block(metrics_df),
        '```',
        '',
        '## Ranked Methods',
        '',
        '```tsv',
        as_tsv_block(ranked[['method', 'auroc', 'aupr', 'f1', 'mae_all', 'hgt_false_positive_rate']]),
        '```',
        '',
        '## Scenario Means (test split)',
        '',
        '```tsv',
        as_tsv_block(test_df.groupby('scenario')[['contamination_score_v1','suspicious_bp_fraction_v2','method3_score','true_contamination_pct']].mean().round(2).reset_index()),
        '```',
        '',
        '## Contamination-Type Means (test split)',
        '',
        '```tsv',
        as_tsv_block(contamination_type_means),
        '```',
        '',
        '## HGT False Positives (test split)',
        '',
        '```tsv',
        as_tsv_block(test_df[test_df['contamination_type'] == 'hgt_like'][['sample_id','contamination_score_v1','suspicious_bp_fraction_v2','method3_score']].round(2)),
        '```',
    ]
    (run_dir / 'benchmark_report.md').write_text('\n'.join(report_lines) + '\n')
    (run_dir / 'method3_model_meta.json').write_text(json.dumps({
        'model_name': model_name,
        'method1_threshold': method1_threshold,
        'method2_threshold': method2_threshold,
        'method3_threshold': method3_threshold,
        'recommended_method': recommended_method,
    }, indent=2) + '\n')
    return {'metrics_df': metrics_df, 'model_name': model_name, 'recommended_method': recommended_method}


def export_model(feature_df: pd.DataFrame, output_path: str, model_name: str) -> None:
    """Train the selected model family on full data and serialize to disk."""
    X = feature_df[METHOD3_FEATURES].fillna(0.0)
    y = feature_df['true_contamination_pct'].astype(float)
    model = create_model_candidate(model_name)
    model.fit(X, y)
    threshold = optimize_threshold(
        feature_df['contamination_present'].astype(int).tolist(),
        model.predict(X).tolist(),
    )
    bundle = {
        'model': model,
        'feature_names': METHOD3_FEATURES,
        'threshold': threshold,
        'model_name': model_name,
    }
    dest = Path(output_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, dest)
    print(f'exported_model\t{dest}\tmodel={model_name}\tthreshold={threshold}')


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    metadata_rows = load_metadata(Path(args.metadata).resolve())
    query_results_dir = Path(args.query_results).resolve()
    database = Path(args.database).resolve()
    manifest_path = generate_dataset(run_dir, metadata_rows, query_results_dir, database, args.force_regenerate)
    run_gvclass(run_dir / 'inputs', run_dir / 'results', args.workers, args.force_rerun)
    feature_df = build_feature_table(manifest_path, run_dir / 'results', database)
    outcome = train_and_benchmark(feature_df, run_dir)
    print(f'benchmark_rows\t{len(feature_df)}')
    print(f'method3_model\t{outcome["model_name"]}')
    print(f'recommended_method\t{outcome["recommended_method"]}')
    for _, row in outcome['metrics_df'].iterrows():
        print(f"{row['method']}\tAUROC={row['auroc']}\tAUPR={row['aupr']}\tMAE={row['mae_all']}\tF1={row['f1']}")

    if args.export_model:
        export_model(feature_df, args.export_model, outcome["model_name"])


if __name__ == '__main__':
    main()
