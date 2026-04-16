from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.hmm_search import (
    _collect_filtered_search_results,
    generate_counts_from_output,
    process_hmmout,
)
from src.core.hmm_search_multi import (
    _dedup_hits,
    dedup_domtbl_file,
)


def _domtbl_line(
    protein_id: str,
    model_name: str,
    seq_score: float,
    dom_score: float,
    target_length: int = 100,
) -> str:
    return (
        f"{protein_id}\t-\t{target_length}\t{model_name}\t-\t300\t1e-20\t{seq_score}\t0.0\t"
        f"1\t1\t1e-20\t1e-20\t{dom_score}\t0.0\t1\t100\t1\t100\t1\t100\t0.90\n"
    )


def test_generate_counts_from_output_deduplicates_same_protein_and_model(
    tmp_path: Path,
) -> None:
    output_file = tmp_path / "models.out.filtered"
    counts_file = tmp_path / "models.counts"
    output_file.write_text(
        "# header\n"
        + _domtbl_line("prot1", "GVOGm0022", 50.0, 50.0)
        + _domtbl_line("prot1", "GVOGm0022", 60.0, 60.0)
        + _domtbl_line("prot2", "GVOGm0022", 55.0, 55.0)
        + _domtbl_line("prot1", "GVOGm0054", 42.0, 42.0)
    )

    generate_counts_from_output(str(output_file), str(counts_file))

    counts = dict(
        line.strip().split("\t")
        for line in counts_file.read_text().splitlines()
        if line.strip()
    )
    assert counts["GVOGm0022"] == "2"
    assert counts["GVOGm0054"] == "1"


def test_collect_filtered_search_results_deduplicates_same_protein_and_model(
    tmp_path: Path,
) -> None:
    output_file = tmp_path / "models.out"
    output_file.write_text(
        "# header\n"
        + _domtbl_line("prot1", "GVOGm0022", 50.0, 50.0)
        + _domtbl_line("prot1", "GVOGm0022", 60.0, 60.0)
        + _domtbl_line("prot2", "GVOGm0022", 55.0, 55.0)
        + _domtbl_line("prot1", "GVOGm0054", 42.0, 42.0)
    )

    model_hits, filtered_lines, total_lines, hit_lines, passed_cutoffs = (
        _collect_filtered_search_results(str(output_file), cutoffs={}, has_ga_cutoffs=False)
    )

    assert total_lines == 5
    assert hit_lines == 4
    assert passed_cutoffs == 3
    assert sorted(model_hits["GVOGm0022"]) == [("prot1", 60.0), ("prot2", 55.0)]
    assert model_hits["GVOGm0054"] == [("prot1", 42.0)]
    assert len([line for line in filtered_lines if not line.startswith("#")]) == 3


def test_process_hmmout_deduplicates_same_protein_and_model_for_og_markers(
    tmp_path: Path,
) -> None:
    hmmout = tmp_path / "models.out"
    filtered_hmmout = tmp_path / "models.out.filtered"
    hmmout.write_text(
        "# header\n"
        + _domtbl_line("genome1|prot1", "OG0001", 50.0, 50.0)
        + _domtbl_line("genome1|prot1", "OG0001", 70.0, 70.0)
        + _domtbl_line("genome1|prot2", "OG0001", 55.0, 55.0)
    )

    count_df, _ = process_hmmout(
        models=["OG0001"],
        hmmout=str(hmmout),
        filtered_hmmout=str(filtered_hmmout),
        cutoffs={"OG0001": (40.0, 0.0)},
    )

    assert int(count_df.loc["genome1", "OG0001"]) == 2
    non_header_lines = [
        line
        for line in filtered_hmmout.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(non_header_lines) == 2
    assert any("\t70.0\t" in line for line in non_header_lines)


# ---------------------------------------------------------------------------
# Phase 1.1 regression tests: multi-HMM dedup in hmm_search_multi.
# ---------------------------------------------------------------------------


def test_dedup_hits_collapses_duplicate_query_model_pairs() -> None:
    """The pure in-memory helper must keep exactly one row per (target, model)."""
    rows = [
        ("prot1", "GVOGm0022", 50.0, 50.0, _domtbl_line("prot1", "GVOGm0022", 50.0, 50.0)),
        # Second domain of the same protein-model pair with higher scores.
        ("prot1", "GVOGm0022", 60.0, 60.0, _domtbl_line("prot1", "GVOGm0022", 60.0, 60.0)),
        ("prot2", "GVOGm0022", 55.0, 55.0, _domtbl_line("prot2", "GVOGm0022", 55.0, 55.0)),
        ("prot1", "GVOGm0054", 42.0, 42.0, _domtbl_line("prot1", "GVOGm0054", 42.0, 42.0)),
    ]

    lines = _dedup_hits(rows)

    assert len(lines) == 3
    # The 60.0 row wins over the 50.0 row for (prot1, GVOGm0022).
    gvog0022_lines = [line for line in lines if "\tGVOGm0022\t" in line]
    assert len(gvog0022_lines) == 2
    assert any("\t60.0\t" in line for line in gvog0022_lines)
    assert not any("\t50.0\t" in line for line in gvog0022_lines)


def test_dedup_domtbl_file_matches_single_hmm_filter_semantics(tmp_path: Path) -> None:
    """``dedup_domtbl_file`` must produce the same shape as the single-HMM
    ``_collect_filtered_search_results`` dedup pass, fixing the bug where the
    multi-HMM engine silently copied its raw per-domain output into
    ``models.out.filtered``."""
    raw_file = tmp_path / "models.out"
    filtered_file = tmp_path / "models.out.filtered"
    raw_file.write_text(
        "# target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value\tscore\tbias\n"
        + _domtbl_line("prot1", "GVOGm0022", 50.0, 50.0)
        # Duplicate (prot1, GVOGm0022) via a second per-domain row — the old
        # shutil.copy2 path would have leaked this into downstream consumers.
        + _domtbl_line("prot1", "GVOGm0022", 60.0, 60.0)
        + _domtbl_line("prot2", "GVOGm0022", 55.0, 55.0)
        + _domtbl_line("prot1", "GVOGm0054", 42.0, 42.0)
    )

    written = dedup_domtbl_file(str(raw_file), str(filtered_file))

    assert written == 3
    contents = filtered_file.read_text()
    data_lines = [line for line in contents.splitlines() if line and not line.startswith("#")]
    assert len(data_lines) == 3
    # The higher-scoring (60.0) domain of prot1/GVOGm0022 survives.
    assert any("prot1\t-\t100\tGVOGm0022\t-\t300\t1e-20\t60.0" in line for line in data_lines)
    assert not any("prot1\t-\t100\tGVOGm0022\t-\t300\t1e-20\t50.0" in line for line in data_lines)
    # Header is preserved verbatim.
    assert contents.splitlines()[0].startswith("#")


def test_dedup_domtbl_file_applies_no_score_cutoffs(tmp_path: Path) -> None:
    """Dedup must never drop a hit based on score magnitude — that guards
    sensitive-mode correctness, where GA/TC/NC cutoffs are intentionally
    disabled upstream."""
    raw_file = tmp_path / "models.out"
    filtered_file = tmp_path / "models.out.filtered"
    # Deliberately emit a below-GA-threshold hit and a mid-range duplicate.
    raw_file.write_text(
        "# header\n"
        + _domtbl_line("prot_low_score", "Mirus_MCP", 3.5, 3.5)
        + _domtbl_line("prot_low_score", "Mirus_MCP", 4.0, 4.0)
    )

    dedup_domtbl_file(str(raw_file), str(filtered_file))

    data_lines = [
        line
        for line in filtered_file.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(data_lines) == 1  # dedup kept the best (4.0), did not filter it out
    assert "\t4.0\t" in data_lines[0]


def test_dedup_domtbl_file_preserves_single_rows(tmp_path: Path) -> None:
    """Non-duplicated input must pass through unchanged."""
    raw_file = tmp_path / "models.out"
    filtered_file = tmp_path / "models.out.filtered"
    raw_file.write_text(
        "# header\n"
        + _domtbl_line("a", "M1", 10.0, 10.0)
        + _domtbl_line("b", "M1", 20.0, 20.0)
        + _domtbl_line("a", "M2", 30.0, 30.0)
    )

    written = dedup_domtbl_file(str(raw_file), str(filtered_file))

    assert written == 3
    data_lines = [
        line
        for line in filtered_file.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(data_lines) == 3


def test_dedup_domtbl_file_rejects_truncated_rows(tmp_path: Path) -> None:
    """Rows with fewer columns than :func:`_format_domain_line` emits must be
    rejected so ``models.out.filtered`` never carries structurally invalid
    lines (Codex-audit finding)."""
    raw_file = tmp_path / "models.out"
    filtered_file = tmp_path / "models.out.filtered"
    # Well-formed 22-column row (from _domtbl_line) + a truncated 15-column
    # row that the previous ``len(parts) >= 15`` check would have preserved.
    truncated_row = (
        "prot_bad\t-\t100\tGVOGm0099\t-\t300\t1e-20\t42.0\t0.0\t"
        "1\t1\t1e-20\t1e-20\t41.0\t0.0\n"
    )
    raw_file.write_text(
        "# header\n"
        + _domtbl_line("prot_good", "GVOGm0022", 50.0, 50.0)
        + truncated_row
    )

    written = dedup_domtbl_file(str(raw_file), str(filtered_file))

    data_lines = [
        line
        for line in filtered_file.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert written == 1
    assert len(data_lines) == 1
    assert "prot_good" in data_lines[0]
    assert "prot_bad" not in filtered_file.read_text()
