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
