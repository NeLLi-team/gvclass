"""Regression tests for required search and taxonomy evidence failures."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from src.core import blast, hmm_search_multi, marker_processing
from src.core.blast import parse_blastp, run_blastp
from src.core.contamination_scoring import ContaminationScorer
from src.core.hmm_search import generate_counts_from_output
from src.core.marker_processing import MarkerProcessor
from src.core.tree_analysis import TreeAnalyzer
from src.utils.error_handling import ProcessingError


def _write_search_fastas(tmp_path: Path) -> tuple[Path, Path]:
    query = tmp_path / "query.faa"
    reference = tmp_path / "reference.faa"
    query.write_text(">query\nMKVLAA\n", encoding="utf-8")
    reference.write_text(">reference\nMKVLAA\n", encoding="utf-8")
    return query, reference


def _search_hit() -> SimpleNamespace:
    search_result = SimpleNamespace(
        identity=lambda: 1.0,
        query_start=0,
        query_end=5,
        target_start=0,
        target_end=5,
    )
    return SimpleNamespace(
        evalue=1e-20,
        query_index=0,
        target_index=0,
        result=search_result,
        query_aln="MKVLAA",
        target_aln="MKVLAA",
        score=42.0,
    )


def _search_fails_after_hit(*_args: object, **_kwargs: object):
    yield _search_hit()
    raise RuntimeError("injected failure after hit")


def _bare_contamination_scorer() -> ContaminationScorer:
    scorer = ContaminationScorer.__new__(ContaminationScorer)
    scorer.labels = {}
    return scorer


def _malformed_m8_record() -> str:
    return "query\tsubject\tbad-identity\t6\t0\t0\t1\t6\t1\t6\t1e-20\t42"


def test_blast_search_failure_before_first_hit_publishes_no_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    query, reference = _write_search_fastas(tmp_path)
    output = tmp_path / "search.m8"
    monkeypatch.setattr(
        blast.pyswrd,
        "search",
        Mock(side_effect=RuntimeError("injected failure before hit")),
    )

    with pytest.raises(ProcessingError, match="injected failure before hit"):
        run_blastp(str(query), str(reference), str(output), threads=1)

    assert not output.exists()
    assert list(tmp_path.glob(".search.m8.*.tmp")) == []


def test_blast_search_failure_after_first_hit_publishes_no_partial_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    query, reference = _write_search_fastas(tmp_path)
    output = tmp_path / "search.m8"
    monkeypatch.setattr(blast.pyswrd, "search", _search_fails_after_hit)

    with pytest.raises(ProcessingError, match="injected failure after hit"):
        run_blastp(str(query), str(reference), str(output), threads=1)

    assert not output.exists()
    assert list(tmp_path.glob(".search.m8.*.tmp")) == []


def test_blast_successful_zero_hits_publishes_empty_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    query, reference = _write_search_fastas(tmp_path)
    output = tmp_path / "search.m8"
    monkeypatch.setattr(blast.pyswrd, "search", lambda *_args, **_kwargs: iter(()))

    run_blastp(str(query), str(reference), str(output), threads=1)

    assert output.exists()
    assert output.read_bytes() == b""


def test_parse_blastp_rejects_malformed_record(tmp_path: Path) -> None:
    output = tmp_path / "malformed.m8"
    output.write_text("query\tsubject\t90\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"malformed\.m8 at line 1"):
        parse_blastp(str(output))


def test_multi_hmm_rejects_missing_requested_file_before_search(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    present_hmm = tmp_path / "present.hmm"
    present_hmm.touch()
    missing_hmm = tmp_path / "missing.hmm"
    query = tmp_path / "query.faa"
    query.write_text(">query\nMKVLAA\n", encoding="utf-8")
    load_sequences = Mock(return_value=[])
    monkeypatch.setattr(hmm_search_multi, "_load_query_sequences", load_sequences)

    with pytest.raises(FileNotFoundError, match=r"missing\.hmm"):
        hmm_search_multi.run_pyhmmer_search_multi(
            [str(present_hmm), str(missing_hmm)],
            str(query),
            str(tmp_path / "models.out"),
            threads=1,
        )

    load_sequences.assert_not_called()


def test_multi_hmm_aborts_when_selected_file_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_hmm = tmp_path / "first.hmm"
    second_hmm = tmp_path / "second.hmm"
    first_hmm.touch()
    second_hmm.touch()
    output = tmp_path / "models.out"
    monkeypatch.setattr(
        hmm_search_multi, "_load_query_sequences", Mock(return_value=[])
    )
    process_hmm = Mock(
        side_effect=[
            ("first", str(first_hmm), [object()], []),
            RuntimeError("injected second HMM failure"),
        ]
    )
    monkeypatch.setattr(hmm_search_multi, "_process_single_hmm_file", process_hmm)

    with pytest.raises(RuntimeError, match="injected second HMM failure"):
        hmm_search_multi.run_pyhmmer_search_multi(
            [str(first_hmm), str(second_hmm)],
            str(tmp_path / "query.faa"),
            str(output),
            threads=1,
        )

    assert process_hmm.call_count == 2
    assert not output.exists()


def test_generate_counts_propagates_missing_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match=r"missing\.filtered"):
        generate_counts_from_output(
            str(tmp_path / "missing.filtered"),
            str(tmp_path / "models.counts"),
        )


def test_generate_counts_rejects_malformed_record(tmp_path: Path) -> None:
    hmm_output = tmp_path / "models.out.filtered"
    counts = tmp_path / "models.counts"
    hmm_output.write_text("protein\t-\t100\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"models\.out\.filtered at line 1"):
        generate_counts_from_output(str(hmm_output), str(counts))

    assert not counts.exists()


def test_generate_counts_propagates_output_failure(tmp_path: Path) -> None:
    hmm_output = tmp_path / "models.out.filtered"
    hmm_output.write_text("protein\t-\t100\tmodel\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError):
        generate_counts_from_output(
            str(hmm_output),
            str(tmp_path / "missing" / "models.counts"),
        )


def test_contamination_empty_m8_is_valid_no_hit_evidence(tmp_path: Path) -> None:
    blast_dir = tmp_path / "blastp_out"
    blast_dir.mkdir()
    (blast_dir / "marker.m8").touch()

    assert _bare_contamination_scorer().collect_best_blast_hits(blast_dir) == {}


def test_contamination_rejects_malformed_m8_record(tmp_path: Path) -> None:
    blast_dir = tmp_path / "blastp_out"
    blast_dir.mkdir()
    (blast_dir / "marker.m8").write_text(
        _malformed_m8_record() + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"marker\.m8 at line 1"):
        _bare_contamination_scorer().collect_best_blast_hits(blast_dir)


def test_contig_features_rejects_malformed_m8_record(tmp_path: Path) -> None:
    blast_dir = tmp_path / "blastp_out"
    blast_dir.mkdir()
    (blast_dir / "marker.m8").write_text(
        _malformed_m8_record() + "\n",
        encoding="utf-8",
    )
    query_faa = tmp_path / "query.faa"
    query_faa.write_text(">query\nMKVLAA\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"marker\.m8 at line 1"):
        _bare_contamination_scorer().collect_contig_features(
            tmp_path / "missing.fna",
            query_faa,
            blast_dir,
            "Imitervirales",
            "Mimiviridae",
        )


def test_contamination_propagates_unreadable_m8_file(tmp_path: Path) -> None:
    blast_dir = tmp_path / "blastp_out"
    blast_dir.mkdir()
    (blast_dir / "unreadable.m8").mkdir()

    with pytest.raises(OSError, match="Failed to read BLAST m8 file"):
        _bare_contamination_scorer().collect_best_blast_hits(blast_dir)


def test_contamination_labels_are_required_and_well_formed(tmp_path: Path) -> None:
    labels_file = tmp_path / "labels.tsv"
    labels_file.write_text("malformed-label\n", encoding="utf-8")
    scorer = _bare_contamination_scorer()
    scorer.labels_file = labels_file

    with pytest.raises(ValueError, match=r"labels\.tsv at line 1"):
        scorer._load_labels()


def test_tree_labels_are_required_and_well_formed(tmp_path: Path) -> None:
    labels_file = tmp_path / "labels.tsv"
    labels_file.write_text("malformed-label\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"labels\.tsv at line 1"):
        TreeAnalyzer(labels_file)


def test_tree_analyzer_rejects_malformed_tree(tmp_path: Path) -> None:
    labels_file = tmp_path / "labels.tsv"
    labels_file.write_text("NCLDV__ref\tNCLDV|p|o|f|g|s\n", encoding="utf-8")
    tree_file = tmp_path / "marker.treefile"
    tree_file.write_text("(query_1:0.1,NCLDV__ref_1", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Failed to parse or search"):
        TreeAnalyzer(labels_file).get_neighbors_from_tree(tree_file, "query")


def test_tree_analyzer_propagates_neighbor_search_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    labels_file = tmp_path / "labels.tsv"
    labels_file.write_text("NCLDV__ref\tNCLDV|p|o|f|g|s\n", encoding="utf-8")
    tree_file = tmp_path / "marker.treefile"
    tree_file.write_text("(query_1:0.1,NCLDV__ref_1:0.2);\n", encoding="utf-8")
    analyzer = TreeAnalyzer(labels_file)
    monkeypatch.setattr(
        analyzer,
        "_find_nearest_neighbor",
        Mock(side_effect=RuntimeError("injected distance failure")),
    )

    with pytest.raises(RuntimeError, match="injected distance failure"):
        analyzer.get_neighbors_from_tree(tree_file, "query")


def test_marker_rerun_removes_stale_skip_before_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    processor = MarkerProcessor("marker", tmp_path / "db", tmp_path / "query")
    alignment = tmp_path / "marker.faa"
    alignment.write_text(
        ">query_1\nAAAA\n>NCLDV__ref_1\nAAAA\n",
        encoding="utf-8",
    )
    assert processor.build_tree(alignment, threads=1) is None
    skip_file = processor.tree_dir / "marker.SKIPPED"
    assert skip_file.exists()
    monkeypatch.setattr(
        marker_processing,
        "run_veryfasttree",
        lambda _path, _threads: "(query_1:0.1,NCLDV__ref_1:0.2,other_1:0.3);",
    )
    alignment.write_text(
        ">query_1\nAAAA\n>NCLDV__ref_1\nAAAA\n>other_1\nAAAA\n",
        encoding="utf-8",
    )

    tree_file = processor.build_tree(alignment, threads=1)

    assert tree_file is not None
    assert tree_file.exists()
    assert not skip_file.exists()


def test_marker_rerun_replaces_stale_state_with_current_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    processor = MarkerProcessor("marker", tmp_path / "db", tmp_path / "query")
    tree_file = processor.tree_dir / "marker.treefile"
    skip_file = processor.tree_dir / "marker.SKIPPED"
    error_file = processor.tree_dir / "marker.ERROR"
    timeout_file = processor.tree_dir / "marker.TIMEOUT"
    tree_file.write_text("stale tree", encoding="utf-8")
    skip_file.write_text("stale skip", encoding="utf-8")
    error_file.write_text("stale error", encoding="utf-8")
    timeout_file.write_text("stale timeout", encoding="utf-8")
    alignment = tmp_path / "marker.faa"
    alignment.write_text(
        ">query_1\nAAAA\n>NCLDV__ref_1\nAAAA\n>other_1\nAAAA\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        marker_processing,
        "run_veryfasttree",
        Mock(side_effect=RuntimeError("current tree failure")),
    )

    with pytest.raises(RuntimeError, match="current tree failure"):
        processor.build_tree(alignment, threads=1)

    assert not tree_file.exists()
    assert not skip_file.exists()
    assert not timeout_file.exists()
    assert "current tree failure" in error_file.read_text(encoding="utf-8")
    assert "stale error" not in error_file.read_text(encoding="utf-8")


def test_tree_analysis_explicit_paths_exclude_unselected_stale_tree(
    tmp_path: Path,
) -> None:
    labels_file = tmp_path / "labels.tsv"
    labels_file.write_text("NCLDV__ref\tNCLDV|p|o|f|g|s\n", encoding="utf-8")
    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    current_tree = tree_dir / "current.treefile"
    stale_tree = tree_dir / "stale.treefile"
    tree_text = "(query_1:0.1,NCLDV__ref_1:0.2,other_1:0.3);\n"
    current_tree.write_text(tree_text, encoding="utf-8")
    stale_tree.write_text(tree_text, encoding="utf-8")
    analyzer = TreeAnalyzer(labels_file)

    selected = analyzer.process_marker_trees(
        tree_dir,
        "query",
        tmp_path / "selected.tree_nn",
        tree_paths=[current_tree],
    )
    compatibility_scan = analyzer.process_marker_trees(
        tree_dir,
        "query",
        tmp_path / "scanned.tree_nn",
    )

    assert set(selected) == {"current"}
    assert set(compatibility_scan) == {"current", "stale"}
