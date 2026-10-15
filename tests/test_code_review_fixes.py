"""Regression tests for bugs identified and fixed during the comprehensive code review."""

from __future__ import annotations

import tarfile
from pathlib import Path

import pandas as pd
import pytest

from src.bin.gvclass_cli import (
    BUILTIN_CONFIG_DEFAULTS,
    build_parser,
    resolve_sensitive_mode,
)
from src.core.contamination_scoring import ContaminationScorer
from src.core.hmm_search_multi import _has_model_cutoffs, _load_hmms, clean_hmm_file
from src.core.reformat import calculate_stats
from src.core.species_tree.protein_gathering import (
    select_query_representative_from_tree,
)
from src.core.tree_analysis import TreeAnalyzer, _is_query_leaf
from src.core.weighted_completeness import WeightedCompletenessCalculator
from src.pipeline.summary_writer import (
    EXTENDED_SUMMARY_COLUMNS,
    FINAL_SUMMARY_COLUMNS,
    INDIVIDUAL_FINAL_SUMMARY_COLUMNS,
    archive_final_summary_extended_files,
    read_individual_summary_results,
    write_final_summary_extended_files,
    write_final_summary_files,
    write_individual_final_summary_file,
)
from src.utils.common import validate_file_path
from src.utils.contig_splitter import split_contigs_from_directory
from src.utils.input_validation import InputValidator


def test_is_query_leaf_exact_and_delimited_matching() -> None:
    """_is_query_leaf should match exact, pipe-delimited, or single-underscore suffixes,
    without matching unrelated reference leaves that happen to share a prefix."""
    assert _is_query_leaf("G1", "G1") is True
    assert _is_query_leaf("G1|c1_1", "G1") is True
    assert _is_query_leaf("G1_1", "G1") is True
    assert _is_query_leaf("G10|c1_1", "G1") is False
    assert _is_query_leaf("G10_1", "G1") is False
    assert _is_query_leaf("NCLDV__G1_1", "NCLDV") is False
    assert _is_query_leaf("NCLDV__QUERY_1", "NCLDV__QUERY") is True


def test_find_nearest_neighbors_and_query_reps_ignore_prefix_colliding_refs(
    tmp_path: Path,
) -> None:
    """Query 'G1' must not match reference 'G10|c_1' as a query leaf."""
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "G10\td_Vir|p_NCLDV|c_Megaviricetes|o_Imitervirales|f_Mimiviridae|g_Mimi|s_Mimi_sp\n"
    )
    tree_file = tmp_path / "OG1.treefile"
    tree_file.write_text("((G1|q_1:0.1,G10|c_1:0.2):0.1,G10|c_2:0.5);\n")

    analyzer = TreeAnalyzer(labels)
    neighbors = analyzer.get_neighbors_from_tree(tree_file, query_id="G1")
    assert list(neighbors.keys()) == ["G1|q_1"]
    assert list(neighbors["G1|q_1"].keys()) == ["G10|c_1"]
    assert abs(neighbors["G1|q_1"]["G10|c_1"] - 0.3) < 1e-6

    rep = select_query_representative_from_tree(
        tree_file,
        "G1",
        keep_prefix="G10",
        query_lengths={"G1|q_1": 5, "G10|c_1": 12},
    )
    assert rep is not None
    assert rep[0] == "G1|q_1"


def test_taxonomy_consensus_and_marker_trees_use_7_rank_indices(
    tmp_path: Path,
) -> None:
    """7-rank labels (domain|phylum|class|order|family|genus|species) must map
    indices 3..6 to order, family, genus, species."""
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "NCLDV__G1\td_Viruses|p_Nucleocytoviricota|c_Megaviricetes|"
        "o_Imitervirales|f_Mimiviridae|g_Mimivirus|s_Mimivirus_bradfordmassiliense\n"
    )
    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    tree_file = tree_dir / "GVOGm0003.treefile"
    tree_file.write_text("(QUERY|c1_1:0.1,NCLDV__G1|p1:0.2);\n")

    analyzer = TreeAnalyzer(labels)
    output_tsv = tmp_path / "nearest_neighbors.tsv"
    all_neighbors = analyzer.process_marker_trees(tree_dir, "QUERY", output_tsv)
    assert "GVOGm0003" in all_neighbors

    assert analyzer.get_taxonomy_consensus(all_neighbors, "order") == "o_Imitervirales"
    assert analyzer.get_taxonomy_consensus(all_neighbors, "family") == "f_Mimiviridae"
    assert analyzer.get_taxonomy_consensus(all_neighbors, "genus") == "g_Mimivirus"
    assert (
        analyzer.get_taxonomy_consensus(all_neighbors, "species")
        == "s_Mimivirus_bradfordmassiliense"
    )

    content = output_tsv.read_text()
    assert "d_Viruses;o_Imitervirales;f_Mimiviridae;g_Mimivirus" in content


def test_clean_hmm_file_preserves_cutoffs_and_has_model_cutoffs_checks_ga(
    tmp_path: Path,
) -> None:
    """clean_hmm_file must retain GA/TC/NC lines, and _has_model_cutoffs must check GA."""
    hmm_with_tc_only = (
        "HMMER3/f [3.3.2 | Nov 2020]\n"
        "NAME  test_tc\n"
        "LENG  1\n"
        "ALPH  amino\n"
        "TC    25.00 25.00;\n"
        "HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y\n"
        "            m->m     m->i     m->d     i->m     i->i     d->m     d->d\n"
        "  COMPO   2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.91696  2.42397  2.77519  2.98518  4.58477  3.61503\n"
        "          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.91696  2.42397  2.77519  2.98518  4.58477  3.61503\n"
        "          0.01656  4.52135  5.29832  0.61904  0.77255  0.00000     *   \n"
        "      1   2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.91696  2.42397  2.77519  2.98518  4.58477  3.61503      1 a - - -\n"
        "          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.91696  2.42397  2.77519  2.98518  4.58477  3.61503\n"
        "          0.01656  4.52135  5.29832  0.61904  0.77255  0.48551  0.95551\n"
        "//\n"
    )
    tc_file = tmp_path / "tc_only.hmm"
    tc_file.write_text(hmm_with_tc_only)
    tc_hmms = _load_hmms(str(tc_file))
    assert len(tc_hmms) == 1
    assert _has_model_cutoffs(tc_hmms[0]) is False

    ga_file = tmp_path / "with_ga.hmm"
    ga_file.write_text(hmm_with_tc_only.replace("TC    25.00 25.00;", "GA    25.00 25.00;"))
    cleaned_ga = clean_hmm_file(str(ga_file))
    cleaned_text = Path(cleaned_ga).read_text()
    assert "GA    25.00 25.00;" in cleaned_text
    ga_hmms = _load_hmms(cleaned_ga)
    assert len(ga_hmms) == 1
    assert _has_model_cutoffs(ga_hmms[0]) is True


def test_protein_to_contig_id_only_strips_numeric_cds_suffix() -> None:
    """protein_to_contig_id should strip trailing _<digits> but preserve non-numeric suffixes."""
    assert ContaminationScorer.protein_to_contig_id("query|contig_1_12") == "contig_1"
    assert ContaminationScorer.protein_to_contig_id("query|contig_A") == "contig_A"
    assert ContaminationScorer.protein_to_contig_id("contig_alpha") == "contig_alpha"


def test_weighted_completeness_fallback_averages_duplicate_marker_rows(
    tmp_path: Path,
) -> None:
    """Fallback rows of one marker in several orders weigh like their mean."""
    header = "order_name\tmarker_name\tpercent_genomes_with_marker\n"
    duplicated = tmp_path / "duplicated.tsv"
    duplicated.write_text(
        header + "OrderA\tOG1\t80.0\nOrderB\tOG1\t90.0\nOrderA\tOG2\t70.0\n"
    )
    averaged = tmp_path / "averaged.tsv"
    averaged.write_text(header + "OrderA\tOG1\t85.0\nOrderA\tOG2\t70.0\n")
    markers = ("OG1", "OG2")

    from_duplicates = WeightedCompletenessCalculator(
        duplicated, enable_adaptive_scaling=False
    ).get_marker_weights("UnknownOrder", markers)
    from_means = WeightedCompletenessCalculator(
        averaged, enable_adaptive_scaling=False
    ).get_marker_weights("UnknownOrder", markers)

    assert from_duplicates == pytest.approx(from_means)


def test_summary_writer_preserves_extended_columns_across_resume(
    tmp_path: Path,
) -> None:
    """<query>.final_summary.tsv must persist every summary and extended column.

    --resume rebuilds the extended table from these files.
    """
    import logging

    assert (
        INDIVIDUAL_FINAL_SUMMARY_COLUMNS[: len(FINAL_SUMMARY_COLUMNS)]
        == FINAL_SUMMARY_COLUMNS
    )
    assert set(EXTENDED_SUMMARY_COLUMNS) <= set(INDIVIDUAL_FINAL_SUMMARY_COLUMNS)

    summary_data = {col: "0" for col in INDIVIDUAL_FINAL_SUMMARY_COLUMNS}
    summary_data["query"] = "sample_q"
    summary_data["cellular_coherent_bp_fraction"] = 12.34
    summary_data["contig_attribution_mode"] = "contig_aware"

    write_individual_final_summary_file(
        tmp_path, "sample_q", summary_data, logging.getLogger(__name__)
    )

    resumed_results = read_individual_summary_results(tmp_path, ["sample_q"])
    write_final_summary_files(resumed_results, tmp_path)
    write_final_summary_extended_files(resumed_results, tmp_path)
    archive_final_summary_extended_files(tmp_path)

    main_df = pd.read_csv(tmp_path / "gvclass_summary.tsv", sep="\t")
    assert list(main_df.columns) == FINAL_SUMMARY_COLUMNS

    tar_path = tmp_path / "gvclass_summary.extended.tar.gz"
    assert tar_path.exists()
    with tarfile.open(tar_path, "r:gz") as tar:
        member = tar.extractfile("gvclass_summary.extended.tsv")
        assert member is not None
        ext_df = pd.read_csv(member, sep="\t")
    assert list(ext_df.columns) == EXTENDED_SUMMARY_COLUMNS
    assert float(ext_df.loc[0, "cellular_coherent_bp_fraction"]) == 12.34
    assert ext_df.loc[0, "contig_attribution_mode"] == "contig_aware"


def test_cli_no_sensitive_flag() -> None:
    """--no-sensitive should resolve sensitive_mode to False."""
    parser = build_parser()
    default_args = parser.parse_args(["example", "-o", "out"])
    assert resolve_sensitive_mode(default_args, BUILTIN_CONFIG_DEFAULTS) is True

    no_sens_args = parser.parse_args(["example", "-o", "out", "--no-sensitive"])
    assert resolve_sensitive_mode(no_sens_args, BUILTIN_CONFIG_DEFAULTS) is False


def test_fa_and_fas_extension_handling(tmp_path: Path) -> None:
    """.fa files should be classified by sequence content and .fas files split properly."""
    fa_file = tmp_path / "genome.fa"
    fa_file.write_text(">contig1\n" + "ATGC" * 6000 + "\n")
    validated = InputValidator.validate_sequence_file(fa_file)
    assert validated == fa_file.resolve()

    fas_dir = tmp_path / "fas_in"
    fas_dir.mkdir()
    (fas_dir / "multi.fas").write_text(">c1\n" + "A" * 100 + "\n>c2\n" + "G" * 100 + "\n")
    split_out = tmp_path / "fas_out"
    out_dir, n_contigs, n_files = split_contigs_from_directory(
        fas_dir, split_out, min_length=50
    )
    assert out_dir == split_out
    assert n_contigs == 2
    assert n_files == 1


def test_validate_path_safety_allows_double_dots_in_filename(tmp_path: Path) -> None:
    """Filenames containing '..' inside the stem are safe if no path component is '..'."""
    f = tmp_path / "genome..v1.fna"
    f.write_text(">c1\nATGC\n")
    resolved = validate_file_path(f)
    assert resolved == f.resolve()


def test_reformat_calculate_stats_counts_lowercase_gc() -> None:
    """calculate_stats should count both uppercase and lowercase G/C bases."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = [SeqRecord(Seq("atgc"), id="c1")]
    stats_df = calculate_stats(records, "lower.fna")
    assert float(stats_df.loc[0, "GCperc"]) == 50.0
