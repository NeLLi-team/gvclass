"""The ``gvclass screen`` command runs both contig screens on every bin of a folder."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from src.bin import gvclass_apptainer, gvclass_cli, gvclass_screen as screen
from src.contig_screen import extract_genes, layout
from tests.contig_screen.conftest import write_bundle

STAGE2_CUT = 0.5
DIGEST = "ab" * 32


class Recorder:
    """Stands in for the two screens: records their arguments, writes a table."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, list[str]]] = []
        self.windows = [("c1", 0, 4000, 0.9)]
        self.contigs = [("c1", 6000, True)]
        self.writes_giant_table = True

    def cellular(self, argv: list[str]) -> int:
        self.calls.append(("cellular", argv))
        out = Path(argv[argv.index("--out") + 1])
        pd.DataFrame(
            self.windows,
            columns=["source_contig", "source_start", "length_bp", "p_context"],
        ).to_csv(out, sep="\t", index=False)
        print("called cellular: 1000 bp (10.0% of the input)")
        return 0

    def giant_virus(self, argv: list[str]) -> int:
        self.calls.append(("giant_virus", argv))
        out = Path(argv[argv.index("--out") + 1])
        if self.writes_giant_table:
            pd.DataFrame(
                self.contigs, columns=["contig", "length_bp", "called_foreign"]
            ).to_csv(out, sep="\t", index=False)
        print("foreign giant-virus contigs: 1, 5000 bp = 5.0% of the input")
        return 0


@pytest.fixture
def patched(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    recorder = Recorder()
    monkeypatch.setattr(screen, "run_cellular", recorder.cellular)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    monkeypatch.setattr(screen, "bundle_present", lambda database: True)
    monkeypatch.setattr(
        screen, "load_addon", lambda database: screen.Addon(DIGEST, lambda: STAGE2_CUT)
    )
    return recorder


def write_bin(path: Path, lengths: dict[str, int]) -> None:
    """Write a bin with one contig of each length."""
    records = "".join(
        f">{name}\n{'ACGT' * (bp // 4)}\n" for name, bp in lengths.items()
    )
    path.write_text(records, encoding="utf-8")


def make_bins(folder: Path, names: list[str]) -> None:
    folder.mkdir()
    for name in names:
        write_bin(folder / name, {"c1": 4000})


def read_summary(out: Path) -> pd.DataFrame:
    return pd.read_csv(out / "contig_screen_summary.tsv", sep="\t")


def test_every_bin_of_a_folder_is_screened_once_by_each_model(tmp_path, patched):
    make_bins(tmp_path / "bins", ["b.fna", "a.fna", "notes.txt"])
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert code == 0
    assert [(m, Path(a[0]).name) for m, a in patched.calls] == [
        ("cellular", "a.fna"),
        ("giant_virus", "a.fna"),
        ("cellular", "b.fna"),
        ("giant_virus", "b.fna"),
    ]
    assert sorted(p.name for p in out.glob("*.tsv")) == [
        "a.cellular.tsv",
        "a.giant_virus.tsv",
        "b.cellular.tsv",
        "b.giant_virus.tsv",
        "contig_screen_summary.tsv",
    ]


def test_the_giant_virus_screen_gets_the_cellular_table_for_the_union(
    tmp_path, patched
):
    make_bins(tmp_path / "bins", ["a.fna"])
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    argv = patched.calls[1][1]
    assert argv[argv.index("--cellular") + 1] == str(out / "a.cellular.tsv")
    assert argv[argv.index("--resources") + 1] == str(tmp_path)


def test_a_single_file_is_screened_and_the_threads_are_passed_on(tmp_path, patched):
    make_bins(tmp_path / "bins", ["a.fna"])

    screen.main(
        [
            str(tmp_path / "bins" / "a.fna"),
            "-o",
            str(tmp_path / "out"),
            "-d",
            str(tmp_path),
            "-t",
            "7",
        ]
    )

    assert len(patched.calls) == 2
    for _, argv in patched.calls:
        assert argv[argv.index("--threads") + 1] == "7"


def test_the_summary_has_one_row_per_bin_with_the_printed_numbers(tmp_path, patched):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    table = read_summary(out)
    assert table["bin"].tolist() == ["a", "b"]
    assert table["cellular_status"].tolist() == ["ok", "ok"]
    assert table["giant_virus_status"].tolist() == ["ok", "ok"]
    text = (out / "a.summary.txt").read_text(encoding="utf-8")
    assert "called cellular: 1000 bp" in text
    assert "foreign giant-virus contigs: 1" in text


def test_the_row_has_the_base_pairs_and_percentages_of_both_screens(tmp_path, patched):
    write_bin_folder = tmp_path / "bins"
    write_bin_folder.mkdir()
    write_bin(write_bin_folder / "a.fna", {"c1": 6000, "c2": 6000, "c3": 4000})
    patched.windows = [("c3", 0, 4000, 0.9)]
    patched.contigs = [("c1", 6000, False), ("c2", 6000, True)]
    out = tmp_path / "out"

    screen.main([str(write_bin_folder), "-o", str(out), "-d", str(tmp_path)])

    row = read_summary(out).iloc[0]
    assert row["input_bp"] == 16000
    assert (row["cellular_bp"], row["cellular_pct"]) == (4000, 25.0)
    assert (row["giant_virus_bp"], row["giant_virus_pct"]) == (6000, 37.5)
    assert (row["flagged_bp"], row["flagged_pct"]) == (10000, 62.5)


def test_the_cellular_cut_of_the_bundle_decides_which_base_pairs_are_called(
    tmp_path, patched, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna"])
    monkeypatch.setattr(
        screen, "load_addon", lambda database: screen.Addon(DIGEST, lambda: 0.95)
    )

    screen.main(
        [str(tmp_path / "bins"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    row = read_summary(tmp_path / "out").iloc[0]
    assert row["cellular_bp"] == 0


def test_a_bin_without_a_contig_of_3_kb_is_not_scored_by_the_cellular_screen(
    tmp_path, patched
):
    folder = tmp_path / "bins"
    folder.mkdir()
    write_bin(folder / "a.fna", {"c1": 2000, "c2": 2800})
    out = tmp_path / "out"

    code = screen.main([str(folder), "-o", str(out), "-d", str(tmp_path)])

    row = read_summary(out).iloc[0]
    assert code == 0
    assert row["cellular_status"] == "not_scored"
    assert row["giant_virus_status"] == "ok"
    assert (
        row[["cellular_bp", "cellular_pct", "flagged_bp", "flagged_pct"]].isna().all()
    )
    assert row["input_bp"] == 4800
    assert [model for model, _ in patched.calls] == ["giant_virus"]
    assert "--cellular" not in patched.calls[0][1]


def test_a_bin_with_one_scored_contig_has_no_giant_virus_numbers(tmp_path, patched):
    make_bins(tmp_path / "bins", ["a.fna"])
    patched.writes_giant_table = False
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    row = read_summary(out).iloc[0]
    assert code == 0
    assert row["giant_virus_status"] == "not_scored"
    assert row["cellular_status"] == "ok"
    assert (row["cellular_bp"], row["cellular_pct"]) == (4000, 100.0)
    assert row[["giant_virus_bp", "giant_virus_pct"]].isna().all()
    assert row[["flagged_bp", "flagged_pct"]].isna().all()


def test_a_failing_screen_is_recorded_and_the_other_bins_still_run(
    tmp_path, monkeypatch, patched
):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])

    def broken(argv: list[str]) -> int:
        patched.calls.append(("cellular", argv))
        return 2 if "a.fna" in argv[0] else patched.cellular(argv)

    monkeypatch.setattr(screen, "run_cellular", broken)
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert code == 1
    table = read_summary(out)
    assert table["cellular_status"].tolist() == ["failed", "ok"]
    assert table["giant_virus_status"].tolist() == ["ok", "ok"]
    assert pd.isna(table.loc[0, "cellular_pct"])
    assert pd.isna(table.loc[0, "flagged_pct"])
    gv_calls = [a for m, a in patched.calls if m == "giant_virus"]
    assert "--cellular" not in gv_calls[0]
    assert "--cellular" in gv_calls[1]


def test_a_folder_without_bins_is_an_error(tmp_path, patched, caplog):
    (tmp_path / "empty").mkdir()

    code = screen.main(
        [str(tmp_path / "empty"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    assert code == 1
    assert "no FASTA" in caplog.text
    assert patched.calls == []


def test_the_database_comes_from_the_option_then_the_environment(
    tmp_path, patched, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna"])
    other = tmp_path / "other"
    other.mkdir()
    monkeypatch.setenv("GVCLASS_DB", str(other))

    screen.main([str(tmp_path / "bins"), "-o", str(tmp_path / "out")])
    screen.main(
        [str(tmp_path / "bins"), "-o", str(tmp_path / "out2"), "-d", str(tmp_path)]
    )

    first, _, second, _ = patched.calls
    assert first[1][first[1].index("--resources") + 1] == str(other)
    assert second[1][second[1].index("--resources") + 1] == str(tmp_path)


def test_the_launcher_hands_the_screen_command_to_the_screen(monkeypatch):
    seen: list[list[str]] = []
    monkeypatch.setattr(screen, "main", lambda argv: seen.append(argv) or 0)
    monkeypatch.setattr("sys.argv", ["gvclass", "screen", "bins", "-o", "out"])

    code = gvclass_cli.main()

    assert code == 0
    assert seen == [["bins", "-o", "out"]]


@pytest.mark.parametrize(
    "flags", [["--apptainer"], ["--image", "x.sif"], ["--image=x.sif"]]
)
def test_the_screen_command_goes_to_the_container_launcher_with_a_container_flag(
    monkeypatch, flags
):
    seen: list[list[str]] = []
    monkeypatch.setattr(gvclass_apptainer, "main", lambda argv: seen.append(argv) or 0)
    monkeypatch.setattr(screen, "main", lambda argv: pytest.fail("ran natively"))
    monkeypatch.setattr("sys.argv", ["gvclass", "screen", "--install", *flags])

    code = gvclass_cli.main()

    assert code == 0
    assert seen == [["screen", "--install", *flags]]


def test_a_missing_screen_bundle_is_reported_before_any_work(
    tmp_path, monkeypatch, caplog
):
    recorder = Recorder()
    monkeypatch.setattr(screen, "run_cellular", recorder.cellular)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    make_bins(tmp_path / "bins", ["a.fna"])

    code = screen.main(
        [str(tmp_path / "bins"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    assert code == 1
    assert "gvclass screen --install" in caplog.text
    assert recorder.calls == []


def test_the_bundle_is_present_when_its_manifest_is(tmp_path, monkeypatch):
    assert not screen.bundle_present(tmp_path)
    (tmp_path / "contig_screen").mkdir()
    (tmp_path / "contig_screen" / "manifest.json").write_text("{}")

    assert screen.bundle_present(tmp_path)


def test_the_addon_gives_its_manifest_digest_and_its_cellular_cut(
    tmp_path, monkeypatch
):
    digest = write_bundle(
        tmp_path / "db", {layout.THRESHOLDS_FILE: b'{"stage1": 0.4, "stage2": 0.3}'}
    )
    monkeypatch.setattr(layout, "MANIFEST_SHA256", digest)

    addon = screen.load_addon(tmp_path / "db")

    assert addon.manifest_sha256 == digest
    assert addon.stage2() == 0.3


def test_a_bundle_that_fails_its_integrity_check_fails_every_bin_and_runs_no_screen(
    tmp_path, patched, monkeypatch, caplog
):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])

    def refuse(database):
        raise layout.IntegrityError("manifest SHA-256 mismatch")

    monkeypatch.setattr(screen, "load_addon", refuse)
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    table = read_summary(out)
    assert code == 1
    assert "mismatch" in caplog.text
    assert patched.calls == []
    assert table["bin"].tolist() == ["a", "b"]
    assert table["cellular_status"].tolist() == ["failed", "failed"]
    assert table["giant_virus_status"].tolist() == ["failed", "failed"]
    assert table[["input_bp", "cellular_pct", "flagged_pct"]].isna().all().all()
    assert "mismatch" in (out / "a.summary.txt").read_text(encoding="utf-8")


def test_a_rejected_bundle_fails_a_bin_too_short_to_score_with_a_failed_row(
    tmp_path, patched, monkeypatch
):
    """The status says why nothing was scored: the add-on, not the bin."""
    make_bins(tmp_path / "bins", ["tiny.fna"])

    def refuse(database):
        raise layout.IntegrityError("manifest SHA-256 mismatch")

    monkeypatch.setattr(screen, "load_addon", refuse)
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    first_line, row = (out / "contig_screen_summary.tsv").read_text().splitlines()
    assert first_line.startswith("bin\tcellular_status\tgiant_virus_status\t")
    assert row.startswith("tiny\tfailed\tfailed\t")


def fake_installer(monkeypatch: pytest.MonkeyPatch) -> list[tuple[Path, str, str]]:
    seen: list[tuple[Path, str, str]] = []

    def fake_install(db_path, name, source):
        seen.append((Path(db_path), name, source["version"]))
        return Path(db_path) / name

    monkeypatch.setattr(screen.DatabaseManager, "install_addon", fake_install)
    monkeypatch.setattr(
        screen.DatabaseManager,
        "CONTIG_SCREEN_SOURCE",
        {
            "version": "v2.1.0",
            "url": "https://example.com/x.tar.gz",
            "sha256": "a" * 64,
        },
    )
    return seen


def test_install_downloads_the_addon_into_the_database(tmp_path, monkeypatch):
    seen = fake_installer(monkeypatch)

    code = screen.main(["--install", "-d", str(tmp_path)])

    assert code == 0
    assert seen == [(tmp_path, "contig_screen", "v2.1.0")]


def test_install_creates_a_database_directory_that_does_not_exist_yet(
    tmp_path, monkeypatch
):
    seen = fake_installer(monkeypatch)

    code = screen.main(["--install", "-d", str(tmp_path / "new" / "db")])

    assert code == 0
    assert seen == [(tmp_path / "new" / "db", "contig_screen", "v2.1.0")]
    assert (tmp_path / "new" / "db").is_dir()


def test_install_without_a_configured_source_is_an_error(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(screen.DatabaseManager, "CONTIG_SCREEN_SOURCE", None)

    code = screen.main(["--install", "-d", str(tmp_path)])

    assert code == 1
    assert "no download source" in caplog.text


def test_an_exception_in_one_screen_fails_that_bin_and_the_batch_goes_on(
    tmp_path, monkeypatch, patched
):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])

    def explodes(argv: list[str]) -> int:
        if "a.fna" in argv[0]:
            msg = "unexpected failure"
            raise ValueError(msg)
        return patched.cellular(argv)

    monkeypatch.setattr(screen, "run_cellular", explodes)
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert code == 1
    assert read_summary(out)["cellular_status"].tolist() == ["failed", "ok"]
    assert "unexpected failure" in (out / "a.summary.txt").read_text(encoding="utf-8")


def test_bins_that_share_a_name_are_refused_before_any_work(tmp_path, patched, caplog):
    make_bins(tmp_path / "bins", ["a.fa", "a.fna", "b.fna"])

    code = screen.main(
        [str(tmp_path / "bins"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    assert code == 1
    assert "same output name" in caplog.text
    assert "a.fa" in caplog.text
    assert patched.calls == []


def test_tables_of_an_earlier_run_do_not_survive_a_run_without_a_result(
    tmp_path, monkeypatch, patched
):
    make_bins(tmp_path / "bins", ["a.fna"])
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.giant_virus.tsv").write_text("stale\n", encoding="utf-8")
    (out / "a.cellular.tsv").write_text("stale\n", encoding="utf-8")
    monkeypatch.setattr(screen, "run_giant_virus", lambda argv: 0)
    patched.writes_giant_table = False

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert not (out / "a.giant_virus.tsv").exists()
    assert "stale" not in (out / "a.cellular.tsv").read_text(encoding="utf-8")


def test_the_current_directory_can_be_screened_without_an_output_option(
    tmp_path, patched, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna"])
    monkeypatch.chdir(tmp_path / "bins")

    code = screen.main([".", "-d", str(tmp_path)])

    assert code == 0
    assert (tmp_path / "bins_screen" / "contig_screen_summary.tsv").exists()


def test_gzipped_fas_files_are_found_in_a_folder(tmp_path):
    make_bins(tmp_path / "bins", ["a.fas.gz", "b.txt"])

    assert [p.name for p in screen.list_bins(tmp_path / "bins")] == ["a.fas.gz"]
    assert screen.bin_name(Path("x/a.fas.gz")) == "a"


def test_gzipped_files_with_an_upper_case_suffix_are_listed_and_readable(tmp_path):
    folder = tmp_path / "bins"
    folder.mkdir()
    with gzip.open(folder / "a.FASTA.GZ", "wt", encoding="utf-8") as handle:
        handle.write(">c1\nACGT\n")

    listed = screen.list_bins(folder)

    assert [path.name for path in listed] == ["a.FASTA.GZ"]
    assert screen.bin_name(listed[0]) == "a"
    with extract_genes.open_text(listed[0]) as handle:
        assert handle.read() == ">c1\nACGT\n"


def test_a_cellular_screen_that_exits_cleanly_without_a_table_is_not_scored(
    tmp_path, patched, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna"])
    monkeypatch.setattr(screen, "run_cellular", lambda argv: 0)
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    row = read_summary(out).iloc[0]
    assert code == 0
    assert row["cellular_status"] == "not_scored"
    assert row["giant_virus_status"] == "ok"
    assert "--cellular" not in patched.calls[0][1]


def test_the_summary_table_shows_missing_numbers_as_nan(tmp_path, patched):
    make_bins(tmp_path / "bins", ["a.fna"])
    patched.writes_giant_table = False
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    header, row = (out / "contig_screen_summary.tsv").read_text().splitlines()
    fields = dict(zip(header.split("\t"), row.split("\t"), strict=True))
    assert fields["giant_virus_pct"] == "NaN"
    assert fields["cellular_pct"] == "100.00"
    assert fields["input_bp"] == "4000"


def test_the_row_carries_the_digest_of_the_manifest_that_the_screens_used(
    tmp_path, patched
):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert read_summary(out)["manifest_sha256"].tolist() == [DIGEST, DIGEST]


@pytest.mark.parametrize(
    ("name", "other"),
    [("1", "2"), ("007", "008"), ("NA", "other"), ("NA", "NaN")],
)
def test_contig_names_that_look_like_numbers_or_na_are_kept_as_text(
    tmp_path, name, other
):
    cellular_table = tmp_path / "cellular.tsv"
    pd.DataFrame(
        [(name, 0, 6000, 0.9)],
        columns=["source_contig", "source_start", "length_bp", "p_context"],
    ).to_csv(cellular_table, sep="\t", index=False)
    giant_table = tmp_path / "giant.tsv"
    pd.DataFrame(
        [(name, 6000, True), (other, 6000, False)],
        columns=["contig", "length_bp", "called_foreign"],
    ).to_csv(giant_table, sep="\t", index=False)

    totals = screen.bin_totals(
        {name: 6000, other: 6000}, cellular_table, giant_table, 0.5
    )

    assert totals["cellular_bp"] == 6000
    assert totals["giant_virus_bp"] == 6000
    assert totals["flagged_bp"] == 6000
    assert totals["cellular_pct"] == 50.0


def test_a_contig_named_na_is_counted_in_the_row_of_the_bin(tmp_path, patched):
    folder = tmp_path / "bins"
    folder.mkdir()
    write_bin(folder / "a.fna", {"NA": 6000, "007": 6000})
    patched.windows = [("NA", 0, 6000, 0.9)]
    patched.contigs = [("007", 6000, True), ("NA", 6000, False)]

    screen.main([str(folder), "-o", str(tmp_path / "out"), "-d", str(tmp_path)])

    row = read_summary(tmp_path / "out").iloc[0]
    assert (row["cellular_bp"], row["giant_virus_bp"], row["flagged_bp"]) == (
        6000,
        6000,
        12000,
    )


def test_the_row_carries_the_digest_of_the_input_file(tmp_path, patched):
    make_bins(tmp_path / "bins", ["a.fna"])
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    expected = hashlib.sha256((tmp_path / "bins" / "a.fna").read_bytes()).hexdigest()
    assert read_summary(out)["input_sha256"].tolist() == [expected]


def test_bins_that_cannot_be_read_get_failed_rows_and_the_batch_goes_on(
    tmp_path, patched
):
    folder = tmp_path / "bins"
    folder.mkdir()
    write_bin(folder / "a.fna", {"c1": 4000})
    (folder / "b.fna").write_text(
        ">c1\n" + "ACGT" * 1000 + "\n>c1\nACGT\n", encoding="utf-8"
    )
    whole = gzip.compress((">c1\n" + "ACGT" * 1000 + "\n").encode())
    (folder / "c.fna.gz").write_bytes(whole[: len(whole) // 2])
    write_bin(folder / "d.fna", {"c1": 4000})
    out = tmp_path / "out"

    code = screen.main([str(folder), "-o", str(out), "-d", str(tmp_path)])

    table = read_summary(out)
    assert code == 1
    assert table["bin"].tolist() == ["a", "b", "c", "d"]
    assert table["cellular_status"].tolist() == ["ok", "failed", "failed", "ok"]
    assert table["giant_virus_status"].tolist() == ["ok", "failed", "failed", "ok"]
    unreadable = table.loc[table["bin"].isin(["b", "c"])]
    assert unreadable[["input_bp", "flagged_pct", "cellular_pct"]].isna().all().all()
    assert [Path(argv[0]).name for _, argv in patched.calls] == [
        "a.fna",
        "a.fna",
        "d.fna",
        "d.fna",
    ]
    assert "error" in (out / "b.summary.txt").read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "thresholds",
    [None, b'{"stage1": 0.4}', b"not json"],
    ids=["unlisted", "missing-stage2", "corrupt"],
)
def test_broken_cellular_thresholds_leave_the_giant_virus_result(
    tmp_path, monkeypatch, thresholds
):
    files = {layout.HMM_FILE: b"hmm"}
    if thresholds is not None:
        files[layout.THRESHOLDS_FILE] = thresholds
    digest = write_bundle(tmp_path / "db", files)
    monkeypatch.setattr(layout, "MANIFEST_SHA256", digest)
    recorder = Recorder()

    def cellular_like_the_real_one(argv: list[str]) -> int:
        recorder.calls.append(("cellular", argv))
        resources = layout.ScreenResources.open(
            Path(argv[argv.index("--resources") + 1])
        )
        layout.load_thresholds(resources)
        return 0

    monkeypatch.setattr(screen, "run_cellular", cellular_like_the_real_one)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    folder = tmp_path / "bins"
    folder.mkdir()
    write_bin(folder / "a.fna", {"c1": 6000, "c2": 6000})
    out = tmp_path / "out"

    code = screen.main([str(folder), "-o", str(out), "-d", str(tmp_path / "db")])

    row = read_summary(out).iloc[0]
    assert code == 1
    assert (row["cellular_status"], row["giant_virus_status"]) == ("failed", "ok")
    assert row["giant_virus_pct"] == 50.0
    assert row[["cellular_pct", "flagged_pct"]].isna().all()
    assert "--cellular" not in recorder.calls[1][1]


def test_the_input_digest_is_taken_before_the_screens_read_the_bin(
    tmp_path, patched, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna"])
    path = tmp_path / "bins" / "a.fna"
    original = hashlib.sha256(path.read_bytes()).hexdigest()
    replaced = {}

    def replaces_the_bin(argv):
        write_bin(path, {"c1": 8000})
        replaced["digest"] = hashlib.sha256(path.read_bytes()).hexdigest()
        return patched.cellular(argv)

    monkeypatch.setattr(screen, "run_cellular", replaces_the_bin)
    out = tmp_path / "out"

    screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert replaced["digest"] != original
    assert read_summary(out)["input_sha256"].tolist() == [original]


def test_duplicate_contig_names_fail_the_bin_with_a_clear_error(tmp_path, patched):
    folder = tmp_path / "bins"
    folder.mkdir()
    (folder / "a.fna").write_text(
        f">c\n{'ACGT' * 4500}\n>c\n{'ACGT' * 1500}\n", encoding="utf-8"
    )
    write_bin(folder / "b.fna", {"c1": 4000})
    out = tmp_path / "out"

    code = screen.main([str(folder), "-o", str(out), "-d", str(tmp_path)])

    table = read_summary(out)
    assert code == 1
    assert table["cellular_status"].tolist() == ["failed", "ok"]
    assert table.loc[0, ["input_bp", "cellular_pct", "flagged_pct"]].isna().all()
    message = (out / "a.summary.txt").read_text(encoding="utf-8")
    assert "duplicate" in message
    assert "'c'" in message
    assert [Path(argv[0]).name for _, argv in patched.calls] == ["b.fna", "b.fna"]


def test_the_summary_is_published_whole_or_not_at_all(tmp_path, patched, monkeypatch):
    make_bins(tmp_path / "bins", ["a.fna"])
    out = tmp_path / "out"
    out.mkdir()
    (out / "contig_screen_summary.tsv").write_text("earlier run\n", encoding="utf-8")

    def disk_full(self, path_or_buf, *args, **kwargs):
        Path(path_or_buf).write_text("bin\tcellular_st", encoding="utf-8")
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(pd.DataFrame, "to_csv", disk_full)

    with pytest.raises(OSError, match="No space"):
        screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert (out / "contig_screen_summary.tsv").read_text() == "earlier run\n"
    assert sorted(path.name for path in out.glob("*.part")) == []
