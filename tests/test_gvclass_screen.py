"""The ``gvclass screen`` command runs both contig screens on every bin of a folder."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.bin import gvclass_cli, gvclass_screen as screen


class Recorder:
    """Stands in for the two screens: records their arguments, writes a table."""

    def __init__(self, fail_on: str | None = None) -> None:
        self.calls: list[tuple[str, list[str]]] = []
        self.fail_on = fail_on

    def cellular(self, argv: list[str]) -> int:
        self.calls.append(("cellular", argv))
        out = Path(argv[argv.index("--out") + 1])
        pd.DataFrame({"contig_id": ["a"], "p_cellular": [0.1]}).to_csv(
            out, sep="\t", index=False
        )
        print("called cellular: 1000 bp (10.0% of the input)")
        return 0

    def giant_virus(self, argv: list[str]) -> int:
        self.calls.append(("giant_virus", argv))
        out = Path(argv[argv.index("--out") + 1])
        pd.DataFrame({"contig": ["a"], "p_foreign": [0.9]}).to_csv(
            out, sep="\t", index=False
        )
        print("foreign giant-virus contigs: 1, 5000 bp = 5.0% of the input")
        return 0


@pytest.fixture
def patched(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    recorder = Recorder()
    monkeypatch.setattr(screen, "run_cellular", recorder.cellular)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    monkeypatch.setattr(screen, "bundle_present", lambda database: True)
    return recorder


def make_bins(folder: Path, names: list[str]) -> None:
    folder.mkdir()
    for name in names:
        (folder / name).write_text(">c1\nACGT\n", encoding="utf-8")


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

    table = pd.read_csv(out / "contig_screen_summary.tsv", sep="\t")
    assert table["bin"].tolist() == ["a", "b"]
    assert table["cellular_status"].tolist() == ["ok", "ok"]
    assert table["giant_virus_status"].tolist() == ["ok", "ok"]
    text = (out / "a.summary.txt").read_text(encoding="utf-8")
    assert "called cellular: 1000 bp" in text
    assert "foreign giant-virus contigs: 1" in text


def test_a_failing_screen_is_recorded_and_the_other_bins_still_run(
    tmp_path, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])
    recorder = Recorder()

    def broken(argv: list[str]) -> int:
        recorder.calls.append(("cellular", argv))
        return 2 if "a.fna" in argv[0] else recorder.cellular(argv)

    monkeypatch.setattr(screen, "run_cellular", broken)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    monkeypatch.setattr(screen, "bundle_present", lambda database: True)
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert code == 1
    table = pd.read_csv(out / "contig_screen_summary.tsv", sep="\t")
    assert table["cellular_status"].tolist() == ["failed", "ok"]
    assert table["giant_virus_status"].tolist() == ["ok", "ok"]
    gv_calls = [a for m, a in recorder.calls if m == "giant_virus"]
    assert "--cellular" not in gv_calls[0]
    assert "--cellular" in gv_calls[1]


def test_a_folder_without_bins_is_an_error(tmp_path, patched, capsys):
    (tmp_path / "empty").mkdir()

    code = screen.main(
        [str(tmp_path / "empty"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    assert code == 1
    assert "no FASTA" in capsys.readouterr().err
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


def test_a_missing_screen_bundle_is_reported_before_any_work(
    tmp_path, monkeypatch, capsys
):
    recorder = Recorder()
    monkeypatch.setattr(screen, "run_cellular", recorder.cellular)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    make_bins(tmp_path / "bins", ["a.fna"])

    code = screen.main(
        [str(tmp_path / "bins"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    assert code == 1
    assert "gvclass screen --install" in capsys.readouterr().err
    assert recorder.calls == []


def test_the_bundle_is_present_when_its_manifest_is(tmp_path):
    assert not screen.bundle_present(tmp_path)
    (tmp_path / "contig_screen").mkdir()
    (tmp_path / "contig_screen" / "manifest.json").write_text("{}")

    assert screen.bundle_present(tmp_path)


def test_install_downloads_the_addon_into_the_database(tmp_path, monkeypatch):
    seen = []

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

    code = screen.main(["--install", "-d", str(tmp_path)])

    assert code == 0
    assert seen == [(tmp_path, "contig_screen", "v2.1.0")]


def test_install_without_a_configured_source_is_an_error(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(screen.DatabaseManager, "CONTIG_SCREEN_SOURCE", None)

    code = screen.main(["--install", "-d", str(tmp_path)])

    assert code == 1
    assert "no download source" in capsys.readouterr().err


def test_an_exception_in_one_screen_fails_that_bin_and_the_batch_goes_on(
    tmp_path, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna", "b.fna"])
    recorder = Recorder()

    def explodes(argv: list[str]) -> int:
        if "a.fna" in argv[0]:
            msg = "no contig of 3 kb"
            raise ValueError(msg)
        return recorder.cellular(argv)

    monkeypatch.setattr(screen, "run_cellular", explodes)
    monkeypatch.setattr(screen, "run_giant_virus", recorder.giant_virus)
    monkeypatch.setattr(screen, "bundle_present", lambda database: True)
    out = tmp_path / "out"

    code = screen.main([str(tmp_path / "bins"), "-o", str(out), "-d", str(tmp_path)])

    assert code == 1
    table = pd.read_csv(out / "contig_screen_summary.tsv", sep="\t")
    assert table["cellular_status"].tolist() == ["failed", "ok"]
    assert "no contig of 3 kb" in (out / "a.summary.txt").read_text(encoding="utf-8")


def test_bins_that_share_a_name_are_refused_before_any_work(tmp_path, patched, capsys):
    make_bins(tmp_path / "bins", ["a.fa", "a.fna", "b.fna"])

    code = screen.main(
        [str(tmp_path / "bins"), "-o", str(tmp_path / "out"), "-d", str(tmp_path)]
    )

    assert code == 1
    error = capsys.readouterr().err
    assert "same output name" in error
    assert "a.fa" in error
    assert patched.calls == []


def test_tables_of_an_earlier_run_do_not_survive_a_run_without_a_result(
    tmp_path, monkeypatch
):
    make_bins(tmp_path / "bins", ["a.fna"])
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.giant_virus.tsv").write_text("stale\n", encoding="utf-8")
    (out / "a.cellular.tsv").write_text("stale\n", encoding="utf-8")
    recorder = Recorder()
    monkeypatch.setattr(screen, "run_cellular", recorder.cellular)
    monkeypatch.setattr(screen, "run_giant_virus", lambda argv: 0)
    monkeypatch.setattr(screen, "bundle_present", lambda database: True)

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
