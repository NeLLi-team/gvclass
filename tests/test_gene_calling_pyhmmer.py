"""The marker search: reported hits only, GA and E-value models, one profile parse."""

from __future__ import annotations

import io
import os
import time
from pathlib import Path
from types import SimpleNamespace

import pyhmmer
import pytest

from src.core import hmm_search
from src.core.hmm_search import (
    COMBINED_HMM,
    digitize_proteins,
    load_profiles,
    read_proteins,
    run_marker_search,
    search,
)
from tests.conftest import RESOURCES_ROOT

EXAMPLE_PROTEINS = Path(__file__).parents[1] / "example" / "PkV-RF01.faa"

PROTEINS = {
    "alpha": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLQDKPEAIRRLVDELAQLGVE",
    "beta": "MSDNLWPERAGHHTCYQVFKILMNSTREDPLAAKWGVYQNHHCFEMSTRDWAGLKKP",
    "gamma": "MGGSPLKVETAIDRNQHWYFCCEEKLPMSTAGRVNDQHILAWYKTFPSEEGCVLMRNI",
}
ALPHABET = pyhmmer.easel.Alphabet.amino()


def _build_hmm(name: str, gathering: float | None = None) -> pyhmmer.plan7.HMM:
    sequence = pyhmmer.easel.TextSequence(
        name=name.encode(), sequence=PROTEINS[name]
    ).digitize(ALPHABET)
    hmm, _, _ = pyhmmer.plan7.Builder(ALPHABET).build(
        sequence, pyhmmer.plan7.Background(ALPHABET)
    )
    if gathering is not None:
        hmm.cutoffs.gathering = (gathering, gathering)
    return hmm


def _write_hmms(path: Path, hmms: list[pyhmmer.plan7.HMM]) -> Path:
    with path.open("wb") as handle:
        for hmm in hmms:
            hmm.write(handle)
    return path


def _write_proteins(path: Path, names: tuple[str, ...] = tuple(PROTEINS)) -> Path:
    path.write_text(
        "".join(f">{name}\n{PROTEINS[name]}\n" for name in names), encoding="utf-8"
    )
    return path


@pytest.fixture
def mixed_hmm(tmp_path: Path) -> Path:
    """Alpha has a reachable GA cutoff, beta an unreachable one, gamma none."""
    return _write_hmms(
        tmp_path / "mixed.hmm",
        [_build_hmm("alpha", 20.0), _build_hmm("beta", 5000.0), _build_hmm("gamma")],
    )


def _found(domains: list[hmm_search.DomainHit]) -> set[tuple[str, str]]:
    return {(row.protein, row.model) for row in domains}


def test_search_sensitive_uses_the_evalue_for_every_model(
    tmp_path: Path, mixed_hmm: Path
) -> None:
    sequences = read_proteins(_write_proteins(tmp_path / "query.faa"))

    domains = search(load_profiles(mixed_hmm), sequences, threads=1, sensitive=True)

    assert _found(domains) == {
        ("alpha", "alpha"),
        ("beta", "beta"),
        ("gamma", "gamma"),
    }


def test_search_gathering_keeps_ga_hits_and_evalue_hits_in_a_mixed_set(
    tmp_path: Path, mixed_hmm: Path
) -> None:
    sequences = read_proteins(_write_proteins(tmp_path / "query.faa"))

    domains = search(load_profiles(mixed_hmm), sequences, threads=1, sensitive=False)

    assert _found(domains) == {("alpha", "alpha"), ("gamma", "gamma")}


def test_search_without_hits_returns_an_empty_list(tmp_path: Path) -> None:
    hmm_path = _write_hmms(tmp_path / "one.hmm", [_build_hmm("alpha")])
    sequences = read_proteins(_write_proteins(tmp_path / "query.faa", ("beta",)))

    domains = search(load_profiles(hmm_path), sequences, threads=1, sensitive=True)

    assert domains == []


def _fake_hit(name: str, *, reported: bool, domain_flags: list[bool]):
    domains = [
        SimpleNamespace(
            reported=flag,
            c_evalue=1e-10,
            i_evalue=1e-10,
            score=40.0 + index,
            bias=0.1,
            env_from=1,
            env_to=50,
            alignment=SimpleNamespace(
                hmm_from=1, hmm_to=50, target_from=2, target_to=49
            ),
        )
        for index, flag in enumerate(domain_flags)
    ]
    return SimpleNamespace(
        name=name,
        length=100,
        reported=reported,
        evalue=1e-12,
        score=60.0,
        bias=0.2,
        domains=domains,
    )


def test_search_drops_unreported_hits_and_domains(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    hits = [
        _fake_hit("kept", reported=True, domain_flags=[True, False, True]),
        _fake_hit("unreported_hit", reported=False, domain_flags=[True]),
        _fake_hit("no_reported_domain", reported=True, domain_flags=[False]),
    ]
    monkeypatch.setattr(
        hmm_search.pyhmmer, "hmmsearch", lambda *_args, **_kwargs: iter([hits])
    )
    model = SimpleNamespace(name="M1", accession=None, M=50)

    domains = search((model,), object(), threads=1, sensitive=True)

    assert [(row.protein, row.domain_index, row.domain_count) for row in domains] == [
        ("kept", 1, 2),
        ("kept", 2, 2),
    ]
    assert [row.domain_score for row in domains] == [40.0, 42.0]


def test_load_profiles_parses_the_file_once_for_two_calls(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hmm_path = _write_hmms(tmp_path / "one.hmm", [_build_hmm("alpha")])
    parse_calls: list[Path] = []
    parse = hmm_search._parse_profiles

    def counting_parse(path: Path):
        parse_calls.append(path)
        return parse(path)

    monkeypatch.setattr(hmm_search, "_parse_profiles", counting_parse)

    first = load_profiles(hmm_path)
    second = load_profiles(str(hmm_path))

    assert first is second
    assert len(parse_calls) == 1


def test_load_profiles_parses_again_after_the_file_changes(tmp_path: Path) -> None:
    hmm_path = _write_hmms(tmp_path / "models.hmm", [_build_hmm("alpha")])
    first = load_profiles(hmm_path)

    _write_hmms(hmm_path, [_build_hmm("alpha"), _build_hmm("beta")])
    os.utime(hmm_path, (time.time() + 5, time.time() + 5))
    second = load_profiles(hmm_path)

    assert len(first) == 1
    assert len(second) == 2


def test_load_profiles_rejects_a_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match=r"missing\.hmm"):
        load_profiles(tmp_path / "missing.hmm")


def test_run_marker_search_searches_a_first_protein_with_ambiguous_residues(
    tmp_path: Path,
) -> None:
    hmm_path = _write_hmms(tmp_path / "one.hmm", [_build_hmm("alpha")])
    faa_path = tmp_path / "ambiguous_first.faa"
    faa_path.write_text(
        ">ambiguous_first\nMGKKNRKKNNNDNNNDNKNDNKNDNKN\n"
        f">matching_protein\n{PROTEINS['alpha']}\n",
        encoding="utf-8",
    )

    filtered = run_marker_search(
        faa_path, hmm_path, tmp_path / "hmmout", threads=1, sensitive=True
    )

    assert "matching_protein" in filtered.read_text(encoding="utf-8")


def test_run_marker_search_forces_the_amino_alphabet(tmp_path: Path) -> None:
    hmm_path = _write_hmms(tmp_path / "one.hmm", [_build_hmm("alpha")])
    faa_path = tmp_path / "acgt_only.faa"
    faa_path.write_text(">acgt_only\n" + "ACGT" * 20 + "\n", encoding="utf-8")

    filtered = run_marker_search(
        faa_path, hmm_path, tmp_path / "hmmout", threads=1, sensitive=True
    )

    assert filtered.read_text(encoding="utf-8").startswith("# target name\taccession\t")


def test_run_marker_search_writes_the_three_tables(
    tmp_path: Path, mixed_hmm: Path
) -> None:
    faa_path = _write_proteins(tmp_path / "query.faa")
    out_dir = tmp_path / "hmmout"

    filtered = run_marker_search(
        faa_path, mixed_hmm, out_dir, threads=1, sensitive=False
    )

    assert filtered == out_dir / "models.out.filtered"
    assert (out_dir / "models.out").is_file()
    assert (out_dir / "models.counts").read_text(encoding="utf-8") == (
        "alpha\t1\ngamma\t1\n"
    )
    assert not (out_dir / "models.score").exists()


def test_run_marker_search_fails_loudly_and_writes_nothing_without_the_hmm_file(
    tmp_path: Path,
) -> None:
    faa_path = _write_proteins(tmp_path / "query.faa")
    out_dir = tmp_path / "hmmout"

    with pytest.raises(FileNotFoundError, match=r"missing\.hmm"):
        run_marker_search(
            faa_path, tmp_path / "missing.hmm", out_dir, threads=1, sensitive=True
        )

    assert not out_dir.exists()


def test_domain_table_coordinates_match_native_pyhmmer(tmp_path: Path) -> None:
    """The written coordinates equal pyhmmer's own domain table at the sequence ends."""
    hmm = _build_hmm("alpha")
    hmm_path = _write_hmms(tmp_path / "one.hmm", [hmm])
    faa_path = _write_proteins(tmp_path / "matching.faa", ("alpha",))
    sequences = digitize_proteins([("alpha", PROTEINS["alpha"])])
    native = io.BytesIO()
    next(pyhmmer.hmmsearch(hmm, sequences, cpus=1)).write(native, format="domains")
    native_row = next(
        line.split()
        for line in native.getvalue().decode().splitlines()
        if not line.startswith("#")
    )

    filtered = run_marker_search(
        faa_path, hmm_path, tmp_path / "hmmout", threads=1, sensitive=True
    )

    observed_row = next(
        line.split()
        for line in filtered.read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    )
    expected = [1, len(PROTEINS["alpha"])] * 3
    assert list(map(int, native_row[15:21])) == expected
    assert list(map(int, observed_row[15:21])) == expected


@pytest.mark.requires_db
@pytest.mark.slow
def test_sensitive_search_of_the_example_writes_only_reported_domains(
    tmp_path: Path,
) -> None:
    """Every written row passes ``E = domE = 1e-5``.

    ``domE`` applies to the conditional E-value, which is column 12.
    """
    hmm_path = RESOURCES_ROOT / COMBINED_HMM
    if not hmm_path.is_file():
        pytest.skip(f"GVClass database not installed at {hmm_path}")

    run_marker_search(EXAMPLE_PROTEINS, hmm_path, tmp_path, threads=1, sensitive=True)

    rows = [
        line.split("\t")
        for line in (tmp_path / "models.out").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    assert rows
    assert max(float(row[6]) for row in rows) <= 1e-5
    assert max(float(row[11]) for row in rows) <= 1e-5
