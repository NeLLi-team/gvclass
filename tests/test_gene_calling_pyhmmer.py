import io
from pathlib import Path

import pyhmmer
import pytest

from src.core.gene_calling import run_hmmsearch, run_pyhmmer_search_simple
from src.core.hmm_search import run_pyhmmer_search
from src.core.hmm_search_multi import run_pyhmmer_search_multi


MODEL_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLQDKPEAIRRLVDELAQLGVE"


@pytest.fixture
def tiny_amino_hmm(tmp_path: Path) -> Path:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(
        name=b"tiny_model",
        sequence=MODEL_PROTEIN,
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        sequence,
        pyhmmer.plan7.Background(alphabet),
    )

    hmm_path = tmp_path / "tiny.hmm"
    with hmm_path.open("wb") as handle:
        hmm.write(handle)
    return hmm_path


def test_run_hmmsearch_handles_ambiguous_first_protein(
    tmp_path: Path,
    tiny_amino_hmm: Path,
) -> None:
    faa_path = tmp_path / "ambiguous_first.faa"
    faa_path.write_text(
        ">ambiguous_first\n"
        "MGKKNRKKNNNDNNNDNKNDNKNDNKN\n"
        ">matching_protein\n"
        f"{MODEL_PROTEIN}\n"
    )

    hit_count, *_ = run_hmmsearch(str(faa_path), str(tiny_amino_hmm))

    assert hit_count >= 1
    assert "matching_protein" in faa_path.with_suffix(".hmmout").read_text()


def test_run_hmmsearch_forces_amino_for_acgt_only_protein(
    tmp_path: Path,
    tiny_amino_hmm: Path,
) -> None:
    faa_path = tmp_path / "acgt_only.faa"
    faa_path.write_text(">acgt_only\n" + "ACGT" * 20 + "\n")

    run_hmmsearch(str(faa_path), str(tiny_amino_hmm))

    header = faa_path.with_suffix(".hmmout").read_text().splitlines()[0]
    assert header.startswith("# target name\taccession\t")


@pytest.mark.parametrize("writer", ["gene_calling", "single_hmm", "multi_hmm"])
def test_domain_coordinates_match_native_pyhmmer(
    tmp_path: Path, tiny_amino_hmm: Path, writer: str
) -> None:
    """All runtime writers preserve native coordinates at sequence boundaries."""
    faa_path = tmp_path / "matching.faa"
    faa_path.write_text(f">matching_protein\n{MODEL_PROTEIN}\n")
    with pyhmmer.plan7.HMMFile(tiny_amino_hmm) as handle:
        hmm = next(handle)
    with pyhmmer.easel.SequenceFile(
        faa_path, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()
    ) as handle:
        sequences = list(handle)
    hits = next(pyhmmer.hmmsearch(hmm, sequences, cpus=1))
    native = io.BytesIO()
    hits.write(native, format="domains")
    native_row = next(
        line.split()
        for line in native.getvalue().decode().splitlines()
        if not line.startswith("#")
    )
    expected = [1, len(MODEL_PROTEIN)] * 3
    assert list(map(int, native_row[15:21])) == expected

    output_path = tmp_path / f"{writer}.domtblout"
    if writer == "gene_calling":
        run_pyhmmer_search_simple(str(tiny_amino_hmm), str(faa_path), str(output_path))
    elif writer == "single_hmm":
        run_pyhmmer_search(str(tiny_amino_hmm), str(faa_path), str(output_path), threads=1)
    else:
        run_pyhmmer_search_multi(
            [str(tiny_amino_hmm)], str(faa_path), str(output_path), threads=1
        )
    observed_row = next(
        line.split()
        for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    )
    assert list(map(int, observed_row[15:21])) == expected
