import gzip
import json
import math
from pathlib import Path

import numpy as np
import pytest

from statgen.reference import load_reference
from statgen.sumstats import load_sumstats, save_sumstats_cache, load_sumstats_cache
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, run_octave, skipif_no_octave

SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"


def _write_gz_tsv(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write(text)


def _valid_sumstats_text(include_optional: bool = True) -> str:
    if include_optional:
        return (
            "chr\tbp\ta1\ta2\tz\tn\tp\tbeta\tse\teaf\tinfo\n"
            "1\t100\tA\tG\t2.5\t1000\t0.01\t0.2\t0.1\t0.4\t0.9\n"
            "1\t200\tC\tT\t1.0\t950\t0.5\t-0.1\t0.2\t0.3\t0.8\n"
            "1\t300\tA\tC\t1.8\t1000\t0\t0.3\t0.15\t0.25\t0.95\n"
            "1\t400\tG\tA\t-1.2\t1000\t-0.1\t0.0\t0.3\t0.2\t0.7\n"
            "X\t100\tA\tG\t3.0\t500\t0.003\t0.5\t0.12\t0.45\t0.85\n"
            "X\t200\tC\tT\t0.5\t500\t0.6\t0.05\t0.25\t0.35\t0.75\n"
            "9\t999\tA\tG\t1.0\t1000\t0.3\t0.1\t0.1\t0.1\t0.1\n"
        )
    return (
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t2.5\t1000\t0.01\n"
        "1\t200\tC\tT\t1.0\t950\t0.5\n"
        "1\t300\tA\tC\t1.8\t1000\t0\n"
        "1\t400\tG\tA\t-1.2\t1000\t-0.1\n"
        "X\t100\tA\tG\t3.0\t500\t0.003\n"
        "X\t200\tC\tT\t0.5\t500\t0.6\n"
    )


def _load_cache_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def test_load_sumstats_alignment_and_accessors(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    assert s.num_snp == 8
    np.testing.assert_allclose(
        s.zvec,
        np.array([2.5, 1.0, 1.8, -1.2, np.nan, 3.0, 0.5, np.nan]),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        s.nvec,
        np.array([1000, 950, 1000, 1000, np.nan, 500, 500, np.nan]),
        equal_nan=True,
    )
    assert math.isclose(s.logpvec[0], 2.0)
    assert math.isclose(s.logpvec[1], -math.log10(0.5))
    assert math.isinf(s.logpvec[2])
    assert np.isnan(s.logpvec[3])
    assert np.isnan(s.logpvec[4])
    assert math.isclose(s.logpvec[5], -math.log10(0.003))
    assert math.isclose(s.logpvec[6], -math.log10(0.6))
    assert np.isnan(s.logpvec[7])

    assert s.beta_vec is not None
    assert s.se_vec is not None
    assert s.eaf_vec is not None
    assert s.info_vec is not None
    assert np.isnan(s.beta_vec[4]) and np.isnan(s.beta_vec[7])


def test_optional_fields_absent_use_sentinel(tmp_path):
    path = tmp_path / "traits_required_only.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=False))
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)
    assert s.beta_vec is None
    assert s.se_vec is None
    assert s.eaf_vec is None
    assert s.info_vec is None


@pytest.mark.parametrize("missing_col", ["chr", "bp", "a1", "a2", "z", "n"])
def test_missing_required_column_fails(tmp_path, missing_col):
    row = {"chr": "1", "bp": "100", "a1": "A", "a2": "G", "z": "1.2", "n": "900", "p": "0.1"}
    cols = [c for c in ("chr", "bp", "a1", "a2", "z", "n", "p") if c != missing_col]
    text = "\t".join(cols) + "\n" + "\t".join(row[c] for c in cols) + "\n"
    path = tmp_path / f"missing_{missing_col}.tsv.gz"
    _write_gz_tsv(path, text)
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match=f"missing required columns: {missing_col}"):
        load_sumstats(path, reference)


def test_strict_required_numeric_validation(tmp_path):
    bad_z = tmp_path / "bad_z.tsv.gz"
    _write_gz_tsv(
        bad_z,
        "chr\tbp\ta1\ta2\tz\tn\n"
        "1\t100\tA\tG\tNA\t1000\n",
    )
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match="z must be finite numeric"):
        load_sumstats(bad_z, reference)

    bad_n = tmp_path / "bad_n.tsv.gz"
    _write_gz_tsv(
        bad_n,
        "chr\tbp\ta1\ta2\tz\tn\n"
        "1\t100\tA\tG\t1.5\tInf\n",
    )
    with pytest.raises(ValueError, match="n must be finite numeric"):
        load_sumstats(bad_n, reference)


def test_p_edge_cases(tmp_path):
    path = tmp_path / "p_edge_cases.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t1.0\t1000\t1\n"
        "1\t200\tC\tT\t1.0\t1000\t-0.1\n"
        "1\t300\tA\tC\t1.0\t1000\t1.1\n"
        "1\t400\tG\tA\t1.0\t1000\t\n"
        "X\t100\tA\tG\t1.0\t1000\t0\n"
        "X\t200\tC\tT\t1.0\t1000\t0.2\n",
    )
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    assert s.logpvec[0] == 0.0
    assert np.isnan(s.logpvec[1])
    assert np.isnan(s.logpvec[2])
    assert np.isnan(s.logpvec[3])
    assert np.isnan(s.logpvec[4])
    assert math.isinf(s.logpvec[5])
    assert math.isclose(s.logpvec[6], -math.log10(0.2))
    assert np.isnan(s.logpvec[7])


def test_allele_flip_does_not_match_reference_key(tmp_path):
    path = tmp_path / "allele_flip.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tG\tA\t7.0\t1000\t0.01\n"  # swapped vs reference 1:100:A:G
        "X\t100\tA\tG\t3.0\t500\t0.01\n",
    )
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    assert np.isnan(s.zvec[0])
    assert np.isnan(s.nvec[0])
    assert np.isnan(s.logpvec[0])
    assert s.zvec[5] == 3.0


def test_select_shards(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)
    x_only = s.select_shards(["X"])
    assert [sh.label for sh in x_only.shards] == ["X"]
    assert x_only.num_snp == 3
    np.testing.assert_allclose(x_only.zvec, np.array([3.0, 0.5, np.nan]), equal_nan=True)


def test_cache_roundtrip_and_subset(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    cache = tmp_path / "sumstats_cache.npz"
    save_sumstats_cache(s, cache)
    loaded = load_sumstats_cache(cache)
    np.testing.assert_allclose(loaded.zvec, s.zvec, equal_nan=True)
    np.testing.assert_allclose(loaded.logpvec, s.logpvec, equal_nan=True)
    assert reference.is_object_compatible(loaded) is True

    loaded_x = load_sumstats_cache(cache, shards=["X"])
    assert [sh.label for sh in loaded_x.shards] == ["X"]
    assert loaded_x.num_snp == 3


def test_cache_validation_errors(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)
    cache = tmp_path / "sumstats_cache.npz"
    save_sumstats_cache(s, cache)

    # Bad schema.
    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["schema"] = "sumstats_cache/bad"
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_schema = tmp_path / "bad_schema.npz"
    np.savez_compressed(bad_schema, **arrays)
    with pytest.raises(ValueError, match="Unsupported sumstats cache schema"):
        load_sumstats_cache(bad_schema)

    # Mismatched shard label/checksum lengths.
    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["shard_checksums"] = meta["shard_checksums"][:-1]
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_lengths = tmp_path / "bad_lengths.npz"
    np.savez_compressed(bad_lengths, **arrays)
    with pytest.raises(ValueError, match="shard_labels and shard_checksums length mismatch"):
        load_sumstats_cache(bad_lengths)

    # Invalid shard request.
    with pytest.raises(ValueError, match="requested shard '2' is not present"):
        load_sumstats_cache(cache, shards=["2"])


def _octave_script(expr: str) -> str:
    matlab_dir = str(MATLAB_DIR)
    fixture_dir = str(FIXTURES_DIR)
    return (
        f"addpath('{matlab_dir}'); "
        f"fixture_dir = '{fixture_dir}'; "
        + expr
    )


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_sumstats_fixture_fails_required_z():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "try; "
        "  s = statgen.load_sumstats([fixture_dir '/sumstats/traits.tsv.gz'], ref); "
        "  fprintf('NOFAIL\\n'); "
        "catch; "
        "  fprintf('FAIL\\n'); "
        "end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert "FAIL" in result.stdout


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_roundtrip(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    cache = tmp_path / "sumstats_cache.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        f"statgen.save_sumstats_cache(s, '{cache}'); "
        f"s2 = statgen.load_sumstats_cache('{cache}'); "
        "fprintf('%d\\n', s2.num_snp); "
        "fprintf('%d\\n', ref.is_object_compatible(s2)); "
        "fprintf('%.6f\\n', s2.zvec(1));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "8"
    assert lines[1] == "1"
    assert lines[2] == "2.500000"


@pytest.mark.octave
@skipif_no_octave
def test_octave_optional_fields_absent_use_empty_vector_sentinel(tmp_path):
    path = tmp_path / "traits_required_only.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=False))
    cache = tmp_path / "sumstats_required_only_cache.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%d %d %d %d\\n', isempty(s.beta_vec), isempty(s.se_vec), isempty(s.eaf_vec), isempty(s.info_vec)); "
        f"statgen.save_sumstats_cache(s, '{cache}'); "
        f"s2 = statgen.load_sumstats_cache('{cache}'); "
        "fprintf('%d %d %d %d\\n', isempty(s2.beta_vec), isempty(s2.se_vec), isempty(s2.eaf_vec), isempty(s2.info_vec));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1 1 1 1"
    assert lines[1] == "1 1 1 1"
