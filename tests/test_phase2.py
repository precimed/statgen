import gzip
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from statgen._utils import allele_hash64
from statgen.reference import load_reference
from statgen.sumstats import create_sumstats, load_sumstats, save_sumstats_cache, load_sumstats_cache
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, matlab_data_lines, matlab_output_text, run_octave, skipif_no_octave

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
            "1\t400\tG\tA\t-1.2\t1000\t0.2\t0.0\t0.3\t0.2\t0.7\n"
            "X\t100\tA\tG\t3.0\t500\t0.003\t0.5\t0.12\t0.45\t0.85\n"
            "X\t200\tC\tT\t0.5\t500\t0.6\t0.05\t0.25\t0.35\t0.75\n"
            "9\t999\tA\tG\t1.0\t1000\t0.3\t0.1\t0.1\t0.1\t0.1\n"
        )
    return (
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\t0.01\n"
        "1\t200\tC\tT\t0.5\n"
        "1\t300\tA\tC\t0\n"
        "1\t400\tG\tA\t0.2\n"
        "X\t100\tA\tG\t0.003\n"
        "X\t200\tC\tT\t0.6\n"
    )


def _genomatch_sumstats_text() -> str:
    return (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tZ\tN\tP\tBETA\tSE\tEAF\tINFO\tDirection\n"
        "1\t100\trs1\tA\tG\t2.5\t1000\t0.01\t0.2\t0.1\t0.4\t0.9\t+\n"
        "1\t200\trs2\tC\tT\t1.0\t950\t0.5\t-0.1\t0.2\t0.3\t0.8\t-\n"
        "X\t100\trsx1\tA\tG\t3.0\t500\t0.003\t0.5\t0.12\t0.45\t0.85\t+\n"
    )


def _mixed_case_vmap_sumstats_text() -> str:
    return (
        "ChR\tBp\tA1\ta2\tZ\tn\tP\n"
        "1\t100\tA\tG\t2.5\t1000\t0.01\n"
        "X\t100\tA\tG\t3.0\t500\t0.003\n"
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
    assert math.isclose(s.logpvec[3], -math.log10(0.2))
    assert np.isnan(s.logpvec[4])
    assert math.isclose(s.logpvec[5], -math.log10(0.003))
    assert math.isclose(s.logpvec[6], -math.log10(0.6))
    assert np.isnan(s.logpvec[7])
    np.testing.assert_array_equal(
        s.is_present,
        np.array([True, True, True, True, False, True, True, False]),
    )

    assert s.beta_vec is not None
    assert s.se_vec is not None
    assert s.eaf_vec is not None
    assert s.info_vec is not None
    assert np.isnan(s.beta_vec[4]) and np.isnan(s.beta_vec[7])


def test_allele_hash64_known_values():
    hashes = allele_hash64(["A", "C", "G", "T", "N", "-", "ACGT", "ATCGGCTA"])
    assert hashes.dtype == np.uint64
    assert [f"{int(x):016X}" for x in hashes] == [
        "0000014300000149",
        "000001450000014B",
        "000001490000014F",
        "000001560000015C",
        "0000015000000156",
        "0000012F00000135",
        "47119B266503B322",
        "1A8AE73A00B4D922",
    ]


def test_load_sumstats_duplicate_matching_key_fails(tmp_path):
    path = tmp_path / "duplicate_key.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t2.5\t1000\t0.01\n"
        "1\t100\tA\tG\t2.6\t1000\t0.02\n",
    )
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match="Ambiguous duplicate sumstats/reference matching key"):
        load_sumstats(path, reference)


def test_load_sumstats_orients_swapped_alleles(tmp_path, capsys):
    path = tmp_path / "swapped.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t1.0\t1000\t0.01\n"
        "1\t200\tT\tC\t1.0\t1000\t0.02\n",
    )
    reference = load_reference(SHARDED_REF, shards=["1"])
    s = load_sumstats(path, reference)
    assert "shard 1: 1 successful allele swap" in capsys.readouterr().err
    assert s.zvec[1] == -1.0

    np.testing.assert_array_equal(
        s.is_present,
        np.array([True, True, False, False, False]),
    )


@pytest.mark.parametrize(
    "text_factory",
    [_genomatch_sumstats_text, _mixed_case_vmap_sumstats_text],
)
def test_load_sumstats_canonicalizes_supported_headers(tmp_path, text_factory):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, text_factory())
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    assert s.zvec[0] == 2.5
    assert s.nvec[0] == 1000
    assert math.isclose(s.logpvec[0], 2.0)
    assert s.zvec[5] == 3.0
    assert s.nvec[5] == 500
    assert math.isclose(s.logpvec[5], -math.log10(0.003))


def test_duplicate_columns_after_column_normalization_fail(tmp_path):
    path = tmp_path / "duplicate.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tCHR\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t1\t100\tA\tG\t1.0\t1000\t0.1\n",
    )
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match="duplicate columns after column normalization: chr"):
        load_sumstats(path, reference)


def test_optional_fields_absent_use_sentinel(tmp_path):
    path = tmp_path / "traits_required_only.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=False))
    reference = load_reference(SHARDED_REF)
    with pytest.warns(RuntimeWarning) as warn:
        s = load_sumstats(path, reference)
    messages = [str(w.message) for w in warn]
    assert any("zvec is absent" in msg for msg in messages)
    assert any("nvec is absent" in msg for msg in messages)
    assert s.zvec is None
    assert s.nvec is None
    assert s.beta_vec is None
    assert s.se_vec is None
    assert s.eaf_vec is None
    assert s.info_vec is None
    cache = tmp_path / "sumstats_required_only_cache.npz"
    s.save_cache(cache)
    with pytest.warns(RuntimeWarning) as cache_warn:
        loaded = load_sumstats_cache(cache)
    cache_messages = [str(w.message) for w in cache_warn]
    assert any("zvec is absent" in msg for msg in cache_messages)
    assert any("nvec is absent" in msg for msg in cache_messages)
    assert loaded.zvec is None
    assert loaded.nvec is None


def test_sumstats_warning_stacklevels(tmp_path):
    path = tmp_path / "traits_required_only.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=False))
    reference = load_reference(SHARDED_REF)

    with warnings.catch_warnings(record=True) as load_warn:
        warnings.simplefilter("always")
        s = load_sumstats(path, reference)
    assert load_warn
    assert all(Path(w.filename).name == "sumstats.py" for w in load_warn)
    assert all(w.lineno == load_sumstats.__code__.co_firstlineno + 3 for w in load_warn)

    cache = tmp_path / "sumstats_required_only_cache.npz"
    s.save_cache(cache)
    with warnings.catch_warnings(record=True) as cache_warn:
        warnings.simplefilter("always")
        load_sumstats_cache(cache)
    assert cache_warn
    assert all(Path(w.filename).name == "sumstats.py" for w in cache_warn)
    assert all(w.lineno == load_sumstats_cache.__code__.co_firstlineno + 42 for w in cache_warn)


@pytest.mark.parametrize("missing_col", ["chr", "bp", "a1", "a2", "p"])
def test_missing_required_column_fails(tmp_path, missing_col):
    row = {"chr": "1", "bp": "100", "a1": "A", "a2": "G", "z": "1.2", "n": "900", "p": "0.1"}
    cols = [c for c in ("chr", "bp", "a1", "a2", "z", "n", "p") if c != missing_col]
    text = "\t".join(cols) + "\n" + "\t".join(row[c] for c in cols) + "\n"
    path = tmp_path / f"missing_{missing_col}.tsv.gz"
    _write_gz_tsv(path, text)
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match=f"missing required columns: {missing_col}"):
        load_sumstats(path, reference)


@pytest.mark.parametrize(
    ("text", "pattern"),
    [
        ("chr\tbp\ta1\ta2\tp\nchr1\t100\tA\tG\t0.1\n", "chr-style labels"),
        ("chr\tbp\ta1\ta2\tp\n1\t100\ta\tG\t0.1\n", "uppercase DNA bases"),
        ("chr\tbp\ta1\ta2\tp\n1\t100\tN\tG\t0.1\n", "uppercase DNA bases"),
        ("chr\tbp\ta1\ta2\tp\n1\t100\tA\tA\t0.1\n", "a1 and a2 must differ"),
    ],
)
def test_sumstats_rejects_invalid_chr_and_alleles(tmp_path, text, pattern):
    path = tmp_path / "bad_sumstats.tsv.gz"
    _write_gz_tsv(path, text)
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match=pattern):
        load_sumstats(path, reference)


def test_sumstats_ignores_y_mt_before_field_validation(tmp_path):
    path = tmp_path / "ignored_y_mt.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tp\n"
        "Y\t1\tN\tN\tNaN\n"
        "MT\t1\tN\tN\tNaN\n"
        "1\t100\tA\tG\t0.1\n",
    )
    reference = load_reference(SHARDED_REF).select_shards(["1"])
    s = load_sumstats(path, reference)
    assert s.num_snp == 5
    assert int(np.count_nonzero(s.is_present)) == 1


def test_optional_z_n_missing_values_warn(tmp_path):
    path = tmp_path / "missing_zn_values.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\tNaN\t1000\t0.1\n"
        "1\t200\tC\tT\t1.5\tInf\t0.2\n",
    )
    reference = load_reference(SHARDED_REF)
    with pytest.warns(RuntimeWarning) as warn:
        s = load_sumstats(path, reference)
    messages = [str(w.message) for w in warn]
    assert any("zvec has missing values" in msg for msg in messages)
    assert any("nvec has missing values" in msg for msg in messages)
    assert np.isnan(s.zvec[0])
    assert np.isnan(s.nvec[1])

    cache = tmp_path / "sumstats_missing_zn_cache.npz"
    s.save_cache(cache)
    with pytest.warns(RuntimeWarning) as cache_warn:
        loaded = load_sumstats_cache(cache)
    cache_messages = [str(w.message) for w in cache_warn]
    assert any("zvec has missing values" in msg for msg in cache_messages)
    assert any("nvec has missing values" in msg for msg in cache_messages)
    np.testing.assert_allclose(loaded.zvec, s.zvec, equal_nan=True)
    np.testing.assert_allclose(loaded.nvec, s.nvec, equal_nan=True)


def test_p_edge_cases(tmp_path):
    path = tmp_path / "p_edge_cases.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t1.0\t1000\t1\n"
        "1\t200\tC\tT\t1.0\t1000\t0.5\n"
        "1\t300\tA\tC\t1.0\t1000\t0.25\n"
        "1\t400\tG\tA\t1.0\t1000\t0.1\n"
        "X\t100\tA\tG\t1.0\t1000\t0\n"
        "X\t200\tC\tT\t1.0\t1000\t0.2\n",
    )
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    assert s.logpvec[0] == 0.0
    assert math.isclose(s.logpvec[1], -math.log10(0.5))
    assert math.isclose(s.logpvec[2], -math.log10(0.25))
    assert math.isclose(s.logpvec[3], 1.0)
    assert np.isnan(s.logpvec[4])
    assert math.isinf(s.logpvec[5])
    assert math.isclose(s.logpvec[6], -math.log10(0.2))
    assert np.isnan(s.logpvec[7])


@pytest.mark.parametrize("bad_p", ["-0.1", "1.1", "Inf"])
def test_invalid_p_values_fail(tmp_path, bad_p):
    path = tmp_path / "bad_p.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        f"1\t100\tA\tG\t1.0\t1000\t{bad_p}\n",
    )
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match=r"p must be finite numeric in \[0, 1\]"):
        load_sumstats(path, reference)


def test_allele_flip_orients_values_to_reference(tmp_path):
    path = tmp_path / "allele_flip.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tG\tA\t7.0\t1000\t0.01\n"  # swapped vs reference 1:100:A:G
        "X\t100\tA\tG\t3.0\t500\t0.01\n",
    )
    reference = load_reference(SHARDED_REF)
    s = load_sumstats(path, reference)

    assert s.zvec[0] == -7.0
    assert s.nvec[0] == 1000
    assert s.logpvec[0] == 2.0
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
    s.save_cache(cache)
    loaded = load_sumstats_cache(cache)
    np.testing.assert_allclose(loaded.zvec, s.zvec, equal_nan=True)
    np.testing.assert_allclose(loaded.nvec, s.nvec, equal_nan=True)
    np.testing.assert_allclose(loaded.logpvec, s.logpvec, equal_nan=True)
    np.testing.assert_array_equal(loaded.is_present, s.is_present)
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


def test_create_sumstats_from_vectors():
    reference = load_reference(SHARDED_REF)
    zvec = np.array([2.5, 1.0, 1.8, -1.2, np.nan, 3.0, 0.5, np.nan], dtype=float)
    nvec = np.array([1000, 950, 1000, 1000, np.nan, 500, 500, np.nan], dtype=float)
    p = np.array([0.01, 0.5, 0.0, 0.2, 1.0, 0.003, 0.6, 0.9], dtype=float)
    beta = np.array([0.2, -0.1, 0.3, 0.0, np.nan, 0.5, 0.05, np.nan], dtype=float)

    s = create_sumstats(reference, p, z=zvec, n=nvec, beta=beta)

    np.testing.assert_allclose(s.zvec, zvec, equal_nan=True)
    np.testing.assert_allclose(s.nvec, nvec, equal_nan=True)
    np.testing.assert_array_equal(s.is_present, np.ones(reference.num_snp, dtype=bool))
    assert math.isclose(s.logpvec[0], 2.0)
    assert math.isclose(s.logpvec[1], -math.log10(0.5))
    assert math.isinf(s.logpvec[2])
    assert math.isclose(s.logpvec[3], -math.log10(0.2))
    assert s.logpvec[4] == 0.0
    assert math.isclose(s.logpvec[5], -math.log10(0.003))
    assert math.isclose(s.logpvec[6], -math.log10(0.6))
    assert math.isclose(s.logpvec[7], -math.log10(0.9))
    np.testing.assert_allclose(s.beta_vec, beta, equal_nan=True)
    assert s.se_vec is None
    assert s.eaf_vec is None
    assert s.info_vec is None


def test_create_sumstats_validation_errors():
    reference = load_reference(SHARDED_REF)
    p = np.full(reference.num_snp, 0.5)
    with pytest.raises(ValueError, match="p length mismatch"):
        create_sumstats(reference, np.array([1.0]))
    with pytest.raises(ValueError, match="z length mismatch"):
        create_sumstats(reference, p, z=np.array([1.0]))
    with pytest.raises(ValueError, match="n\\[0\\] must be finite numeric or NaN"):
        create_sumstats(reference, p, n=np.r_[np.inf, np.zeros(reference.num_snp - 1)])
    with pytest.raises(ValueError, match=r"p\[0\] must be finite numeric in \[0, 1\]"):
        create_sumstats(reference, np.r_[np.inf, np.zeros(reference.num_snp - 1)])
    with pytest.raises(ValueError, match=r"p\[0\] must be finite numeric in \[0, 1\]"):
        create_sumstats(reference, np.r_[-0.1, np.zeros(reference.num_snp - 1)])
    with pytest.raises(ValueError, match="beta\\[0\\] must be finite numeric or NaN"):
        create_sumstats(
            reference,
            p,
            beta=np.r_[np.inf, np.zeros(reference.num_snp - 1)],
        )


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
def test_octave_allele_hash64_matches_python():
    alleles = ["A", "C", "G", "T", "N", "-", "ACGT", "ATCGGCTA"]
    expected = [f"{int(x):016X}" for x in allele_hash64(alleles)]
    script = _octave_script(
        "h = statgen.internal.allele_hash64({'A','C','G','T','N','-','ACGT','ATCGGCTA'}); "
        "for i = 1:numel(h); fprintf('%s\\n', dec2hex(h(i), 16)); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines() == expected


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_sumstats_duplicate_matching_key_fails(tmp_path):
    path = tmp_path / "duplicate_key.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t777\tACGT\tATCGGCTA\t2.5\t1000\t0.01\n"
        "1\t777\tACGT\tATCGGCTA\t2.6\t1000\t0.02\n",
    )
    reference = tmp_path / "reference.bim"
    reference.write_text("1\trs777\t0\t777\tACGT\tATCGGCTA\n")
    script = _octave_script(
        f"ref = statgen.load_reference('{reference}'); "
        f"try; statgen.load_sumstats('{path}', ref); fprintf('NOFAIL\\n'); catch ME; fprintf('%s\\n', ME.message); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert "a1_hash64=5121044840502768418" in result.stdout
    assert "a2_hash64=1912595228078168354" in result.stdout


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_roundtrip(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    cache = tmp_path / "sumstats_cache.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        f"s.save_cache('{cache}'); "
        f"s2 = statgen.load_sumstats_cache('{cache}'); "
        "fprintf('%d\\n', s2.num_snp); "
        "fprintf('%d\\n', ref.is_object_compatible(s2)); "
        "fprintf('%.6f\\n', s2.zvec(1)); "
        f"L = load('{cache}'); "
        "fprintf('%d %d %d\\n', isfield(L, 'metadata'), isfield(L, 'zvec'), isfield(L, 'cache_shards'));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "8"
    assert lines[1] == "1"
    assert lines[2] == "2.500000"
    assert lines[3] == "1 1 0"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_canonicalizes_supported_headers(tmp_path):
    path = tmp_path / "traits_genomatch.tsv.gz"
    _write_gz_tsv(path, _genomatch_sumstats_text())
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%.6f\\n', s.zvec(1)); "
        "fprintf('%.6f\\n', s.nvec(1)); "
        "fprintf('%.6f\\n', s.zvec(6)); "
        "fprintf('%.6f\\n', s.beta_vec(1));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == ["2.500000", "1000.000000", "3.000000", "0.200000"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_gzip_uses_statgen_scratch(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=False))
    scratch = tmp_path / "scratch"
    script = _octave_script(
        f"setenv('STATGEN_SCRATCH', '{scratch}'); "
        f"ref = statgen.load_reference('{SHARDED_REF}'); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%d\\n', s.num_snp); "
        f"d = dir('{scratch}'); "
        "names = {d.name}; "
        "names = names(~strcmp(names, '.') & ~strcmp(names, '..')); "
        "fprintf('%d\\n', numel(names)); "
        "setenv('STATGEN_SCRATCH', '');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "8"
    assert lines[1] == "0"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_gzip_parse_error_cleans_scratch(tmp_path):
    path = tmp_path / "bad_numeric_token.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\tinvalid\n",
    )
    scratch = tmp_path / "scratch"
    script = _octave_script(
        f"setenv('STATGEN_SCRATCH', '{scratch}'); "
        f"ref = statgen.load_reference('{SHARDED_REF}'); "
        f"ok = 0; try; statgen.load_sumstats('{path}', ref); catch; ok = 1; end; "
        f"d = dir('{scratch}'); "
        "names = {d.name}; "
        "names = names(~strcmp(names, '.') & ~strcmp(names, '..')); "
        "fprintf('%d %d\\n', ok, numel(names)); "
        "setenv('STATGEN_SCRATCH', '');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 0"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_gzip_default_scratch_does_not_use_fixture_dir():
    source_dir = FIXTURES_DIR / "sumstats"
    before = {p.name for p in source_dir.glob("sumstats_gunzip_*") if p.is_dir()}
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "s = statgen.load_sumstats([fixture_dir '/sumstats/traits.tsv.gz'], ref); "
        "fprintf('%d\\n', s.num_snp);"
    )
    result = run_octave(script)
    after = {p.name for p in source_dir.glob("sumstats_gunzip_*") if p.is_dir()}
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "8"
    assert after == before


@pytest.mark.octave
@skipif_no_octave
def test_octave_missing_required_columns_fail(tmp_path):
    base_row = {"chr": "1", "bp": "100", "a1": "A", "a2": "G", "z": "1.2", "n": "900", "p": "0.1"}
    missing_cols = ["chr", "bp", "a1", "a2", "p"]
    paths = []
    for missing_col in missing_cols:
        cols = [c for c in ("chr", "bp", "a1", "a2", "z", "n", "p") if c != missing_col]
        text = "\t".join(cols) + "\n" + "\t".join(base_row[c] for c in cols) + "\n"
        p = tmp_path / f"missing_{missing_col}.tsv.gz"
        _write_gz_tsv(p, text)
        paths.append(str(p))

    paths_expr = "{" + ", ".join(f"'{p}'" for p in paths) + "}"
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"paths = {paths_expr}; "
        "ok = zeros(1, numel(paths)); "
        "for i = 1:numel(paths); "
        "  try; "
        "    statgen.load_sumstats(paths{i}, ref); "
        "    ok(i) = 0; "
        "  catch; "
        "    ok(i) = 1; "
        "  end; "
        "end; "
        "fprintf('%d %d %d %d %d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_rejects_invalid_chr_and_alleles(tmp_path):
    cases = {
        "chr_style": (
            "chr\tbp\ta1\ta2\tp\nchr1\t100\tA\tG\t0.1\n",
            "chr-style labels",
        ),
        "lowercase": (
            "chr\tbp\ta1\ta2\tp\n1\t100\ta\tG\t0.1\n",
            "uppercase DNA bases",
        ),
        "non_dna": (
            "chr\tbp\ta1\ta2\tp\n1\t100\tN\tG\t0.1\n",
            "uppercase DNA bases",
        ),
        "same_alleles": (
            "chr\tbp\ta1\ta2\tp\n1\t100\tA\tA\t0.1\n",
            "a1 and a2 must differ",
        ),
    }
    paths = []
    expected = []
    for name, (text, msg) in cases.items():
        path = tmp_path / f"{name}.tsv.gz"
        _write_gz_tsv(path, text)
        paths.append(str(path))
        expected.append(msg)

    paths_expr = "{" + ", ".join(f"'{p}'" for p in paths) + "}"
    expected_expr = "{" + ", ".join(f"'{m}'" for m in expected) + "}"
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"paths = {paths_expr}; "
        f"expected = {expected_expr}; "
        "ok = zeros(1, numel(paths)); "
        "for i = 1:numel(paths); "
        "  try; "
        "    statgen.load_sumstats(paths{i}, ref); "
        "  catch ME; "
        "    ok(i) = ~isempty(strfind(ME.message, expected{i})); "
        "  end; "
        "end; "
        "fprintf('%d %d %d %d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_ignores_y_mt_before_field_validation(tmp_path):
    path = tmp_path / "ignored_y_mt.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tp\n"
        "Y\t1\tN\tN\tNaN\n"
        "MT\t1\tN\tN\tNaN\n"
        "1\t100\tA\tG\t0.1\n",
    )
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim'], {'1'}); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%d %d\\n', s.num_snp, sum(s.is_present));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "5 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_orients_swapped_alleles(tmp_path):
    path = tmp_path / "swapped.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t1.0\t1000\t0.01\n"
        "1\t200\tT\tC\t1.0\t1000\t0.02\n",
    )
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim'], {'1'}); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%d ', s.is_present); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    combined = matlab_output_text(result.stdout + result.stderr)
    assert "shard 1: 1 successful allele swap(s)" in combined
    assert "remain unmatched" not in combined
    assert matlab_data_lines(result.stdout)[-1] == "1 1 0 0 0"


@pytest.mark.octave
@skipif_no_octave
def test_octave_optional_z_n_missing_values_warn(tmp_path):
    path = tmp_path / "missing_zn_values.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\tNaN\t1000\t0.1\n"
        "1\t200\tC\tT\t1.5\tInf\t0.2\n",
    )
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%d %d\\n', isnan(s.zvec(1)), isnan(s.nvec(2)));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_explicit_nan_tokens(tmp_path):
    optional_nan = tmp_path / "optional_nan.tsv.gz"
    _write_gz_tsv(
        optional_nan,
        "chr\tbp\ta1\ta2\tz\tn\tp\tbeta\tse\teaf\tinfo\n"
        "1\t100\tA\tG\t1.0\t1000\t0.1\tNaN\tNaN\tNaN\tNaN\n",
    )
    required_nan = tmp_path / "required_nan.tsv.gz"
    _write_gz_tsv(
        required_nan,
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\tNaN\n",
    )
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{optional_nan}', ref); "
        "ok = [isnan(s.beta_vec(1)), isnan(s.se_vec(1)), isnan(s.eaf_vec(1)), isnan(s.info_vec(1))]; "
        f"missing = statgen.load_sumstats('{required_nan}', ref); required_missing = all(isnan(missing.logpvec)); "
        "fprintf('%d %d %d %d %d\\n', ok, required_missing);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_p_edge_cases(tmp_path):
    path = tmp_path / "p_edge_cases.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tA\tG\t1.0\t1000\t1\n"
        "1\t200\tC\tT\t1.0\t1000\t0.5\n"
        "1\t300\tA\tC\t1.0\t1000\t0.25\n"
        "1\t400\tG\tA\t1.0\t1000\t0.1\n"
        "X\t100\tA\tG\t1.0\t1000\t0\n"
        "X\t200\tC\tT\t1.0\t1000\t0.2\n",
    )
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "v = s.logpvec; "
        "ok = zeros(1, 8); "
        "ok(1) = (v(1) == 0); "
        "ok(2) = abs(v(2) - (-log10(0.5))) < 1e-12; "
        "ok(3) = abs(v(3) - (-log10(0.25))) < 1e-12; "
        "ok(4) = abs(v(4) - 1.0) < 1e-12; "
        "ok(5) = isnan(v(5)); "
        "ok(6) = isinf(v(6)) && v(6) > 0; "
        "ok(7) = abs(v(7) - (-log10(0.2))) < 1e-12; "
        "ok(8) = isnan(v(8)); "
        "fprintf('%d %d %d %d %d %d %d %d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1 1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sumstats_is_present_from_load(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\\n', s.is_present);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1 0 1 1 0"


@pytest.mark.octave
@skipif_no_octave
def test_octave_invalid_p_values_fail(tmp_path):
    paths = []
    for i, bad_p in enumerate(["-0.1", "1.1", "invalid"]):
        path = tmp_path / f"bad_p_{i}.tsv.gz"
        _write_gz_tsv(
            path,
            "chr\tbp\ta1\ta2\tz\tn\tp\n"
            f"1\t100\tA\tG\t1.0\t1000\t{bad_p}\n",
        )
        paths.append(str(path))
    paths_expr = "{" + ", ".join(f"'{p}'" for p in paths) + "}"
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"paths = {paths_expr}; "
        "ok = zeros(1, numel(paths)); "
        "for i = 1:numel(paths); "
        "  try; statgen.load_sumstats(paths{i}, ref); catch; ok(i) = 1; end; "
        "end; "
        "fprintf('%d %d %d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_allele_flip_orients_values_to_reference(tmp_path):
    path = tmp_path / "allele_flip.tsv.gz"
    _write_gz_tsv(
        path,
        "chr\tbp\ta1\ta2\tz\tn\tp\n"
        "1\t100\tG\tA\t7.0\t1000\t0.01\n"
        "X\t100\tA\tG\t3.0\t500\t0.01\n",
    )
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "ok1 = s.zvec(1) == -7; "
        "ok2 = s.nvec(1) == 1000; "
        "ok3 = s.logpvec(1) == 2; "
        "ok4 = abs(s.zvec(6) - 3.0) < 1e-12; "
        "fprintf('%d %d %d %d\\n', ok1, ok2, ok3, ok4);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_cache_validation_errors(tmp_path):
    path = tmp_path / "traits.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=True))
    cache = tmp_path / "sumstats_cache.mat"
    bad_schema = tmp_path / "sumstats_bad_schema.mat"
    bad_lengths = tmp_path / "sumstats_bad_lengths.mat"
    bad_zero_shards = tmp_path / "sumstats_zero_shards.mat"
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        f"statgen.save_sumstats_cache(s, '{cache}'); "
        f"L = load('{cache}'); "
        "metadata = L.metadata; zvec = L.zvec; nvec = L.nvec; logpvec = L.logpvec; "
        "beta_vec = L.beta_vec; se_vec = L.se_vec; eaf_vec = L.eaf_vec; info_vec = L.info_vec; "
        "metadata.schema = 'sumstats_cache/bad'; "
        f"save('{bad_schema}', 'metadata', 'zvec', 'nvec', 'logpvec', 'beta_vec', 'se_vec', 'eaf_vec', 'info_vec'); "
        f"ok1 = 0; try; statgen.load_sumstats_cache('{bad_schema}'); catch; ok1 = 1; end; "
        f"L2 = load('{cache}'); "
        "metadata = L2.metadata; zvec = L2.zvec; nvec = L2.nvec; logpvec = L2.logpvec; "
        "beta_vec = L2.beta_vec; se_vec = L2.se_vec; eaf_vec = L2.eaf_vec; info_vec = L2.info_vec; "
        "metadata.shard_checksums = metadata.shard_checksums(1:end-1); "
        f"save('{bad_lengths}', 'metadata', 'zvec', 'nvec', 'logpvec', 'beta_vec', 'se_vec', 'eaf_vec', 'info_vec'); "
        f"ok2 = 0; try; statgen.load_sumstats_cache('{bad_lengths}'); catch; ok2 = 1; end; "
        f"ok3 = 0; try; statgen.load_sumstats_cache('{cache}', {{'2'}}); catch; ok3 = 1; end; "
        f"L3 = load('{cache}'); "
        "metadata = L3.metadata; zvec = L3.zvec; nvec = L3.nvec; logpvec = L3.logpvec; "
        "beta_vec = L3.beta_vec; se_vec = L3.se_vec; eaf_vec = L3.eaf_vec; info_vec = L3.info_vec; "
        "metadata.n_shards = 0; metadata.shard_labels = {}; metadata.shard_checksums = {}; metadata.shard_start0 = []; metadata.shard_stop0 = []; "
        f"save('{bad_zero_shards}', 'metadata', 'zvec', 'nvec', 'logpvec', 'beta_vec', 'se_vec', 'eaf_vec', 'info_vec'); "
        f"ok4 = 0; try; statgen.load_sumstats_cache('{bad_zero_shards}'); catch ME; ok4 = ~isempty(strfind(ME.message, 'n_shards must be positive')); end; "
        "fprintf('%d %d %d %d\\n', ok1, ok2, ok3, ok4);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_optional_fields_absent_use_empty_vector_sentinel(tmp_path):
    path = tmp_path / "traits_required_only.tsv.gz"
    _write_gz_tsv(path, _valid_sumstats_text(include_optional=False))
    cache = tmp_path / "sumstats_required_only_cache.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "fprintf('%d %d %d %d %d %d\\n', isempty(s.zvec), isempty(s.nvec), isempty(s.beta_vec), isempty(s.se_vec), isempty(s.eaf_vec), isempty(s.info_vec)); "
        f"statgen.save_sumstats_cache(s, '{cache}'); "
        f"s2 = statgen.load_sumstats_cache('{cache}'); "
        "fprintf('%d %d %d %d %d %d\\n', isempty(s2.zvec), isempty(s2.nvec), isempty(s2.beta_vec), isempty(s2.se_vec), isempty(s2.eaf_vec), isempty(s2.info_vec));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "1 1 1 1 1 1"
    assert lines[1] == "1 1 1 1 1 1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_create_sumstats_from_vectors():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "z = [2.5; 1.0; 1.8; -1.2; NaN; 3.0; 0.5; NaN]; "
        "n = [1000; 950; 1000; 1000; NaN; 500; 500; NaN]; "
        "p = [0.01; 0.5; 0; 0.2; 1.0; 0.003; 0.6; 0.9]; "
        "s = statgen.create_sumstats(ref, p, z, n); "
        "ok = zeros(1,8); "
        "ok(1) = (s.logpvec(1) == 2); "
        "ok(2) = abs(s.logpvec(2) - (-log10(0.5))) < 1e-12; "
        "ok(3) = isinf(s.logpvec(3)); "
        "ok(4) = abs(s.logpvec(4) - (-log10(0.2))) < 1e-12; "
        "ok(5) = (s.logpvec(5) == 0); "
        "ok(6) = abs(s.logpvec(6) - (-log10(0.003))) < 1e-12; "
        "ok(7) = abs(s.logpvec(7) - (-log10(0.6))) < 1e-12; "
        "ok(8) = abs(s.logpvec(8) - (-log10(0.9))) < 1e-12; "
        "fprintf('%d %d %d %d %d %d %d %d\\n', ok); "
        "fprintf('%d\\n', isempty(s.beta_vec) && isempty(s.se_vec) && isempty(s.eaf_vec) && isempty(s.info_vec));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "1 1 1 1 1 1 1 1"
    assert lines[1] == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_create_sumstats_rejects_inf_inputs():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "z = zeros(ref.num_snp, 1); "
        "n = zeros(ref.num_snp, 1); "
        "p = zeros(ref.num_snp, 1); p(1) = Inf; "
        "b = zeros(ref.num_snp, 1); b(1) = Inf; "
        "ok1 = 0; try; statgen.create_sumstats(ref, p, z, n); catch; ok1 = 1; end; "
        "p = zeros(ref.num_snp, 1); "
        "ok2 = 0; try; statgen.create_sumstats(ref, p, z, n, b); catch; ok2 = 1; end; "
        "fprintf('%d %d\\n', ok1, ok2);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "1 1"
