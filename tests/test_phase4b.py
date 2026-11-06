import numpy as np
import pytest
from scipy import sparse

from statgen.ld import _write_ld_npz_distribution, fast_prune, load_ld
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, matlab_data_lines, run_octave, skipif_no_octave


SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
LD_PY = FIXTURES_DIR / "ld/python"


def _expected_multiply(ld, values, chrX_sex=None):
    assert values.ndim == 1
    out = []
    for group in ld.shard_groups:
        if group[0].chr == "X":
            sex = ld.default_chrX_sex if chrX_sex is None else chrX_sex
            shard = {s.sex: s for s in group}[sex]
        else:
            shard = group[0]
        n = shard.num_snp
        part = values[:n]
        values = values[n:]
        out.append((shard.ld_r.toarray() ** 2) @ part)
    return np.concatenate(out)


def _synthetic_specs(reference):
    return [
        {
            "reference_shard": reference.shards[0],
            "a1freq": [0.10, 0.20, 0.30, 0.40, 0.50],
            "ld_pairs": np.array([[0, 1, 0.25], [2, 3, -0.50]], dtype=np.float64),
        },
        {
            "reference_shard": reference.shards[1],
            "sex": "female",
            "a1freq": [0.22, 0.33, 0.44],
            "ld_pairs": np.array([[0, 2, 0.70]], dtype=np.float64),
        },
    ]


def _float_line(line):
    return np.array([float(x) for x in line.split()], dtype=np.float64)


def test_ld_panel_a1freq_chrx_selection_and_select_shards():
    ld = load_ld(LD_PY, default_chrX_sex="male")

    np.testing.assert_allclose(
        ld.a1freq("female"),
        [0.30, 0.40, 0.20, 0.35, 0.15, 0.28, 0.32, 0.22],
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        ld.a1freq(),
        [0.30, 0.40, 0.20, 0.35, 0.15, 0.22, 0.28, 0.18],
        rtol=1e-6,
    )
    with pytest.raises(ValueError, match="chrX sex"):
        ld.a1freq("unknown")

    chr1 = ld.select_shards(["1"])
    assert chr1.default_chrX_sex == "male"
    assert chr1.num_snp == 5
    assert [s.label for s in chr1.shards] == ["1"]
    np.testing.assert_allclose(chr1.a1freq("unknown"), [0.30, 0.40, 0.20, 0.35, 0.15])

    chrx = ld.select_shards(["X"])
    assert chrx.default_chrX_sex == "male"
    assert chrx.shards[0].sex == "male"


def test_ld_panel_multiply_r2_vector_matrix_sparse_and_dtype():
    reference = load_reference(SHARDED_REF)
    ld = load_ld(LD_PY)
    vec = np.arange(1, reference.num_snp + 1, dtype=np.float64)

    np.testing.assert_allclose(ld.multiply_r2(vec), _expected_multiply(ld, vec), rtol=1e-7)
    np.testing.assert_allclose(
        ld.multiply_r2(vec, chrX_sex="male"),
        _expected_multiply(ld, vec, chrX_sex="male"),
        rtol=1e-7,
    )

    mat32 = np.column_stack([vec, vec * 2]).astype(np.float32)
    got32 = ld.multiply_r2(mat32)
    assert got32.shape == mat32.shape
    assert got32.dtype == np.float32
    np.testing.assert_allclose(got32[:, 0], ld.multiply_r2(vec).astype(np.float32), rtol=1e-6)

    mat64 = np.column_stack([vec, np.array([2, -1, 0, 3, 5, 7, -2, 4], dtype=np.float64)])
    got64 = ld.multiply_r2(mat64)
    assert got64.shape == mat64.shape
    assert got64.dtype == np.float64
    np.testing.assert_allclose(got64[:, 0], _expected_multiply(ld, mat64[:, 0]), rtol=1e-7)
    np.testing.assert_allclose(got64[:, 1], _expected_multiply(ld, mat64[:, 1]), rtol=1e-7)

    sparse_in = sparse.csr_matrix(mat32)
    sparse_out = ld.multiply_r2(sparse_in)
    assert sparse.issparse(sparse_out)
    assert sparse_out.dtype == np.float32
    np.testing.assert_allclose(sparse_out.toarray(), got32, rtol=1e-6)

    chr1 = ld.select_shards(["1"])
    np.testing.assert_allclose(chr1.multiply_r2(vec[:5], chrX_sex="male"), chr1.multiply_r2(vec[:5]))
    np.testing.assert_allclose(chr1.multiply_r2(vec[:5], chrX_sex="unknown"), chr1.multiply_r2(vec[:5]))

    with pytest.raises(ValueError, match="vector length mismatch"):
        ld.multiply_r2(vec[:-1])
    with pytest.raises(ValueError, match="vector or 2D matrix"):
        ld.multiply_r2(np.zeros((reference.num_snp, 1, 1), dtype=np.float32))


def test_fast_prune_stable_order_and_chrx_override():
    reference = load_reference(SHARDED_REF)
    ld = load_ld(LD_PY)
    logp = np.array([9, 9, 7, 1, 2, 6, 5, 4], dtype=np.float64)

    np.testing.assert_allclose(
        fast_prune(logp, ld, r2_threshold=0.35),
        [9, np.nan, 7, np.nan, 2, 6, np.nan, 4],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        fast_prune(logp, ld, r2_threshold=0.35, chrX_sex="male"),
        [9, np.nan, 7, np.nan, 2, 6, 5, 4],
        equal_nan=True,
    )

    out32 = fast_prune(logp.astype(np.float32), ld, r2_threshold=0.35)
    assert out32.dtype == np.float32
    out_int = fast_prune(logp.astype(np.int64), ld, r2_threshold=0.35)
    assert out_int.dtype == np.float64

    chr1 = ld.select_shards(["1"])
    np.testing.assert_allclose(
        fast_prune(logp[:5], chr1, r2_threshold=0.35, chrX_sex="male"),
        fast_prune(logp[:5], chr1, r2_threshold=0.35),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        fast_prune(logp[:5], chr1, r2_threshold=0.35, chrX_sex="unknown"),
        fast_prune(logp[:5], chr1, r2_threshold=0.35),
        equal_nan=True,
    )

    all_nan = np.full(reference.num_snp, np.nan, dtype=np.float64)
    np.testing.assert_allclose(fast_prune(all_nan, ld), all_nan, equal_nan=True)

    with pytest.raises(ValueError, match="vector length mismatch"):
        fast_prune(logp[:-1], ld)


def test_phase5a_generated_artifacts_support_phase4b_operations(tmp_path):
    reference = load_reference(SHARDED_REF)
    root = tmp_path / "ld_python"
    _write_ld_npz_distribution(root, _synthetic_specs(reference))
    ld = load_ld(root)

    np.testing.assert_allclose(ld.a1freq(), [0.10, 0.20, 0.30, 0.40, 0.50, 0.22, 0.33, 0.44])
    vec = np.arange(1, reference.num_snp + 1, dtype=np.float64)
    np.testing.assert_allclose(ld.multiply_r2(vec), _expected_multiply(ld, vec))
    np.testing.assert_allclose(
        fast_prune(np.array([5, 4, 3, 2, 1, 6, 5, 4], dtype=float), ld, 0.4),
        [5, 4, 3, 2, 1, 6, 5, np.nan],
        equal_nan=True,
    )


@pytest.mark.octave
@skipif_no_octave
def test_octave_ld_panel_operations_match_fixture_values():
    script = (
        "statgen.set_verbosity('quiet'); "
        f"ld = statgen.load_ld('{FIXTURES_DIR / 'ld/matlab'}'); "
        "M = (1:ld.num_snp)'; "
        "R = ld.multiply_r2(M); "
        "Rm = ld.multiply_r2(M, 'male'); "
        "P = statgen.fast_prune([9;9;7;1;2;6;5;4], ld, 0.35); "
        "Pm = statgen.fast_prune([9;9;7;1;2;6;5;4], ld, 0.35, 'male'); "
        "sub = ld.select_shards({'1'}); "
        "sf = sub.a1freq('unknown'); "
        "fprintf('%.2f ', ld.a1freq()); fprintf('\\n'); "
        "fprintf('%.2f ', R); fprintf('\\n'); "
        "fprintf('%.2f ', Rm); fprintf('\\n'); "
        "fprintf('%.2f ', P); fprintf('\\n'); "
        "fprintf('%.2f ', Pm); fprintf('\\n'); "
        "fprintf('%d,%d,%.2f\\n', sub.num_snp, numel(sub.shards), sf(1));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    np.testing.assert_allclose(
        _float_line(lines[0]),
        [0.30, 0.40, 0.20, 0.35, 0.15, 0.28, 0.32, 0.22],
        atol=0.005,
    )
    np.testing.assert_allclose(
        _float_line(lines[1]),
        [2.89, 3.81, 3.09, 6.95, 6.96, 8.96, 10.515, 8.86],
        atol=0.01,
    )
    np.testing.assert_allclose(
        _float_line(lines[2]),
        [2.89, 3.81, 3.09, 6.95, 6.96, 8.12, 10.43, 9.42],
        atol=0.01,
    )
    np.testing.assert_allclose(
        _float_line(lines[3]),
        [9, np.nan, 7, np.nan, 2, 6, np.nan, 4],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        _float_line(lines[4]),
        [9, np.nan, 7, np.nan, 2, 6, 5, 4],
        equal_nan=True,
    )
    assert lines[5] == "5,1,0.30"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_ld_can_drop_raw_ld_r_but_keep_operations():
    script = (
        "statgen.set_verbosity('quiet'); "
        f"ld = statgen.load_ld('{FIXTURES_DIR / 'ld/matlab'}', 'female', false); "
        "M = (1:ld.num_snp)'; "
        "R = ld.multiply_r2(M); "
        "P = statgen.fast_prune([9;9;7;1;2;6;5;4], ld, 0.35); "
        "fprintf('%d,%d\\n', isempty(ld.shard_groups{1}{1}.ld_r), isempty(ld.shard_groups{1}{1}.ld_r2)); "
        "fprintf('%.2f %.2f %.2f\\n', R(1), R(6), R(7)); "
        "fprintf('%.2f %.2f %.2f\\n', P(1), P(6), P(8));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "1,0"
    np.testing.assert_allclose(_float_line(lines[1]), [2.89, 8.96, 10.515], atol=0.01)
    np.testing.assert_allclose(_float_line(lines[2]), [9, 6, 4], atol=0.01)
