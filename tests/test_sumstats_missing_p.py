"""Missing-p semantics at the source and factory boundaries in all runtimes."""

import json
import warnings

import numpy as np
import pytest

from statgen import get_verbosity, set_verbosity
from statgen.reference import load_reference
from statgen.sumstats import create_sumstats, load_sumstats, load_sumstats_cache
from tests.conftest import (
    FIXTURES_DIR, R_PACKAGE_DIR, matlab_output_text, run_octave, run_rscript,
    skipif_no_octave, skipif_no_rscript,
)

REFERENCE = FIXTURES_DIR / "reference/sharded/1.bim"
FIELDS = ("zvec", "nvec", "beta_vec", "se_vec", "eaf_vec", "info_vec")
HEADER = "chr\tbp\ta1\ta2\tp\tz\tn\tbeta\tse\teaf\tinfo\n"


def _r_script(expr):
    sources = " ".join(
        f"source({json.dumps(str(path))});"
        for path in sorted((R_PACKAGE_DIR / "R").glob("*.R"))
    )
    return sources + expr


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "missing.tsv"
    path.write_text(
        HEADER
        + "1\t100\tA\tG\tNaN\t9\t9\t9\t9\t0.9\t9\n"
        + "1\t100\tA\tG\tNA\t9\t9\t9\t9\t0.9\t9\n"
        + "1\t100\tG\tA\t0.01\t2\t100\t0.2\t0.1\t0.3\t0.9\n"
        + "1\t200\tC\tT\tNA\tNaN\tNaN\t9\t9\t0.9\t9\n"
        + "1\t300\tA\tC\t\t9\t9\t9\t9\t0.9\t9\n"
        + "1\t400\tG\tA\t0\t4\t100\t0.4\t0.1\t0.4\t0.9\n"
    )
    return path


def test_missing_p_before_matching_and_reporting(source, tmp_path, capsys):
    ref = load_reference(REFERENCE)
    previous = get_verbosity()
    try:
        set_verbosity("info")
        with warnings.catch_warnings(record=True) as caught:
            stats = load_sumstats(source, ref)
        assert not caught  # Missing z/n only occur on dropped rows.
        assert "dropped 4 row(s) with missing p" in capsys.readouterr().err
        np.testing.assert_array_equal(stats.is_present, [True, False, False, True, False])
        assert stats.logpvec[0] == 2 and np.isposinf(stats.logpvec[3])
        assert stats.zvec[0] == -2 and stats.beta_vec[0] == -0.2
        assert stats.eaf_vec[0] == 0.7
        for field in FIELDS:
            assert np.isnan(getattr(stats, field)[[1, 2, 4]]).all()
        cache = tmp_path / "stats.npz"
        stats.save_cache(cache)
        restored = load_sumstats_cache(cache)
        for field in (*FIELDS, "logpvec", "is_present"):
            np.testing.assert_equal(getattr(restored, field), getattr(stats, field))
        set_verbosity("quiet")
        load_sumstats(source, ref)
        assert capsys.readouterr().err == ""
    finally:
        set_verbosity(previous)


@pytest.mark.parametrize("present_columns", [True, False])
def test_all_missing_p_preserves_optional_presence(tmp_path, present_columns):
    path = tmp_path / "all_missing.tsv"
    path.write_text(
        HEADER + "1\t100\tA\tG\tNA\t1\t2\t3\t4\t0.5\t6\n"
        if present_columns else "chr\tbp\ta1\ta2\tp\n1\t100\tA\tG\tNA\n"
    )
    with pytest.warns(RuntimeWarning) as caught:
        stats = load_sumstats(path, load_reference(REFERENCE))
    assert any("no matching source rows" in str(w.message) for w in caught)
    assert not stats.is_present.any()
    for field in FIELDS:
        value = getattr(stats, field)
        if present_columns:
            assert value.shape == (5,) and np.isnan(value).all()
        else:
            assert value is None


def test_factory_missing_p_masks_copies_and_roundtrips(tmp_path):
    ref = load_reference(REFERENCE)
    p = np.array([np.nan, 0, 1, 0.1, np.nan])
    original = np.arange(5, dtype=float)
    stats = create_sumstats(ref, p, **dict.fromkeys(("z", "n", "beta", "se", "eaf", "info"), original))
    np.testing.assert_equal(original, np.arange(5))
    np.testing.assert_equal(p, [np.nan, 0, 1, 0.1, np.nan])
    np.testing.assert_equal(stats.is_present, [False, True, True, True, False])
    for field in FIELDS:
        np.testing.assert_equal(getattr(stats, field), [np.nan, 1, 2, 3, np.nan])
    cache = tmp_path / "created.npz"
    stats.save_cache(cache)
    restored = load_sumstats_cache(cache)
    np.testing.assert_equal(restored.logpvec, stats.logpvec)
    for field in FIELDS:
        np.testing.assert_equal(getattr(restored, field), getattr(stats, field))
    assert create_sumstats(ref, p).zvec is None


@pytest.mark.octave
@skipif_no_octave
def test_octave_missing_p_load_factory_and_cache(source, tmp_path):
    cache = tmp_path / "stats.mat"
    result = run_octave(
        f"ref = statgen.load_reference('{REFERENCE}'); "
        f"s = statgen.load_sumstats('{source}', ref); "
        "assert(isequal(s.is_present, logical([1;0;0;1;0]))); "
        "assert(s.logpvec(1) == 2 && isinf(s.logpvec(4))); "
        "assert(s.zvec(1) == -2 && s.beta_vec(1) == -0.2 && s.eaf_vec(1) == 0.7); "
        "fields = {'zvec','nvec','beta_vec','se_vec','eaf_vec','info_vec'}; "
        "for k = 1:numel(fields); v = s.(fields{k}); assert(all(isnan(v([2 3 5])))); end; "
        f"s.save_cache('{cache}'); t = statgen.load_sumstats_cache('{cache}'); "
        "for k = 1:numel(fields); assert(isequaln(t.(fields{k}), s.(fields{k}))); end; "
        "assert(isequaln(t.logpvec, s.logpvec)); "
        "p = [NaN;0;1;0.1;NaN]; v = (1:5)'; "
        "s = statgen.create_sumstats(ref, p, v, v, v, v, v, v); "
        "assert(isequal(s.is_present, logical([0;1;1;1;0]))); "
        "for k = 1:numel(fields); assert(isequaln(s.(fields{k}), [NaN;2;3;4;NaN])); end; "
        f"s.save_cache('{cache}'); t = statgen.load_sumstats_cache('{cache}'); "
        "for k = 1:numel(fields); assert(isequaln(t.(fields{k}), s.(fields{k}))); end; "
        "assert(isequaln(t.logpvec, s.logpvec)); assert(isequal(v, (1:5)')); "
        "s = statgen.create_sumstats(ref, p); assert(isempty(s.zvec)); "
        "statgen.set_verbosity('quiet'); "
        f"statgen.load_sumstats('{source}', ref);"
    )
    assert result.returncode == 0, result.stderr
    assert (result.stdout + result.stderr).count("dropped 4 row(s) with missing p") == 1
    assert "missing values among present" not in result.stdout + result.stderr


@pytest.mark.octave
@skipif_no_octave
def test_octave_all_missing_p_preserves_optional_presence(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text(HEADER + "1\t100\tA\tG\tNA\t1\t2\t3\t4\t0.5\t6\n")
    result = run_octave(
        f"ref = statgen.load_reference('{REFERENCE}'); "
        f"s = statgen.load_sumstats('{path}', ref); "
        "assert(~any(s.is_present)); "
        "fields = {'zvec','nvec','beta_vec','se_vec','eaf_vec','info_vec'}; "
        "for k = 1:numel(fields); v = s.(fields{k}); assert(numel(v) == 5 && all(isnan(v))); end;"
    )
    assert result.returncode == 0, result.stderr
    output = matlab_output_text(result.stdout + result.stderr)
    assert "no matching source rows" in output
    assert "is absent" not in output


@pytest.mark.r
@skipif_no_rscript
def test_r_missing_p_load_and_cache(source, tmp_path):
    cache = tmp_path / "stats.rds"
    result = run_rscript(_r_script(
        f"ref <- load_reference('{REFERENCE}'); s <- load_sumstats('{source}', ref); "
        "stopifnot(identical(is_present(s), c(TRUE,FALSE,FALSE,TRUE,FALSE)), "
        "logpvec(s)[1] == 2, is.infinite(logpvec(s)[4]), zvec(s)[1] == -2, "
        "beta_vec(s)[1] == -0.2, eaf_vec(s)[1] == 0.7); "
        "fields <- list(zvec,nvec,beta_vec,se_vec,eaf_vec,info_vec); "
        "for (f in fields) stopifnot(all(is.nan(f(s)[c(2,3,5)]))); "
        f"save_sumstats_cache(s, '{cache}'); t <- load_sumstats_cache('{cache}'); "
        "for (f in c(fields, list(logpvec))) stopifnot(identical(f(t), f(s))); "
        "set_verbosity('quiet'); "
        f"load_sumstats('{source}', ref);"
    ))
    assert result.returncode == 0, result.stderr
    assert result.stderr.count("dropped 4 row(s) with missing p") == 1
    assert "missing values among present" not in result.stderr


@pytest.mark.parametrize("field,token", [
    ("chr", "chr1"), ("bp", "NA"), ("bp", "1.5"),
    ("a1", "a"), ("a2", "A"),
    *[(field, "invalid") for field in ("bp", "p", "z", "n", "beta", "se", "eaf", "info")],
])
def test_missing_p_does_not_bypass_source_validation(tmp_path, field, token):
    columns = HEADER.strip().split("\t")
    row = ["1", "100", "A", "G", "NA", "1", "2", "3", "4", "0.5", "6"]
    row[columns.index(field)] = token
    path = tmp_path / "bad.tsv"
    path.write_text(HEADER + "\t".join(row) + "\n")
    with pytest.raises(ValueError):
        load_sumstats(path, load_reference(REFERENCE))


@pytest.mark.parametrize("runtime", [
    pytest.param("octave", marks=[pytest.mark.octave, skipif_no_octave]),
    pytest.param("r", marks=[pytest.mark.r, skipif_no_rscript]),
])
def test_other_runtimes_validate_before_dropping(tmp_path, runtime):
    paths = []
    for field, token in (
        [("chr", "chr1"), ("bp", "NA"), ("bp", "1.5"), ("a1", "a"), ("a2", "A")]
        + [(field, "invalid") for field in ("bp", "p", "z", "n", "beta", "se", "eaf", "info")]
    ):
        row = ["1", "100", "A", "G", "NA", "1", "2", "3", "4", "0.5", "6"]
        row[HEADER.strip().split("\t").index(field)] = token
        path = tmp_path / f"bad_{len(paths)}.tsv"
        path.write_text(HEADER + "\t".join(row) + "\n")
        paths.append(str(path))
    if runtime == "octave":
        path_list = "{" + ",".join(f"'{path}'" for path in paths) + "}"
        result = run_octave(
            f"ref = statgen.load_reference('{REFERENCE}'); paths = {path_list}; "
            "for k = 1:numel(paths); failed = false; "
            "try; statgen.load_sumstats(paths{k}, ref); catch; failed = true; end; "
            "assert(failed, paths{k}); end;"
        )
    else:
        path_list = "c(" + ",".join(json.dumps(path) for path in paths) + ")"
        result = run_rscript(_r_script(
            f"ref <- load_reference('{REFERENCE}'); paths <- {path_list}; "
            "for (path in paths) stopifnot(inherits(try(load_sumstats(path, ref), silent=TRUE), 'try-error'));"
        ))
    assert result.returncode == 0, result.stderr
