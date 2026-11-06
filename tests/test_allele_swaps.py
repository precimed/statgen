"""Stage 1: shared independent oracle, boundary failures, and cache migration."""
import json
import shutil

import numpy as np
import pandas as pd
import pytest

from statgen import get_verbosity, set_verbosity
from statgen._utils import allele_hash64
from statgen._variant_match import match_shard_numeric, numeric_variant_keys
from statgen.genotype import load_genotype, load_genotype_cache
from statgen.reference import load_reference, load_reference_cache
from statgen.sumstats import load_sumstats, load_sumstats_cache
from tests.conftest import FIXTURES_DIR

ROOT = FIXTURES_DIR / "allele_swaps"
EXPECTED = json.loads((ROOT / "expected.json").read_text())


def _keys(bp, a1, a2):
    return numeric_variant_keys(bp, allele_hash64(a1), allele_hash64(a2))


def _as_float(values):
    return np.array(values, dtype=float)


@pytest.fixture(autouse=True)
def restore_verbosity():
    previous = get_verbosity()
    set_verbosity("info")
    yield
    set_verbosity(previous)


@pytest.mark.parametrize("enabled", [False, True])
def test_shared_swap_oracle(enabled, tmp_path, capsys):
    expected = EXPECTED["options"][str(enabled).lower()]
    reference = load_reference(ROOT / "reference.bim")
    with pytest.warns(RuntimeWarning, match="missing"):
        stats = load_sumstats(ROOT / "sumstats.tsv.gz", reference, allow_indel_swaps=enabled)
    stats_reporting = capsys.readouterr().err
    genotype = load_genotype(ROOT / "genotype/@", reference, allow_indel_swaps=enabled)
    genotype_reporting = capsys.readouterr().err
    np.testing.assert_array_equal(genotype.source_row0, expected["source_row0"])
    np.testing.assert_array_equal(genotype.is_allele_swapped, expected["is_allele_swapped"])
    present = np.array(expected["source_row0"]) >= 0
    np.testing.assert_array_equal(genotype.is_present, present)
    np.testing.assert_array_equal(stats.is_present, present)
    for field, accessor in [("z", "zvec"), ("n", "nvec"), ("beta", "beta_vec"),
                            ("se", "se_vec"), ("eaf", "eaf_vec"), ("info", "info_vec")]:
        np.testing.assert_allclose(getattr(stats, accessor),
                                   _as_float(expected["aligned_sumstats"][field]), equal_nan=True)
    p = _as_float(expected["aligned_sumstats"]["p"])
    np.testing.assert_allclose(stats.logpvec, -np.log10(p), equal_nan=True)
    # Verify the portable TSV mapping independently through the numeric matcher.
    source = pd.read_csv(ROOT / "sumstats.tsv.gz", sep="\t", dtype={"chr": str})
    offsets = {"1": 0, "X": 13}
    for shard in reference.shards:
        src_indices = np.flatnonzero(source["chr"] == shard.label)
        result = match_shard_numeric(
            _keys(shard.bp, shard.a1, shard.a2),
            _keys(source["bp"][src_indices], source["a1"][src_indices], source["a2"][src_indices]),
            shard.label, allow_indel_swaps=enabled,
        )
        rows = result.source_rows.copy()
        found = rows >= 0
        rows[found] = src_indices[rows[found]]
        start = offsets[shard.label]
        np.testing.assert_array_equal(rows, expected["sumstats_source_row0"][start:start + shard.num_snp])
        counts = expected["reporting"][shard.label]
        assert np.count_nonzero(result.is_allele_swapped) == counts["swaps"]
        assert result.num_strand_ambiguous == counts["strand_ambiguous"]
        assert result.num_disabled_indels == counts["disabled_indels"]
        for output in (stats_reporting, genotype_reporting):
            assert f"shard {shard.label}: {counts['swaps']} successful allele swap(s)" in output
            if counts["strand_ambiguous"]:
                assert f"shard {shard.label}: {counts['strand_ambiguous']} strand-ambiguous A/T or C/G swap(s)" in output
            if counts["disabled_indels"]:
                assert f"shard {shard.label}: {counts['disabled_indels']} indel row(s)" in output
    assert "0 strand-ambiguous" not in stats_reporting
    assert "0 indel row(s)" not in stats_reporting
    assert str(ROOT / "sumstats.tsv.gz") in stats_reporting
    assert str(ROOT / "genotype/1.bim") in genotype_reporting

    # Repeated out-of-order SNP requests exercise flips and chrX sample expansion.
    requested = np.r_[np.flatnonzero(present)[::-1], [13, 1, 13]]
    raw = _as_float(expected["raw_counts"])[:, requested]
    scaled = _as_float(expected["scaled_counts"])[:, requested]
    int_expected = np.where(np.isnan(raw), -1, raw).astype(np.int8)
    np.testing.assert_array_equal(genotype.fetch_genotypes_int8(requested), int_expected)
    np.testing.assert_allclose(genotype.fetch_genotypes(requested), raw, equal_nan=True)
    np.testing.assert_allclose(genotype.fetch_genotypes(requested, haploid_mode="ploidy_scaled"),
                               scaled, equal_nan=True)
    np.testing.assert_array_equal(genotype.fetch_genotypes_int8(requested), int_expected)
    with pytest.raises(ValueError, match="not present"):
        genotype.fetch_genotypes_int8([8])
    np.testing.assert_array_equal(genotype.is_subject_present("X"), [False, True, True, False])
    assert reference.is_object_compatible(genotype)
    cache = tmp_path / "genotype.npz"
    genotype.save_cache(cache)
    with np.load(cache) as data:
        assert json.loads(bytes(data["_meta"]))["schema"] == "genotype_cache/0.2"
        assert data["is_allele_swapped"].dtype == bool
    loaded = load_genotype_cache(cache)
    np.testing.assert_array_equal(loaded.is_allele_swapped, genotype.is_allele_swapped)
    np.testing.assert_array_equal(loaded.fetch_genotypes_int8(requested), int_expected)
    np.testing.assert_allclose(loaded.fetch_genotypes(requested, haploid_mode="ploidy_scaled"),
                               scaled, equal_nan=True)
    # Remove source metadata and BEDs from the cached paths; overrides still use
    # stored source orientation, mappings, ploidy, and sample reconciliation.
    arrays = _cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]))
    meta["bed_paths"] = [str(tmp_path / f"gone{label}.bed") for label in ("1", "X")]
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    np.savez(cache, **arrays)
    loaded = load_genotype_cache(cache)
    np.testing.assert_array_equal(
        loaded.fetch_genotypes_int8(requested, bed_path=str(ROOT / "genotype/@.bed")), int_expected)
    for subset in (loaded.select_shards(["X"]), load_genotype_cache(cache, shards=["X"])):
        np.testing.assert_array_equal(subset.is_allele_swapped, genotype.is_allele_swapped[13:])
        req = [0, 0, 1] + ([2] if enabled else [])
        np.testing.assert_allclose(
            subset.fetch_genotypes(req, bed_path=str(ROOT / "genotype/@.bed"), haploid_mode="ploidy_scaled"),
            _as_float(expected["scaled_counts"])[:, np.array(req) + 13], equal_nan=True)
    stats_cache = tmp_path / "stats.npz"
    stats.save_cache(stats_cache)
    with np.load(stats_cache) as data:
        assert json.loads(bytes(data["_meta"]))["schema"] == "sumstats_cache/0.1"
        assert "is_allele_swapped" not in data.files
    with pytest.warns(RuntimeWarning, match="missing"):
        roundtrip = load_sumstats_cache(stats_cache)
    np.testing.assert_allclose(roundtrip.zvec, stats.zvec, equal_nan=True)
    np.testing.assert_allclose(roundtrip.eaf_vec, stats.eaf_vec, equal_nan=True)


def _cache_arrays(path):
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


@pytest.mark.parametrize("corruption", ["missing", "shape", "dtype", "absent"])
def test_swap_cache_validation(corruption, tmp_path):
    ref = load_reference(ROOT / "reference.bim")
    cache = tmp_path / "cache.npz"
    load_genotype(ROOT / "genotype/@", ref).save_cache(cache)
    arrays = _cache_arrays(cache)
    if corruption == "missing":
        del arrays["is_allele_swapped"]
    elif corruption == "shape":
        arrays["is_allele_swapped"] = arrays["is_allele_swapped"][:, None]
    elif corruption == "dtype":
        arrays["is_allele_swapped"] = arrays["is_allele_swapped"].astype(float)
    else:
        arrays["is_allele_swapped"][8] = True
    np.savez(cache, **arrays)
    with pytest.raises(ValueError, match="is_allele_swapped"):
        load_genotype_cache(cache)


def test_legacy_cache_has_false_swaps_without_source_parsing(monkeypatch, tmp_path):
    import statgen.genotype as module

    def fail(*args, **kwargs):
        raise AssertionError("legacy cache must not parse or rematch sources")

    monkeypatch.setattr(module, "_parse_bim", fail)
    monkeypatch.setattr(module, "_parse_fam", fail)
    monkeypatch.setattr(module, "_parse_ploidy", fail)
    monkeypatch.setattr(module, "match_shard_numeric", fail)
    cache = FIXTURES_DIR / "cache/genotype_python_genotype_cache_0_1.npz"
    panel = load_genotype_cache(cache)
    assert not panel.is_allele_swapped.any()
    np.testing.assert_array_equal(panel.fetch_genotypes_int8([0, 5, 0]), np.full((4, 3), 2, np.int8))
    migrated = tmp_path / "migrated.npz"
    panel.save_cache(migrated)
    np.testing.assert_array_equal(load_genotype_cache(migrated).source_row0, panel.source_row0)
    assert not load_genotype_cache(migrated).is_allele_swapped.any()


@pytest.mark.parametrize("alleles", [("A", "G"), ("A", "T"), ("C", "G")])
def test_reference_rejects_reciprocal_snvs(alleles, tmp_path):
    a, b = alleles
    path = tmp_path / "reference.bim"
    path.write_text("1\tlong\t0\t50\tAC\tGT\n1\tsnv\t0\t75\tA\tC\n"
                    f"1\ta\t0\t100\t{a}\t{b}\n1\tb\t0\t100\t{b}\t{a}\n")
    with pytest.raises(ValueError, match="reference shard 1: reciprocal SNV keys are not allowed at bp=100"):
        load_reference(path)


@pytest.mark.parametrize("enabled", [False, True])
def test_matcher_competing_candidates_and_exact_precedence(enabled):
    ref = _keys([100], ["A"], ["G"])
    for a1, a2 in [(["A", "A"], ["G", "G"]), (["G", "G"], ["A", "A"])]:
        with pytest.raises(ValueError, match="Ambiguous duplicate"):
            match_shard_numeric(ref, _keys([100, 100], a1, a2), "1", allow_indel_swaps=enabled)
    # Duplicated reversed rows cannot make a retained exact match ambiguous.
    result = match_shard_numeric(ref, _keys([100] * 3, ["G", "A", "G"], ["A", "G", "A"]), "1")
    np.testing.assert_array_equal(result.source_rows, [1])
    assert not result.is_allele_swapped.any()
    # Duplicates outside the reference slot need no validation.
    result = match_shard_numeric(ref, _keys([200, 200], ["A", "A"], ["G", "G"]), "1")
    np.testing.assert_array_equal(result.source_rows, [-1])


def test_disabled_indel_count_counts_unused_rows_only():
    ref = _keys([100, 200, 200, 300, 400],
                ["A", "A", "AT", "A", "AC"], ["AG", "AT", "A", "G", "GT"])
    src = _keys([100, 100, 200, 300, 400, 999],
                ["AG", "AG", "A", "G", "GT", "AT"], ["A", "A", "AT", "A", "AC", "A"])
    disabled = match_shard_numeric(ref, src, "1")
    assert disabled.num_disabled_indels == 2  # Rows, not distinct reversed keys.
    np.testing.assert_array_equal(disabled.source_rows, [-1, 2, -1, 3, -1])
    np.testing.assert_array_equal(disabled.is_allele_swapped, [False, False, False, True, False])
    with pytest.raises(ValueError, match="Ambiguous duplicate"):
        match_shard_numeric(ref, src, "1", allow_indel_swaps=True)


def test_sumstats_unmatched_duplicates_and_absent_optional_fields(tmp_path, capsys):
    reference = load_reference(ROOT / "reference.bim", shards=["1"])
    path = tmp_path / "source.tsv"
    path.write_text("chr\tbp\ta1\ta2\tp\n1\t200\tG\tA\t0\n"
                    "1\t100\tA\tG\t0.1\n1\t100\tG\tA\t0.2\n"
                    "1\t100\tG\tA\t0.3\n1\t9999\tA\tC\t0.1\n1\t9999\tA\tC\t0.1\n")
    set_verbosity("quiet")
    with pytest.warns(RuntimeWarning) as recorded:
        stats = load_sumstats(path, reference)
    assert len(recorded) == 2  # Unrelated absent-Z/N warnings remain at quiet.
    assert stats.zvec is None and stats.eaf_vec is None
    assert stats.is_present[1] and np.isinf(stats.logpvec[1])
    assert capsys.readouterr().err == ""
    # Source validation is not bypassed for unmatched duplicate rows.
    path.write_text(path.read_text().replace("9999\tA\tC\t0.1", "9999\tA\tC\t1.1"))
    with pytest.raises(ValueError, match="p must be finite"):
        load_sumstats(path, reference)


def test_quiet_suppresses_swap_reports_only(capsys):
    ref = load_reference(ROOT / "reference.bim")
    set_verbosity("quiet")
    with pytest.warns(RuntimeWarning, match="missing"):
        load_sumstats(ROOT / "sumstats.tsv.gz", ref)
    load_genotype(ROOT / "genotype/@", ref)
    assert capsys.readouterr().err == ""


def test_cache_reference_matching_does_not_decode_or_revalidate(monkeypatch, tmp_path):
    from statgen.reference import ReferenceShard

    ref = load_reference(ROOT / "reference.bim")
    cache = tmp_path / "reference.npz"
    ref.save_cache(cache)
    cached = load_reference_cache(cache)

    def fail(*args):
        raise AssertionError("matching must trust numeric cached reference payloads")

    monkeypatch.setattr(ReferenceShard, "a1", property(fail))
    monkeypatch.setattr(ReferenceShard, "a2", property(fail))
    monkeypatch.setattr("statgen.reference._validate_reference_keys", fail)
    with pytest.warns(RuntimeWarning, match="missing"):
        stats = load_sumstats(ROOT / "sumstats.tsv.gz", cached)
    assert stats.zvec[1] == -3
    assert load_genotype(ROOT / "genotype/@", cached).is_allele_swapped[1]


def test_matcher_preserves_low_hash_bits():
    # These values coincide when rounded to double, but are distinct uint64s.
    values = np.array([2**60, 2**60 + 1], dtype=np.uint64)
    ref = numeric_variant_keys([100, 100], values, [1000, 1000])
    src = numeric_variant_keys([100, 100], values[::-1], [1000, 1000])
    result = match_shard_numeric(ref, src, "1")
    np.testing.assert_array_equal(result.source_rows, [1, 0])
    assert not result.is_allele_swapped.any()


def test_nonsharded_genotype_keeps_physical_rows(tmp_path):
    # Reuse the shared autosomal BED and prepend an ignored Y row. The physical
    # row and sidecar offsets must survive projection onto the reference shard.
    base = tmp_path / "source"
    source = ROOT / "genotype/1"
    shutil.copyfile(source.with_suffix(".fam"), base.with_suffix(".fam"))
    base.with_suffix(".bim").write_text("Y\tignored\t0\t1\tA\tG\n" + source.with_suffix(".bim").read_text())
    bed = source.with_suffix(".bed").read_bytes()
    base.with_suffix(".bed").write_bytes(bed[:3] + b"\x00" + bed[3:])
    base.with_suffix(".ploidy").write_text("0\t0\n" + "2\t2\n" * 14)
    ref = load_reference(ROOT / "reference.bim", shards=["1"])
    panel = load_genotype(base, ref)
    assert panel.source_row0[1] == 3
    np.testing.assert_array_equal(panel.fetch_genotypes_int8([1])[:, 0], [2, 0, 1, -1])
    np.testing.assert_array_equal(panel.ploidy_male[panel.is_present], 2)


def test_sumstats_alignment_is_independent_of_source_order(tmp_path):
    reference = load_reference(ROOT / "reference.bim")
    source = pd.read_csv(ROOT / "sumstats.tsv.gz", sep="\t", dtype={"chr": str})
    path = tmp_path / "permuted.tsv.gz"
    source.sample(frac=1, random_state=17).to_csv(path, sep="\t", index=False, na_rep="NA")
    with pytest.warns(RuntimeWarning, match="missing"):
        original = load_sumstats(ROOT / "sumstats.tsv.gz", reference, allow_indel_swaps=True)
    with pytest.warns(RuntimeWarning, match="missing"):
        permuted = load_sumstats(path, reference, allow_indel_swaps=True)
    for accessor in ("logpvec", "zvec", "nvec", "beta_vec", "se_vec", "eaf_vec", "info_vec"):
        np.testing.assert_allclose(getattr(permuted, accessor), getattr(original, accessor), equal_nan=True)


def test_swapped_missing_statistics_remain_missing(tmp_path):
    reference = load_reference(ROOT / "reference.bim", shards=["1"])
    path = tmp_path / "missing.tsv"
    path.write_text("chr\tbp\ta1\ta2\tp\tz\tn\tbeta\tse\teaf\tinfo\n"
                    "1\t200\tG\tA\t0.1\tNA\tNA\tNA\tNA\tNA\tNA\n")
    with pytest.warns(RuntimeWarning, match="missing"):
        stats = load_sumstats(path, reference)
    assert stats.is_present[1]
    assert stats.logpvec[1] == 1
    for accessor in ("zvec", "nvec", "beta_vec", "se_vec", "eaf_vec", "info_vec"):
        assert np.isnan(getattr(stats, accessor)[1])


def test_swap_state_is_independent_of_missing_hardcalls(tmp_path):
    for suffix in (".bim", ".fam", ".bed"):
        shutil.copyfile((ROOT / "genotype/1").with_suffix(suffix), (tmp_path / "source").with_suffix(suffix))
    bed = tmp_path / "source.bed"
    payload = bytearray(bed.read_bytes())
    payload[3 + 2] = 0x55  # All four subjects are missing at swapped source row 2.
    bed.write_bytes(payload)
    reference = load_reference(ROOT / "reference.bim", shards=["1"])
    genotype = load_genotype(tmp_path / "source", reference)
    assert genotype.is_present[1] and genotype.is_allele_swapped[1]
    np.testing.assert_array_equal(genotype.fetch_genotypes_int8([1, 1]), np.full((4, 2), -1, np.int8))


@pytest.mark.parametrize("value", [None, 0, 1, 2, 0.0, "yes", [], np.array(True), np.array([True])])
@pytest.mark.parametrize("loader,path", [
    (load_sumstats, ROOT / "sumstats.tsv.gz"),
    (load_genotype, ROOT / "genotype/@"),
])
def test_loaders_reject_non_boolean_indel_option(loader, path, value):
    reference = load_reference(ROOT / "reference.bim")
    with pytest.raises(ValueError, match="allow_indel_swaps must be a boolean scalar"):
        loader(path, reference, allow_indel_swaps=value)


@pytest.mark.parametrize("enabled", [False, True])
def test_loaders_accept_numpy_boolean_indel_option(enabled):
    reference = load_reference(ROOT / "reference.bim")
    expected = EXPECTED["options"][str(enabled).lower()]
    genotype = load_genotype(ROOT / "genotype/@", reference, allow_indel_swaps=np.bool_(enabled))
    np.testing.assert_array_equal(genotype.is_allele_swapped, expected["is_allele_swapped"])
    with pytest.warns(RuntimeWarning):
        stats = load_sumstats(ROOT / "sumstats.tsv.gz", reference, allow_indel_swaps=np.bool_(enabled))
    np.testing.assert_allclose(stats.eaf_vec, _as_float(expected["aligned_sumstats"]["eaf"]), equal_nan=True)
