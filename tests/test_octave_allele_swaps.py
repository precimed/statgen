"""Stage 2: MATLAB/Octave against the independent portable swap oracle."""
import json
import shutil

import numpy as np
import pytest
from scipy.io import loadmat

from statgen.genotype import load_genotype
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, run_octave, skipif_no_octave

pytestmark = [pytest.mark.octave, skipif_no_octave]
ROOT = FIXTURES_DIR / "allele_swaps"
EXPECTED = json.loads((ROOT / "expected.json").read_text())
CACHE = FIXTURES_DIR / "cache/genotype_matlab_genotype_cache_0_2.mat"


def _run(script):
    result = run_octave(script, timeout=60)
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"
    return result


@pytest.mark.parametrize("enabled", [False, True])
def test_shared_oracle_and_roundtrips(enabled, tmp_path):
    expected = EXPECTED["options"][str(enabled).lower()]
    option = str(enabled).lower()
    output = tmp_path / "aligned.mat"
    cache = tmp_path / "genotype.mat"
    stats_cache = tmp_path / "sumstats.mat"
    present = np.array(expected["source_row0"]) >= 0
    requested = np.r_[np.flatnonzero(present)[::-1], [1, 13, 1]]
    indices = " ".join(str(i + 1) for i in requested)
    result = _run(
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}'); "
        f"g = statgen_load_genotype('{ROOT / 'genotype/@'}', ref, 'allow_indel_swaps', {option}); "
        f"s = statgen_load_sumstats('{ROOT / 'sumstats.tsv.gz'}', ref, 'allow_indel_swaps', {option}); "
        "source_row0 = g.source_row0; swapped = g.is_allele_swapped; present = g.is_present; "
        "assert(islogical(swapped) && isequal(size(swapped), [ref.num_snp, 1])); "
        f"requested = [{indices}]; "
        "raw = g.fetch_genotypes_int8(requested); assert(isa(raw, 'int8')); "
        "floating = g.fetch_genotypes(requested); "
        "scaled = g.fetch_genotypes(requested, 'ploidy_scaled'); "
        "again = g.fetch_genotypes_int8(requested); assert(isequal(raw, again)); "
        f"g.save_cache('{cache}'); h = statgen.load_genotype_cache('{cache}'); "
        "assert(isequal(h.is_allele_swapped, swapped)); "
        "assert(isequal(h.source_row0, source_row0)); "
        "assert(isequal(h.fetch_genotypes_int8(requested), raw)); "
        "assert(isequaln(h.fetch_genotypes(requested, 'ploidy_scaled'), scaled)); "
        f"sx = statgen.load_genotype_cache('{cache}', {{'X'}}); gx = g.select_shards({{'X'}}); "
        "assert(isequal(sx.is_allele_swapped, swapped(14:16))); "
        "assert(isequal(sx.is_allele_swapped, gx.is_allele_swapped)); "
        "assert(isequal(sx.fetch_genotypes_int8([1 2 1]), gx.fetch_genotypes_int8([1 2 1]))); "
        "failed = false; try; g.fetch_genotypes_int8(9); catch ME; failed = ~isempty(strfind(ME.message, 'not present')); end; assert(failed); "
        f"s.save_cache('{stats_cache}'); sc = statgen.load_sumstats_cache('{stats_cache}'); "
        "assert(isequaln(sc.zvec, s.zvec) && isequaln(sc.eaf_vec, s.eaf_vec)); "
        "logp = s.logpvec; z = s.zvec; n = s.nvec; beta = s.beta_vec; se = s.se_vec; eaf = s.eaf_vec; info = s.info_vec; "
        f"save('{output}', 'source_row0', 'swapped', 'present', 'raw', 'floating', 'scaled', 'logp', 'z', 'n', 'beta', 'se', 'eaf', 'info', '-v7');"
    )
    data = loadmat(output, mat_dtype=True)
    np.testing.assert_array_equal(data["source_row0"].ravel(), expected["source_row0"])
    np.testing.assert_array_equal(data["swapped"].ravel(), expected["is_allele_swapped"])
    assert data["swapped"].dtype == bool
    np.testing.assert_array_equal(data["present"].ravel(), present)
    raw = np.array(expected["raw_counts"], dtype=float)[:, requested]
    np.testing.assert_array_equal(data["raw"], np.where(np.isnan(raw), -1, raw).astype(np.int8))
    np.testing.assert_allclose(data["floating"], raw, equal_nan=True)
    np.testing.assert_allclose(data["scaled"], np.array(expected["scaled_counts"], dtype=float)[:, requested], equal_nan=True)
    for field in ("z", "n", "beta", "se", "eaf", "info"):
        np.testing.assert_allclose(data[field].ravel(), np.array(expected["aligned_sumstats"][field], dtype=float), equal_nan=True)
    p = np.array(expected["aligned_sumstats"]["p"], dtype=float)
    np.testing.assert_allclose(data["logp"].ravel(), -np.log10(p), equal_nan=True)
    # Cross-runtime agreement is additional evidence; the oracle above is primary.
    python = load_genotype(ROOT / "genotype/@", load_reference(ROOT / "reference.bim"), allow_indel_swaps=enabled)
    np.testing.assert_array_equal(data["raw"], python.fetch_genotypes_int8(requested))
    combined = result.stdout + result.stderr
    for label, counts in expected["reporting"].items():
        for path in (ROOT / "sumstats.tsv.gz", ROOT / f"genotype/{label}.bim"):
            prefix = f"{path}: shard {label}: "
            assert prefix + f"{counts['swaps']} successful allele swap(s)" in combined
            if counts["strand_ambiguous"]:
                assert prefix + f"{counts['strand_ambiguous']} strand-ambiguous A/T or C/G swap(s)" in combined
            if counts["disabled_indels"]:
                assert prefix + f"{counts['disabled_indels']} indel row(s)" in combined
    assert "0 strand-ambiguous" not in combined
    assert "0 indel row(s)" not in combined
    assert "remain unmatched" not in combined


def test_matcher_ambiguity_precedence_counts_and_hash_precision():
    _run(
        "h = statgen.internal.allele_hash64({'A'; 'G'; 'AG'; 'AC'; 'GT'}); "
        "match = @(bp, a1, a2, enabled) statgen.internal.match_shard_numeric(100, h(1), h(2), bp, a1, a2, '1', 'sumstats', enabled); "
        "for enabled = [false true]; "
        "for swapped = [false true]; "
        "a = h(1 + swapped); b = h(2 - swapped); "
        "failed = false; try; match([100;100], [a;a], [b;b], enabled); catch ME; failed = strcmp(ME.identifier, 'statgen:match'); end; assert(failed); end; "
        "[loc, sw, amb, disabled] = match([100;100;100;999;999], h([2;1;2;1;1]), h([1;2;1;2;2]), enabled); "
        "assert(isequal(loc, 2) && ~sw && amb == 0 && disabled == 0); end; "
        "[loc, sw, amb, disabled] = statgen.internal.match_shard_numeric([100;200;200;300;400], h([1;1;3;1;4]), h([3;3;1;2;5]), "
        "[100;100;200;300;400;999], h([3;3;1;2;5;3]), h([1;1;3;1;4;1]), '1', 'sumstats', false); "
        "assert(isequal(loc, [0;3;0;4;0]) && isequal(sw, logical([0;0;0;1;0])) && disabled == 2); "
        "failed = false; try; statgen.internal.match_shard_numeric(100,h(1),h(3),[100;100],h([3;3]),h([1;1]),'1','sumstats',true); catch ME; failed = strcmp(ME.identifier, 'statgen:match'); end; assert(failed); "
        # Preserve low hash bits above the exact-double bound. Both pairs are
        # multi-base/unclassified synthetic hashes; exact matching must differ.
        "high = bitshift(uint64(1), 60); a = high + uint64(1); b = high + uint64(2); "
        "assert(double(a) == double(b)); "
        "[loc, sw] = statgen.internal.match_shard_numeric([100;100], [a;b], [b;a], [100;100], [b;a], [a;b], '1', 'sumstats', true); "
        "assert(isequal(loc, [2;1]) && ~any(sw)); "
        "[loc, sw] = statgen.internal.match_shard_numeric([], uint64([]), uint64([]), 100, h(1), h(2), '1'); assert(isempty(loc) && isempty(sw));"
    )


@pytest.mark.parametrize("alleles", [("A", "G"), ("A", "T"), ("C", "G")])
def test_reference_rejects_reciprocal_snvs(alleles, tmp_path):
    a, b = alleles
    path = tmp_path / "reference.bim"
    path.write_text("1\tlong\t0\t50\tAC\tGT\n1\tsnv\t0\t75\tA\tC\n"
                    f"1\ta\t0\t100\t{a}\t{b}\n1\tb\t0\t100\t{b}\t{a}\n")
    _run(f"failed = false; try; statgen.load_reference('{path}'); catch ME; failed = ~isempty(strfind(ME.message, 'reference shard 1: reciprocal SNV keys are not allowed at bp=100')); end; assert(failed);")


def test_cache_loaded_reference_matching_does_not_decode_or_revalidate(tmp_path):
    path = tmp_path / "reference.bim"
    path.write_text("1\ta\t0\t100\tA\tG\n1\tb\t0\t100\tC\tT\n")
    source = tmp_path / "source.tsv"
    source.write_text("chr\tbp\ta1\ta2\tp\tz\tn\n1\t100\tA\tG\t0.1\t1\t100\n1\t100\tG\tA\t0.2\t2\t100\n")
    cache = tmp_path / "reference.mat"
    _run(
        f"ref = statgen.load_reference('{path}'); statgen.save_reference_cache(ref, '{cache}'); "
        f"payload = load('{cache}'); "
        "payload.a1_hash64(2) = payload.a2_hash64(1); payload.a2_hash64(2) = payload.a1_hash64(1); "
        "payload.a1_text_by_shard = {'invalid'}; payload.a2_text_by_shard = {'invalid'}; "
        f"save('{cache}', '-struct', 'payload', '-v7'); ref = statgen.load_reference_cache('{cache}'); "
        f"s = statgen.load_sumstats('{source}', ref); assert(isequal(s.zvec, [1;2]));"
    )


def test_quiet_preserves_unrelated_warnings_and_missing_optional_fields(tmp_path):
    source = tmp_path / "source.tsv"
    source.write_text("chr\tbp\ta1\ta2\tp\n1\t200\tG\tA\t0\n1\t100\tA\tG\t0.1\n"
                      "1\t100\tG\tA\t0.2\n1\t100\tG\tA\t0.3\n1\t9999\tA\tC\t0.1\n1\t9999\tA\tC\t0.1\n")
    result = _run(
        "statgen.set_verbosity('quiet'); "
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}'); "
        f"s = statgen.load_sumstats('{source}', ref); "
        "assert(s.is_present(2) && isinf(s.logpvec(2)) && isempty(s.zvec) && isempty(s.eaf_vec)); "
        f"g = statgen.load_genotype('{ROOT / 'genotype/@'}', ref);"
    )
    combined = result.stdout + result.stderr
    assert "successful allele swap" not in combined
    assert "indel row(s)" not in combined
    assert "strand-ambiguous" not in combined
    assert "zvec" in combined and "nvec" in combined  # Warnings still emitted.


def test_swapped_missing_optional_values(tmp_path):
    source = tmp_path / "missing.tsv"
    source.write_text("chr\tbp\ta1\ta2\tp\tz\tn\tbeta\tse\teaf\tinfo\n1\t200\tG\tA\t0.1\tNA\tNA\tNA\tNA\tNA\tNA\n")
    _run(
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}', {{'1'}}); "
        f"s = statgen.load_sumstats('{source}', ref); "
        "assert(s.is_present(2) && s.logpvec(2) == 1); "
        "assert(all(isnan([s.zvec(2), s.nvec(2), s.beta_vec(2), s.se_vec(2), s.eaf_vec(2), s.info_vec(2)])));"
    )


def test_cache_fixture_legacy_migration_and_swapped_overrides(tmp_path):
    source = tmp_path / "source"
    shutil.copytree(ROOT / "genotype", source)
    cache = tmp_path / "genotype.mat"
    migrated = tmp_path / "migrated.mat"
    script = (
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}'); "
        f"g = statgen.load_genotype('{source / '@'}', ref, 'allow_indel_swaps', true); "
        f"fixture = statgen.load_genotype_cache('{CACHE}'); "
        "assert(isequal(fixture.is_allele_swapped, g.is_allele_swapped)); "
        "assert(isequal(fixture.source_row0, g.source_row0)); "
        "requested = [find(g.is_present)' 2 5 14 16 2]; "
        "expected = g.fetch_genotypes_int8(requested); "
        "assert(isequal(fixture.fetch_genotypes_int8(requested), expected)); "
        f"g.save_cache('{cache}'); "
    )
    for file in source.iterdir():
        script += f"delete('{file}'); "
    script += (
        f"h = statgen.load_genotype_cache('{cache}'); "
        f"assert(isequal(h.fetch_genotypes_int8(requested, '{ROOT / 'genotype/@.bed'}'), expected)); "
        f"assert(isequaln(h.fetch_genotypes(requested, '{ROOT / 'genotype/@.bed'}', 'ploidy_scaled'), fixture.fetch_genotypes(requested, 'ploidy_scaled'))); "
        f"old = statgen.load_genotype_cache('{FIXTURES_DIR / 'cache/genotype_matlab_genotype_cache_0_1.mat'}'); "
        "assert(~any(old.is_allele_swapped)); "
        "assert(isequal(old.fetch_genotypes_int8([1 6 1]), int8(2 * ones(4,3)))); "
        f"old.save_cache('{migrated}'); migrated = statgen.load_genotype_cache('{migrated}'); "
        "assert(~any(migrated.is_allele_swapped) && isequal(migrated.source_row0, old.source_row0));"
    )
    _run(script)


@pytest.mark.parametrize("mutation", [
    "payload = rmfield(payload, 'is_allele_swapped')",
    "payload.is_allele_swapped = payload.is_allele_swapped(1:end-1)",
    "payload.is_allele_swapped = reshape(payload.is_allele_swapped, 4, 4)",
    "payload.is_allele_swapped = double(payload.is_allele_swapped)",
    "payload.is_allele_swapped(9) = true",
])
def test_cache_rejects_invalid_swap_flags(mutation, tmp_path):
    cache = tmp_path / "bad.mat"
    _run(
        f"payload = load('{CACHE}'); {mutation}; save('{cache}', '-struct', 'payload', '-v7'); "
        f"failed = false; try; statgen.load_genotype_cache('{cache}'); catch ME; failed = ~isempty(strfind(ME.message, 'is_allele_swapped')); end; assert(failed);"
    )


def test_swap_state_independent_of_missing_hardcalls(tmp_path):
    for suffix in (".bim", ".fam", ".bed"):
        shutil.copyfile((ROOT / "genotype/1").with_suffix(suffix), (tmp_path / "source").with_suffix(suffix))
    bed = tmp_path / "source.bed"
    payload = bytearray(bed.read_bytes())
    payload[3 + 2] = 0x55
    bed.write_bytes(payload)
    _run(
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}', {{'1'}}); "
        f"g = statgen.load_genotype('{tmp_path / 'source'}', ref); "
        "assert(g.is_present(2) && g.is_allele_swapped(2)); "
        "assert(isequal(g.fetch_genotypes_int8([2 2]), int8(-ones(4,2))));"
    )


@pytest.mark.parametrize("args", ["'allow_indel_swaps'", "'unknown', true", "'allow_indel_swaps', 1", "'allow_indel_swaps', [true false]"])
def test_source_option_requires_logical_scalar(args):
    _run(
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}'); "
        f"failed = false; try; statgen.load_genotype('{ROOT / 'genotype/@'}', ref, {args}); catch ME; failed = strcmp(ME.identifier, 'statgen:arg'); end; assert(failed); "
        f"failed = false; try; statgen.load_sumstats('{ROOT / 'sumstats.tsv.gz'}', ref, {args}); catch ME; failed = strcmp(ME.identifier, 'statgen:arg'); end; assert(failed);"
    )


def test_nonsharded_physical_rows_and_ploidy_survive_reference_projection(tmp_path):
    # An unselected supported chromosome precedes X in one physical BED file.
    # The mapping must use physical rows rather than per-reference-shard rows.
    from tests.conftest import write_plink_bed_calls

    prefix = tmp_path / "combined"
    prefix.with_suffix(".bim").write_text(
        "1\tignored\t0\t1\tA\tC\n"
        "X\tsnv\t0\t100\tG\tA\n"
        "X\texact\t0\t200\tC\tT\n"
        "X\tindel\t0\t300\tAT\tA\n"
    )
    shutil.copyfile(ROOT / "genotype/1.fam", prefix.with_suffix(".fam"))
    prefix.with_suffix(".ploidy").write_text("2\t2\n1\t2\n2\t2\n0\t1\n")
    write_plink_bed_calls(prefix.with_suffix(".bed"), [[2, 2, 2, 2], [0, 1, 2, -1], [2, 1, 0, -1], [0, 1, 2, -1]])
    cache = tmp_path / "combined.mat"
    _run(
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}', {{'X'}}); "
        f"g = statgen.load_genotype('{prefix}', ref, 'allow_indel_swaps', true); "
        "assert(isequal(g.source_row0, [1;2;3]) && isequal(g.is_allele_swapped, logical([1;0;1]))); "
        "assert(isequal(g.ploidy_male, [1;2;0]) && isequal(g.ploidy_female, [2;2;1])); "
        "assert(isequal(g.fetch_genotypes_int8([1 3 1]), int8([2 2 2;1 1 1;0 0 0;-1 -1 -1]))); "
        "assert(isequaln(g.fetch_genotypes([1 3], 'ploidy_scaled'), [1 0;1 0.5;0 0;NaN NaN])); "
        f"g.save_cache('{cache}'); h = statgen.load_genotype_cache('{cache}'); "
        f"assert(isequal(h.fetch_genotypes_int8([1 3 1], '{prefix.with_suffix('.bed')}'), g.fetch_genotypes_int8([1 3 1])));"
    )


@pytest.mark.parametrize("verbosity", ["info", "quiet"])
@pytest.mark.parametrize("source_row, should_warn", [
    ("X\t100\tA\tG", True),       # No source rows on the reference shard.
    ("1\t9999\tA\tC", True),      # Source rows exist, but none match.
    ("1\t500\tAG\tA", True),      # Disabled indel fallback leaves no match.
    ("1\t200\tG\tA", False),      # An automatic SNV swap counts as a match.
])
def test_no_matching_sumstats_rows_warning(source_row, should_warn, verbosity, tmp_path):
    source = tmp_path / "source.tsv"
    source.write_text("chr\tbp\ta1\ta2\tp\tz\tn\n" + source_row + "\t0.1\t1\t100\n")
    result = _run(
        f"statgen.set_verbosity('{verbosity}'); "
        f"ref = statgen.load_reference('{ROOT / 'reference.bim'}', {{'1'}}); "
        f"s = statgen.load_sumstats('{source}', ref); "
        f"assert(any(s.is_present) == {str(not should_warn).lower()}); "
        "if ~any(s.is_present); assert(all(isnan(s.logpvec))); end;"
    )
    warning = f"{source}: no matching source rows for reference shard '1'"
    assert (warning in result.stderr + result.stdout) == should_warn
