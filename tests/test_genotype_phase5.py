import numpy as np
import pytest

from tests.conftest import (
    FIXTURES_DIR,
    GENOTYPE_CHR1_CALLS,
    copy_genotype_shard_files,
    matlab_data_lines,
    run_octave,
    skipif_no_octave,
    write_plink_bed_calls,
)


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_fixture_mixed_shard_repeated_and_float_wrapper():
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        "gi = g.fetch_genotypes_int8([6 1 5 6 3]); "
        "gf = g.fetch_genotypes([6 1 5 6 3]); "
        "empty_i = g.fetch_genotypes_int8([]); "
        "empty_f = g.fetch_genotypes([]); "
        "fprintf('%d %d\\n', size(gi, 1), size(gi, 2)); "
        "fprintf('%d ', gi); fprintf('\\n'); "
        "fprintf('%s %.0f\\n', class(gf), sum(isnan(gf(:)))); "
        "fprintf('%d %d %s %d %d\\n', size(empty_i, 1), size(empty_i, 2), class(empty_i), size(empty_f, 1), size(empty_f, 2));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "4 5",
        "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2",
        "double 0",
        "4 0 int8 4 0",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_decodes_all_two_bit_states_and_preserves_order(tmp_path):
    for suffix in (".bim", ".fam"):
        (tmp_path / f"chr1{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes()
        )
    write_plink_bed_calls(tmp_path / "chr1.bed", GENOTYPE_CHR1_CALLS)

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}', {{'1'}}); "
        f"g = statgen.load_genotype('{tmp_path / 'chr1'}', ref); "
        "gi = g.fetch_genotypes_int8([2 1 2 5]); "
        "gf = g.fetch_genotypes([2 1 2 5]); "
        "gf(isnan(gf)) = -9; "
        "fprintf('%d ', gi); fprintf('\\n'); "
        "fprintf('%.0f ', gf); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "0 1 2 -1 2 -1 1 0 0 1 2 -1 -1 2 0 1",
        "0 1 2 -9 2 -9 1 0 0 1 2 -9 -9 2 0 1",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_genotypes_ploidy_scaled_maps_haploid_calls():
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        "raw = g.fetch_genotypes(6); "
        "scaled = g.fetch_genotypes(6, 'ploidy_scaled'); "
        "fprintf('%.0f ', raw); fprintf('\\n'); "
        "fprintf('%.0f ', scaled); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "2 2 2 2",
        "1 2 1 2",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_genotypes_ploidy_scaled_rejects_unknown_sex(tmp_path):
    copy_genotype_shard_files(tmp_path, "X")
    (tmp_path / "X.fam").write_text(
        "FAM1\tIND1\t0\t0\t0\t-9\n"
        "FAM1\tIND2\t0\t0\t2\t-9\n"
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM2\tIND4\t0\t0\t2\t-9\n"
    )
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}', {{'X'}}); "
        f"g = statgen.load_genotype('{tmp_path / 'X'}', ref); "
        "raw = g.fetch_genotypes(1); "
        "err = 0; try; g.fetch_genotypes(1, 'ploidy_scaled'); catch; err = 1; end; "
        "fprintf('%d %.0f\\n', err, raw(1));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 2"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_genotypes_ploidy_scaled_allows_unknown_sex_absent_from_chrx(tmp_path):
    copy_genotype_shard_files(tmp_path, "1")
    copy_genotype_shard_files(tmp_path, "X")
    (tmp_path / "1.fam").write_text(
        "FAM1\tIND1\t0\t0\t1\t-9\n"
        "FAM1\tIND2\t0\t0\t0\t-9\n"
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM2\tIND4\t0\t0\t2\t-9\n"
    )
    (tmp_path / "X.fam").write_text(
        "FAM1\tIND1\t0\t0\t1\t-9\n"
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM2\tIND4\t0\t0\t2\t-9\n"
    )
    write_plink_bed_calls(
        tmp_path / "X.bed",
        np.array(
            [
                [2, 2, 2],
                [2, 2, 2],
                [2, 2, 2],
            ],
            dtype=np.int8,
        ),
    )
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{tmp_path / '@'}', ref); "
        "scaled = g.fetch_genotypes(6, 'ploidy_scaled'); "
        "scaled(isnan(scaled)) = -9; "
        "fprintf('%.0f ', scaled); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 -9 1 2"]


@pytest.mark.octave
@skipif_no_octave
@pytest.mark.parametrize("num_sample", [5, 6, 7, 8])
def test_octave_fetch_decodes_partial_final_bed_byte_for_all_sample_count_modulo_4(tmp_path, num_sample):
    (tmp_path / "mod4.bim").write_text(
        "1\tmod4_rs1\t0\t100\tA\tG\n"
        "1\tmod4_rs2\t0\t200\tC\tT\n"
        "1\tmod4_rs3\t0\t300\tG\tA\n",
        encoding="utf-8",
    )
    (tmp_path / "mod4.fam").write_text(
        "".join(f"F{i}\tI{i}\t0\t0\t{1 if i % 2 else 2}\t-9\n" for i in range(1, num_sample + 1)),
        encoding="utf-8",
    )
    base = np.array(
        [
            [0, 1, 2, -1, 0, 1, 2, -1],
            [2, -1, 1, 0, 2, -1, 1, 0],
            [-1, 2, 0, 1, -1, 2, 0, 1],
        ],
        dtype=np.int8,
    )
    calls = base[:, :num_sample]
    write_plink_bed_calls(tmp_path / "mod4.bed", calls)
    expected = " ".join(str(int(x)) for x in calls[[2, 0, 1, 2]].T.reshape(-1, order="F"))

    script = (
        f"ref = statgen.load_reference('{tmp_path / 'mod4.bim'}'); "
        f"g = statgen.load_genotype('{tmp_path / 'mod4'}', ref); "
        "gi = g.fetch_genotypes_int8([3 1 2 3]); "
        "fprintf('%d ', gi); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == [expected]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_chrx_subset_expands_to_panel_sample_axis(tmp_path):
    copy_genotype_shard_files(tmp_path, "1")
    copy_genotype_shard_files(tmp_path, "X")
    (tmp_path / "X.fam").write_text(
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t1\t-9\n"
    )
    write_plink_bed_calls(
        tmp_path / "X.bed",
        np.array(
            [
                [0, 1],
                [2, -1],
                [1, 0],
            ],
            dtype=np.int8,
        ),
    )

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{tmp_path / '@'}', ref); "
        "gi = g.fetch_genotypes_int8(6); "
        "fprintf('%d ', gi); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 -1 0 -1"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_bed_path_overrides_and_validation_failures(tmp_path):
    bad_magic = tmp_path / "bad_magic.bed"
    bad_magic.write_bytes(b"\x00\x00\x01" + b"\x00" * 5)
    sample_major = tmp_path / "sample_major.bed"
    sample_major.write_bytes(b"\x6c\x1b\x00" + b"\x00" * 5)
    trailing = tmp_path / "trailing.bed"
    trailing.write_bytes(b"\x6c\x1b\x01" + b"\x00" * 6)

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"gs = statgen.load_genotype('{G_SHARDED}', ref); "
        f"gn = statgen.load_genotype('{G_NONSHARDED}', ref); "
        f"gs.save_cache('{tmp_path / 's.mat'}'); gn.save_cache('{tmp_path / 'n.mat'}'); "
        f"hs = statgen.load_genotype_cache('{tmp_path / 's.mat'}', [], '{FIXTURES_DIR / 'genotype/sharded/@.bed'}'); "
        f"hn = statgen.load_genotype_cache('{tmp_path / 'n.mat'}', [], '{FIXTURES_DIR / 'genotype/nonsharded/all.bed'}'); "
        "ok = isequal(hs.fetch_genotypes_int8([1 6]), int8(2 * ones(4, 2))); "
        "ok = ok && isequal(hn.fetch_genotypes_int8([1 6]), int8(2 * ones(4, 2))); "
        f"e1 = 0; try; statgen.load_genotype_cache('{tmp_path / 'n.mat'}', [], '{tmp_path / '@.bed'}'); catch; e1 = 1; end; "
        f"e2 = 0; try; statgen.load_genotype_cache('{tmp_path / 's.mat'}', [], '{tmp_path / 'flat.bed'}'); catch ME; e2 = strcmp(ME.identifier, 'statgen:cache') && ~isempty(strfind(ME.message, 'load_genotype_cache: flat bed_path override requires equal source_num_snp and source_num_sample')); end; "
        f"ref1 = statgen.load_reference('{REF_SHARDED}', {{'1'}}); "
        f"g1 = statgen.load_genotype('{G_SHARDED}', ref1); g1.save_cache('{tmp_path / 'one.mat'}'); "
        f"e3 = 0; try; hx = statgen.load_genotype_cache('{tmp_path / 'one.mat'}', [], '{bad_magic}'); hx.fetch_genotypes_int8(1); catch; e3 = 1; end; "
        f"e4 = 0; try; hx = statgen.load_genotype_cache('{tmp_path / 'one.mat'}', [], '{sample_major}'); hx.fetch_genotypes_int8(1); catch; e4 = 1; end; "
        f"e5 = 0; try; statgen.load_genotype_cache('{tmp_path / 'one.mat'}', [], '{trailing}'); catch; e5 = 1; end; "
        "fprintf('%d %d %d %d %d %d\\n', ok, e1, e2, e3, e4, e5);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1 1 1 1 1"]
