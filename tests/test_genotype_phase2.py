import math

import numpy as np
import pytest

from statgen.genotype import load_genotype, load_genotype_cache
from statgen.reference import load_reference
from tests.conftest import (
    FIXTURES_DIR,
    GENOTYPE_CHR1_CALLS,
    GENOTYPE_CHR1_BIM,
    copy_genotype_shard_files,
    copy_sharded_genotype,
    write_plink_bim,
    write_plink_bed_calls,
    write_plink_fam,
)


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


def test_fetch_committed_fixture_one_multi_cross_shard_and_repeated():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)

    one = panel.fetch_genotypes_int8([0])
    assert one.shape == (4, 1)
    assert one.tolist() == [[2], [2], [2], [2]]

    empty = panel.fetch_genotypes_int8([])
    assert empty.shape == (4, 0)
    assert empty.dtype == np.int8

    empty_float = panel.fetch_genotypes([])
    assert empty_float.shape == (4, 0)
    assert empty_float.dtype == np.float64

    requested = [5, 0, 4, 5, 2]
    geno = panel.fetch_genotypes_int8(requested)
    assert geno.shape == (4, len(requested))
    assert np.array_equal(geno, np.full((4, len(requested)), 2, dtype=np.int8))


def test_fetch_decodes_all_two_bit_states_and_preserves_order(tmp_path):
    for suffix in (".bim", ".fam"):
        (tmp_path / f"chr1{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes()
        )
    write_plink_bed_calls(tmp_path / "chr1.bed", GENOTYPE_CHR1_CALLS)

    ref = load_reference(REF_SHARDED, shards=["1"])
    panel = load_genotype(tmp_path / "chr1", ref)

    requested = [1, 0, 1, 4]
    expected = GENOTYPE_CHR1_CALLS[requested].T
    geno_int8 = panel.fetch_genotypes_int8(requested)
    assert np.array_equal(geno_int8, expected)

    geno = panel.fetch_genotypes(requested)
    assert geno.dtype == np.float64
    assert math.isnan(geno[3, 0])
    expected_float = expected.astype(np.float64)
    expected_float[expected == -1] = np.nan
    assert np.array_equal(geno, expected_float, equal_nan=True)


def test_load_genotype_orients_swapped_alleles(tmp_path, capsys):
    copy_genotype_shard_files(tmp_path, "1")
    rows = list(GENOTYPE_CHR1_BIM)
    rows[1] = ("1", "rs1002", 0, 200, "T", "C")
    write_plink_bim(tmp_path / "1.bim", rows)

    ref = load_reference(REF_SHARDED, shards=["1"])
    panel = load_genotype(tmp_path / "1", ref)
    assert "shard 1: 1 successful allele swap" in capsys.readouterr().err

    np.testing.assert_array_equal(
        panel.is_present,
        np.array([True, True, True, True, True]),
    )
    assert panel.source_row0[1] == 1
    assert panel.is_allele_swapped[1]


def test_fetch_genotypes_ploidy_scaled_maps_haploid_calls():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)

    raw = panel.fetch_genotypes([5])
    scaled = panel.fetch_genotypes([5], haploid_mode="ploidy_scaled")

    assert raw[:, 0].tolist() == [2.0, 2.0, 2.0, 2.0]
    assert scaled[:, 0].tolist() == [1.0, 2.0, 1.0, 2.0]


def test_fetch_genotypes_ploidy_scaled_rejects_unknown_sex_when_ploidy_differs(tmp_path):
    copy_genotype_shard_files(tmp_path, "X")
    write_plink_fam(
        tmp_path / "X.fam",
        [
            ("FAM1", "IND1", 0, 0, 0, -9),
            ("FAM1", "IND2", 0, 0, 2, -9),
            ("FAM2", "IND3", 0, 0, 1, -9),
            ("FAM2", "IND4", 0, 0, 2, -9),
        ],
    )
    ref = load_reference(REF_SHARDED, shards=["X"])
    panel = load_genotype(tmp_path / "X", ref)

    assert panel.fetch_genotypes([0]).shape == (4, 1)
    with pytest.raises(ValueError, match="requires known FAM sex"):
        panel.fetch_genotypes([0], haploid_mode="ploidy_scaled")


def test_fetch_genotypes_ploidy_scaled_allows_unknown_sex_absent_from_chrx(tmp_path):
    copy_genotype_shard_files(tmp_path, "1")
    copy_genotype_shard_files(tmp_path, "X")
    write_plink_fam(
        tmp_path / "1.fam",
        [
            ("FAM1", "IND1", 0, 0, 1, -9),
            ("FAM1", "IND2", 0, 0, 0, -9),
            ("FAM2", "IND3", 0, 0, 1, -9),
            ("FAM2", "IND4", 0, 0, 2, -9),
        ],
    )
    write_plink_fam(
        tmp_path / "X.fam",
        [
            ("FAM1", "IND1", 0, 0, 1, -9),
            ("FAM2", "IND3", 0, 0, 1, -9),
            ("FAM2", "IND4", 0, 0, 2, -9),
        ],
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
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(str(tmp_path / "@"), ref)

    scaled = panel.fetch_genotypes([5], haploid_mode="ploidy_scaled")
    assert np.array_equal(
        scaled[:, 0],
        np.array([1.0, np.nan, 1.0, 2.0]),
        equal_nan=True,
    )


def test_fetch_genotypes_rejects_unknown_haploid_mode():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)

    with pytest.raises(ValueError, match="haploid_mode"):
        panel.fetch_genotypes([0], haploid_mode="normalized")


@pytest.mark.parametrize("num_sample", [5, 6, 7, 8])
def test_fetch_decodes_partial_final_bed_byte_for_all_sample_count_modulo_4(tmp_path, num_sample):
    rows = [
        ("1", "mod4_rs1", 0, 100, "A", "G"),
        ("1", "mod4_rs2", 0, 200, "C", "T"),
        ("1", "mod4_rs3", 0, 300, "G", "A"),
    ]
    (tmp_path / "mod4.bim").write_text(
        "".join(f"{chrom}\t{snp}\t{cm}\t{bp}\t{a1}\t{a2}\n" for chrom, snp, cm, bp, a1, a2 in rows),
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

    ref = load_reference(tmp_path / "mod4.bim")
    panel = load_genotype(tmp_path / "mod4", ref)
    requested = [2, 0, 1, 2]

    got = panel.fetch_genotypes_int8(requested)
    assert got.shape == (num_sample, len(requested))
    np.testing.assert_array_equal(got, calls[requested].T)


def test_fetch_chrx_subset_expands_to_panel_sample_axis(tmp_path):
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

    ref = load_reference(REF_SHARDED)
    panel = load_genotype(str(tmp_path / "@"), ref)

    geno = panel.fetch_genotypes_int8([5])
    assert geno[:, 0].tolist() == [1, -1, 0, -1]


def test_absent_snp_fails_before_reading_bed_payload(tmp_path):
    extra_bim = (
        "1\trs1001\t0\t100\tA\tG\n"
        "1\trs1002\t0\t200\tC\tT\n"
        "1\trs1003\t0\t300\tA\tC\n"
        "1\trs1004\t0\t400\tG\tA\n"
        "1\trs1005\t0\t500\tT\tC\n"
        "1\trs1006\t0\t600\tA\tC\n"
    )
    (tmp_path / "extra.bim").write_text(extra_bim)
    ref = load_reference(tmp_path / "extra.bim")
    source = tmp_path / "source"
    source.mkdir()
    copy_sharded_genotype(source)
    panel = load_genotype(str(source / "@"), ref)
    (source / "1.bed").unlink()

    with pytest.raises(ValueError, match="not present in the genotype source"):
        panel.fetch_genotypes_int8([5])


def test_absent_snp_error_reports_panel_global_index_for_later_shard(tmp_path):
    (tmp_path / "1.bim").write_bytes((FIXTURES_DIR / "reference/sharded/1.bim").read_bytes())
    (tmp_path / "X.bim").write_text(
        (FIXTURES_DIR / "reference/sharded/X.bim").read_text()
        + "X\trsX004\t0\t400\tG\tA\n"
    )
    ref = load_reference(str(tmp_path / "@.bim"))
    source = tmp_path / "source"
    source.mkdir()
    copy_sharded_genotype(source)
    panel = load_genotype(str(source / "@"), ref)
    (source / "X.bed").unlink()

    with pytest.raises(ValueError, match=r"requested SNP 8 in shard 'X' is not present"):
        panel.fetch_genotypes_int8([8])


def test_flat_bed_override_works_for_nonsharded_panel(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_NONSHARDED, ref)

    cache = tmp_path / "g.npz"
    panel.save_cache(cache)
    geno = load_genotype_cache(cache, bed_path=FIXTURES_DIR / "genotype/nonsharded/all.bed").fetch_genotypes_int8([0, 5])
    assert np.array_equal(geno, np.full((4, 2), 2, dtype=np.int8))


def test_nonsharded_fetches_chrx_snps_from_shared_bed_without_override():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_NONSHARDED, ref)

    assert panel.source_layout == "non_sharded"
    assert panel.is_subject_present("X").tolist() == [True, True, True, True]
    assert panel.shards[1].source_subject_row0.tolist() == [0, 1, 2, 3]

    geno = panel.fetch_genotypes_int8([5, 7, 0])
    assert np.array_equal(geno, np.full((4, 3), 2, dtype=np.int8))


def test_at_bed_override_works_for_sharded_panel_and_single_shard_subset(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)
    pattern = str(FIXTURES_DIR / "genotype/sharded/@.bed")

    cache = tmp_path / "g.npz"
    panel.save_cache(cache)
    geno = load_genotype_cache(cache, bed_path=pattern).fetch_genotypes_int8([0, 5])
    assert np.array_equal(geno, np.full((4, 2), 2, dtype=np.int8))

    chrx = load_genotype_cache(cache, shards=["X"], bed_path=pattern)
    x_geno = chrx.fetch_genotypes_int8([0])
    assert np.array_equal(x_geno, np.full((4, 1), 2, dtype=np.int8))


def test_at_override_against_nonsharded_panel_fails_before_reading(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_NONSHARDED, ref)

    cache = tmp_path / "g.npz"
    panel.save_cache(cache)
    with pytest.raises(ValueError, match="@ override incompatible with non-sharded"):
        load_genotype_cache(cache, bed_path=str(tmp_path / "@.bed"))


def test_flat_override_against_unequal_sharded_sources_fails_before_reading(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)

    cache = tmp_path / "g.npz"
    panel.save_cache(cache)
    with pytest.raises(ValueError, match="equal source_num_snp and source_num_sample across selected shards"):
        load_genotype_cache(cache, bed_path=tmp_path / "flat.bed")


def test_bed_override_size_and_header_validation(tmp_path):
    ref = load_reference(REF_SHARDED, shards=["1"])
    panel = load_genotype(G_SHARDED, ref)
    cache = tmp_path / "g.npz"
    panel.save_cache(cache)
    expected_size = 3 + 5

    bad_magic = tmp_path / "bad_magic.bed"
    bad_magic.write_bytes(b"\x00\x00\x01" + b"\x00" * (expected_size - 3))
    with pytest.raises(ValueError, match="invalid PLINK BED magic"):
        load_genotype_cache(cache, bed_path=bad_magic).fetch_genotypes_int8([0])

    sample_major = tmp_path / "sample_major.bed"
    sample_major.write_bytes(b"\x6c\x1b\x00" + b"\x00" * (expected_size - 3))
    with pytest.raises(ValueError, match="SNP-major mode"):
        load_genotype_cache(cache, bed_path=sample_major).fetch_genotypes_int8([0])

    trailing = tmp_path / "trailing.bed"
    trailing.write_bytes(b"\x6c\x1b\x01" + b"\x00" * (expected_size - 2))
    with pytest.raises(ValueError, match="incorrectly sized"):
        load_genotype_cache(cache, bed_path=trailing)
