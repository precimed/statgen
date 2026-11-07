import hashlib
import numpy as np
import pytest
from pathlib import Path

from statgen.reference import (
    ReferenceShard,
    ReferencePanel,
    load_reference,
    load_reference_cache,
    save_reference_cache,
)
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, run_octave, skipif_no_octave, write_plink_bim

SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
NONSHARDED = FIXTURES_DIR / "reference/nonsharded/all.bim"

# Fixture data (matches generate.py)
CHR1_ROWS = [
    ("1", "rs1001", 0, 100, "A", "G"),
    ("1", "rs1002", 0, 200, "C", "T"),
    ("1", "rs1003", 0, 300, "A", "C"),
    ("1", "rs1004", 0, 400, "G", "A"),
    ("1", "rs1005", 0, 500, "T", "C"),
]
CHRX_ROWS = [
    ("X", "rsX001", 0, 100, "A", "G"),
    ("X", "rsX002", 0, 200, "C", "T"),
    ("X", "rsX003", 0, 300, "A", "C"),
]


def _checksum(rows):
    text = "".join(f"{r[0]}:{r[3]}:{r[4]}:{r[5]}\n" for r in rows)
    return hashlib.md5(text.encode()).hexdigest()


CHR1_CHECKSUM = _checksum(CHR1_ROWS)
CHRX_CHECKSUM = _checksum(CHRX_ROWS)


# ---------------------------------------------------------------------------
# Python: loading
# ---------------------------------------------------------------------------

def test_sharded_shard_labels():
    panel = load_reference(SHARDED)
    assert [s.label for s in panel.shards] == ["1", "X"]


def test_sharded_num_snp():
    panel = load_reference(SHARDED)
    assert panel.shards[0].num_snp == 5
    assert panel.shards[1].num_snp == 3
    assert panel.num_snp == 8


def test_sharded_checksums():
    panel = load_reference(SHARDED)
    assert panel.shards[0].checksum == CHR1_CHECKSUM
    assert panel.shards[1].checksum == CHRX_CHECKSUM


def test_reference_allele_hashes_are_uint64_and_panel_wide():
    panel = load_reference(SHARDED)
    assert panel.a1_hash64.dtype == np.uint64
    assert panel.a2_hash64.dtype == np.uint64
    assert panel.a1_hash64.shape == (panel.num_snp,)
    assert panel.a2_hash64.shape == (panel.num_snp,)
    assert np.array_equal(panel.a1_hash64[:5], panel.shards[0].a1_hash64)
    assert np.array_equal(panel.a2_hash64[5:], panel.shards[1].a2_hash64)


def test_reference_variant_type_accessors_are_hash_backed():
    shard = ReferenceShard(
        "1",
        ["1"] * 6,
        ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
        [100, 200, 300, 400, 500, 600],
        ["A", "T", "C", "A", "AC", "A"],
        ["T", "A", "G", "C", "G", "AT"],
    )
    panel = ReferencePanel([shard])

    expected_single_base = np.array([True, True, True, True, False, False])
    expected_ambiguous = np.array([True, True, True, False, False, False])
    assert np.array_equal(shard.is_single_nucleotide_variant, expected_single_base)
    assert np.array_equal(panel.is_single_nucleotide_variant, expected_single_base)
    assert np.array_equal(shard.is_strand_ambiguous, expected_ambiguous)
    assert np.array_equal(panel.is_strand_ambiguous, expected_ambiguous)


def test_validate_checksums_passes_for_loaded_reference():
    panel = load_reference(SHARDED)
    assert panel.validate_checksums() is True


def test_nonsharded_split_by_chr():
    panel = load_reference(NONSHARDED)
    assert [s.label for s in panel.shards] == ["1", "X"]
    assert panel.shards[0].num_snp == 5
    assert panel.shards[1].num_snp == 3


def test_nonsharded_split_checksums_match_sharded():
    sharded = load_reference(SHARDED)
    split = load_reference(NONSHARDED)
    for s1, s2 in zip(sharded.shards, split.shards):
        assert s1.checksum == s2.checksum


def test_load_reference_shards_subset():
    panel = load_reference(SHARDED, shards=["X"])
    assert [s.label for s in panel.shards] == ["X"]
    assert panel.num_snp == 3


# ---------------------------------------------------------------------------
# Python: accessors
# ---------------------------------------------------------------------------

def test_shard_offsets():
    panel = load_reference(SHARDED)
    offsets = panel.shard_offsets
    assert offsets[0] == {"shard_label": "1", "start0": 0, "stop0": 5}
    assert offsets[1] == {"shard_label": "X", "start0": 5, "stop0": 8}


def test_panel_chr_vector():
    panel = load_reference(SHARDED)
    assert list(panel.chr) == ["1"] * 5 + ["X"] * 3


def test_panel_bp_vector():
    panel = load_reference(SHARDED)
    assert list(panel.bp[:5]) == [100, 200, 300, 400, 500]
    assert list(panel.bp[5:]) == [100, 200, 300]


def test_panel_snp_vector():
    panel = load_reference(SHARDED)
    assert list(panel.snp[:3]) == ["rs1001", "rs1002", "rs1003"]
    assert list(panel.snp[5:]) == ["rsX001", "rsX002", "rsX003"]


def test_panel_a1_a2():
    panel = load_reference(SHARDED)
    assert panel.a1[0] == "A"
    assert panel.a2[0] == "G"


def test_chrX_shard_accessible():
    panel = load_reference(SHARDED)
    x_shard = panel.shards[1]
    assert x_shard.label == "X"
    assert list(x_shard.chr) == ["X", "X", "X"]


def test_shard_does_not_store_cm():
    panel = load_reference(SHARDED)
    assert hasattr(panel.shards[0], "bp")
    assert not hasattr(panel.shards[0], "cm")


def test_select_shards():
    panel = load_reference(SHARDED)
    subset = panel.select_shards(["X"])
    assert [s.label for s in subset.shards] == ["X"]
    assert list(subset.bp) == [100, 200, 300]


# ---------------------------------------------------------------------------
# Python: is_object_compatible
# ---------------------------------------------------------------------------

class _FakeShard:
    def __init__(self, label, num_snp, checksum=None):
        self.label = label
        self.num_snp = num_snp
        if checksum is not None:
            self.checksum = checksum


class _FakePanel:
    def __init__(self, shards):
        self.shards = shards


def test_compat_same_panel():
    p = load_reference(SHARDED)
    q = load_reference(SHARDED)
    assert p.is_object_compatible(q) is True


def test_compat_no_shards_attr():
    p = load_reference(SHARDED)
    assert p.is_object_compatible(object()) is False


def test_compat_shard_count_mismatch():
    p = load_reference(SHARDED)
    fake = _FakePanel([_FakeShard("1", 5)])  # only one shard instead of two
    assert p.is_object_compatible(fake) is False


def test_compat_row_count_mismatch():
    p = load_reference(SHARDED)
    # chr1 has 5 SNPs, chrX has 3 — swap the counts
    fake = _FakePanel([_FakeShard("1", 3), _FakeShard("X", 5)])
    assert p.is_object_compatible(fake) is False


def test_compat_checksum_mismatch():
    p = load_reference(SHARDED)
    fake = _FakePanel([
        _FakeShard("1", 5, checksum="deadbeef" * 4),
        _FakeShard("X", 3, checksum=CHRX_CHECKSUM),
    ])
    assert p.is_object_compatible(fake) is False


def test_compat_no_checksum_passes():
    p = load_reference(SHARDED)
    # Shards with correct num_snp but no checksum attribute → compatible
    fake = _FakePanel([_FakeShard("1", 5), _FakeShard("X", 3)])
    assert p.is_object_compatible(fake) is True


def test_compat_shard_label_mismatch():
    p = load_reference(SHARDED)
    fake = _FakePanel([_FakeShard("X", 5), _FakeShard("1", 3)])
    assert p.is_object_compatible(fake) is False


def test_compat_does_not_raise(caplog):
    import logging
    p = load_reference(SHARDED)
    fake = _FakePanel([_FakeShard("1", 99), _FakeShard("X", 3)])
    with caplog.at_level(logging.WARNING, logger="statgen.reference"):
        result = p.is_object_compatible(fake)
    assert result is False  # no exception raised


# ---------------------------------------------------------------------------
# Python: cache roundtrip
# ---------------------------------------------------------------------------

def test_cache_roundtrip(tmp_path):
    panel = load_reference(SHARDED)
    cache = tmp_path / "ref.npz"
    save_reference_cache(panel, cache)
    loaded = load_reference_cache(cache)

    assert loaded.num_snp == panel.num_snp
    assert len(loaded.shards) == len(panel.shards)
    for s1, s2 in zip(panel.shards, loaded.shards):
        assert s1.label == s2.label
        assert s1.checksum == s2.checksum
        assert list(s1.chr) == list(s2.chr)
        assert list(s1.snp) == list(s2.snp)
        assert list(s1.bp) == list(s2.bp)
        assert list(s1.a1) == list(s2.a1)
        assert list(s1.a2) == list(s2.a2)
        assert not hasattr(s2, "cm")

def test_cache_shards_subset(tmp_path):
    panel = load_reference(SHARDED)
    cache = tmp_path / "ref.npz"
    panel.save_cache(cache)
    loaded = load_reference_cache(cache, shards=["X"])
    assert [s.label for s in loaded.shards] == ["X"]
    assert loaded.num_snp == 3


def test_cache_is_binary(tmp_path):
    panel = load_reference(SHARDED)
    cache = tmp_path / "ref.npz"
    save_reference_cache(panel, cache)
    with open(cache, "rb") as f:
        magic = f.read(2)
    assert magic == b"PK"  # ZIP magic bytes — confirms binary npz, not JSON
    with np.load(cache, allow_pickle=False) as data:
        assert not any(k.endswith("_cm") for k in data.files)


def test_cache_wrong_schema(tmp_path):
    # Write a valid npz with a bad schema in _meta
    import json
    bad_meta = json.dumps({"schema": "bad/99", "shard_labels": [], "shard_checksums": []}).encode()
    bad_npz = tmp_path / "bad.npz"
    np.savez_compressed(bad_npz, _meta=np.frombuffer(bad_meta, dtype=np.uint8))
    with pytest.raises(ValueError, match="Unsupported"):
        load_reference_cache(bad_npz)


def test_cache_missing_shard_payload_rejected(tmp_path):
    import json
    bad_meta = json.dumps({
        "schema": "reference_cache/0.1",
        "shard_labels": ["1"],
        "shard_checksums": [CHR1_CHECKSUM],
    }).encode()
    bad_npz = tmp_path / "bad.npz"
    np.savez_compressed(bad_npz, _meta=np.frombuffer(bad_meta, dtype=np.uint8))
    with pytest.raises(ValueError, match="a1_hash64/a2_hash64"):
        load_reference_cache(bad_npz)


def test_cache_load_trusts_checksum_until_explicit_validation(tmp_path):
    panel = load_reference(SHARDED)
    cache = tmp_path / "ref.npz"
    save_reference_cache(panel, cache)

    with np.load(cache, allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    arrays["s0_a1_text"] = arrays["s0_a1_text"].copy()
    arrays["s0_a1_text"][0] = ord("T")
    bad_cache = tmp_path / "bad_payload.npz"
    np.savez_compressed(bad_cache, **arrays)

    loaded = load_reference_cache(bad_cache)
    assert loaded.shards[0].checksum == CHR1_CHECKSUM
    with pytest.raises(ValueError, match="Reference checksum mismatch for shard 1"):
        loaded.validate_checksums()


# ---------------------------------------------------------------------------
# Python: BIM validation
# ---------------------------------------------------------------------------

def test_bim_wrong_column_count(tmp_path):
    bad = tmp_path / "bad.bim"
    bad.write_text("1\trs1\t0\t100\tA\n")  # 5 columns
    with pytest.raises(ValueError, match="6 whitespace-delimited"):
        load_reference(bad)


def test_bim_empty_chr(tmp_path):
    bad = tmp_path / "bad.bim"
    bad.write_text("\trs1\t0\t100\tA\tG\n")
    with pytest.raises(ValueError, match="chr"):
        load_reference(bad)


def test_bim_bad_bp(tmp_path):
    bad = tmp_path / "bad.bim"
    bad.write_text("1\trs1\t0\tnot_an_int\tA\tG\n")
    with pytest.raises(ValueError, match="bp"):
        load_reference(bad)


def test_bim_unsorted_bp_fails(tmp_path):
    bad = tmp_path / "unsorted_bp.bim"
    bad.write_text(
        "1\trs1\t0\t200\tA\tG\n"
        "1\trs2\t0\t100\tC\tT\n"
    )
    with pytest.raises(ValueError, match=r"\(chr_rank, bp\)"):
        load_reference(bad)


def test_bim_unsorted_chr_fails(tmp_path):
    bad = tmp_path / "unsorted_chr.bim"
    bad.write_text(
        "X\trsX\t0\t100\tA\tG\n"
        "1\trs1\t0\t200\tC\tT\n"
    )
    with pytest.raises(ValueError, match=r"\(chr_rank, bp\)"):
        load_reference(bad)


def test_bim_noncanonical_chr_fails(tmp_path):
    bad = tmp_path / "noncanonical_chr.bim"
    bad.write_text("chr1\trs1\t0\t100\tA\tG\n")
    with pytest.raises(ValueError, match="chr-style"):
        load_reference(bad)


def test_bim_ignores_y_mt_rows(tmp_path):
    bad = tmp_path / "with_y_mt.bim"
    bad.write_text(
        "1\trs1\t0\t100\tA\tG\n"
        "Y\trsY\t0\t150\tA\tC\n"
        "MT\trsMT\t0\t180\tG\tA\n"
        "X\trsX\t0\t200\tC\tT\n"
    )
    panel = load_reference(bad)
    assert [s.label for s in panel.shards] == ["1", "X"]
    assert panel.num_snp == 2


def test_bim_bad_allele_syntax_fails(tmp_path):
    bad = tmp_path / "bad_allele.bim"
    bad.write_text("1\trs1\t0\t100\tAa\tG\n")
    with pytest.raises(ValueError, match="uppercase DNA"):
        load_reference(bad)


def test_bim_equal_a1_a2_fails(tmp_path):
    bad = tmp_path / "same_alleles.bim"
    bad.write_text("1\trs1\t0\t100\tA\tA\n")
    with pytest.raises(ValueError, match="a1 and a2 must differ"):
        load_reference(bad)


def test_bim_multibase_alleles_are_valid(tmp_path):
    path = tmp_path / "multibase.bim"
    path.write_text(
        "1\trs1\t0\t100\tAC\tG\n"
        "1\trs2\t0\t200\tA\tGTT\n"
    )
    panel = load_reference(path)
    assert list(panel.a1) == ["AC", "A"]
    assert list(panel.a2) == ["G", "GTT"]


def test_bim_allele_order_not_validated_when_chr_bp_sorted(tmp_path):
    bad = tmp_path / "unsorted_a1.bim"
    bad.write_text(
        "1\trs1\t0\t100\tC\tA\n"
        "1\trs2\t0\t100\tA\tG\n"
    )
    panel = load_reference(bad)
    assert panel.num_snp == 2


def test_bim_duplicate_tuple_fails(tmp_path):
    bad = tmp_path / "dupe_tuple.bim"
    bad.write_text(
        "1\trs1\t0\t100\tA\tG\n"
        "1\trs2\t0\t100\tA\tG\n"
    )
    with pytest.raises(ValueError, match="duplicate"):
        load_reference(bad)


def test_bim_uses_dataframe_reader(monkeypatch):
    import statgen.reference as ref_mod

    called = {"read_csv": False}
    original = ref_mod.pd.read_csv

    def _wrapped_read_csv(*args, **kwargs):
        called["read_csv"] = True
        return original(*args, **kwargs)

    monkeypatch.setattr(ref_mod.pd, "read_csv", _wrapped_read_csv)
    panel = load_reference(NONSHARDED)

    assert called["read_csv"] is True
    assert panel.num_snp == 8


def test_cache_load_skips_source_bim_parsing(monkeypatch, tmp_path):
    import statgen.reference as ref_mod

    panel = load_reference(SHARDED)
    cache = tmp_path / "ref.npz"
    save_reference_cache(panel, cache)

    def _fail(*_args, **_kwargs):
        raise AssertionError("source BIM parsing should not run on cache load")

    monkeypatch.setattr(ref_mod, "_parse_bim", _fail)
    loaded = load_reference_cache(cache)
    assert loaded.num_snp == panel.num_snp


def test_bim_sharded_no_match(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_reference(str(tmp_path / "@.bim"))


def test_invalid_shard_subset_errors():
    panel = load_reference(SHARDED)
    with pytest.raises(ValueError, match="canonical subsequence order"):
        panel.select_shards(["X", "1"])
    with pytest.raises(ValueError, match="duplicate"):
        panel.select_shards(["1", "1"])
    with pytest.raises(ValueError, match="not present"):
        panel.select_shards(["2"])


# ---------------------------------------------------------------------------
# Octave: cross-language parity
# ---------------------------------------------------------------------------

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
def test_octave_sharded_num_snp():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "fprintf('%d\\n', ref.num_snp);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "8"


@pytest.mark.octave
@skipif_no_octave
def test_octave_sharded_shard_labels():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "for i = 1:numel(ref.shards); fprintf('%s\\n', ref.shards{i}.label); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines == ["1", "X"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_checksums_match_python():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "for i = 1:numel(ref.shards); fprintf('%s\\n', ref.shards{i}.checksum); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == CHR1_CHECKSUM
    assert lines[1] == CHRX_CHECKSUM


@pytest.mark.octave
@skipif_no_octave
def test_octave_reference_checksum_ignores_matlab_num2str_padding(tmp_path):
    rows = [
        ("20", "rs1", 0, 82090, "A", "C"),
        ("20", "rs2", 0, 102060, "C", "T"),
        ("20", "rs3", 0, 10000000, "G", "A"),
    ]
    bim_path = tmp_path / "reference_chr20.bim"
    write_plink_bim(bim_path, rows)
    expected = _checksum(rows)
    script = (
        f"ref = statgen.load_reference('{bim_path}'); "
        "fprintf('%s\\n', ref.shards{1}.checksum);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines == [expected]


@pytest.mark.octave
@skipif_no_octave
def test_octave_reference_variant_type_accessors():
    script = _octave_script(
        "shard = statgen.ReferenceShard('1', "
        "{'1'; '1'; '1'; '1'; '1'; '1'}, "
        "{'rs1'; 'rs2'; 'rs3'; 'rs4'; 'rs5'; 'rs6'}, "
        "[100; 200; 300; 400; 500; 600], "
        "{'A'; 'T'; 'C'; 'A'; 'AC'; 'A'}, "
        "{'T'; 'A'; 'G'; 'C'; 'G'; 'AT'}); "
        "ref = statgen.ReferencePanel({shard}); "
        "fprintf('%d ', ref.is_single_nucleotide_variant); fprintf('\\n'); "
        "fprintf('%d ', ref.is_strand_ambiguous); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0].strip() == "1 1 1 1 0 0"
    assert lines[1].strip() == "1 1 1 0 0 0"


@pytest.mark.octave
@skipif_no_octave
def test_octave_reference_loader_equal_a1_a2_fails(tmp_path):
    bad = tmp_path / "same_alleles.bim"
    bad.write_text("1\trs1\t0\t100\tA\tA\n")
    script = _octave_script(
        "ok = 0; "
        "try; "
        f"  statgen.load_reference('{bad}'); "
        "catch ME; "
        "  ok = ~isempty(strfind(ME.message, 'a1 and a2 must differ')); "
        "end; "
        "fprintf('%d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_nonsharded_split():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/nonsharded/all.bim']); "
        "fprintf('%d\\n', numel(ref.shards)); "
        "for i = 1:numel(ref.shards); fprintf('%s %d\\n', ref.shards{i}.label, ref.shards{i}.num_snp); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "2"
    assert lines[1] == "1 5"
    assert lines[2] == "X 3"


@pytest.mark.octave
@skipif_no_octave
def test_octave_shard_offsets():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "off = ref.shard_offsets; "
        "for i = 1:numel(off); "
        "  fprintf('%s %d %d\\n', off(i).shard_label, off(i).start0, off(i).stop0); "
        "end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1 0 5"
    assert lines[1] == "X 5 8"


@pytest.mark.octave
@skipif_no_octave
def test_octave_bp_vector():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "fprintf('%d\\n', ref.bp);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    bp_vals = [int(x) for x in result.stdout.strip().splitlines()]
    assert bp_vals == [100, 200, 300, 400, 500, 100, 200, 300]


@pytest.mark.octave
@skipif_no_octave
def test_octave_cache_roundtrip(tmp_path):
    cache_path = str(tmp_path / "ref_cache.mat")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"statgen.save_reference_cache(ref, '{cache_path}'); "
        f"ref2 = statgen.load_reference_cache('{cache_path}'); "
        f"fprintf('%d\\n', ref2.num_snp); "
        f"for i = 1:numel(ref2.shards); "
        f"  fprintf('%s %s\\n', ref2.shards{{i}}.label, ref2.shards{{i}}.checksum); "
        f"end; "
        f"cache_ok = 1; "
        f"for i = 1:numel(ref.shards); "
        f"  s1 = ref.shards{{i}}; s2 = ref2.shards{{i}}; "
        f"  cache_ok = cache_ok && isequal(s1.chr, s2.chr) && isequal(s1.snp, s2.snp) "
        f"           && isequal(s1.bp, s2.bp) && isequal(s1.a1, s2.a1) && isequal(s1.a2, s2.a2) "
        f"           && isequal(s1.a1_hash64, s2.a1_hash64) && isequal(s1.a2_hash64, s2.a2_hash64); "
        f"end; "
        f"display_ok = 0; try; txt = evalc('ref2'); display_ok = ~isempty(strfind(txt, 'ReferencePanel')); catch; end; "
        f"shard_display_ok = 0; try; txt = evalc('ref2.shards{{1}}'); shard_display_ok = ~isempty(strfind(txt, 'ReferenceShard')); catch; end; "
        f"fprintf('%d\\n', cache_ok); "
        f"fprintf('%d\\n', ref2.validate_checksums()); "
        f"fprintf('%d\\n', display_ok); "
        f"fprintf('%d\\n', shard_display_ok); "
        f"s_cache = load('{cache_path}'); "
        f"fprintf('%d %d %d %d %d %d %d\\n', isfield(s_cache, 'metadata'), isfield(s_cache, 'chr'), isfield(s_cache, 'snp'), isfield(s_cache, 'a1'), isfield(s_cache, 'snp_text_by_shard'), isfield(s_cache, 'a1_hash64'), isfield(s_cache, 'a2_hash64')); "
        f"fprintf('%d %d %d\\n', numel(s_cache.snp_text_by_shard), s_cache.metadata.shard_start0(1), s_cache.metadata.shard_stop0(end));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "8"
    assert lines[1] == f"1 {CHR1_CHECKSUM}"
    assert lines[2] == f"X {CHRX_CHECKSUM}"
    assert lines[3] == "1"
    assert lines[4] == "1"
    assert lines[5] == "1"
    assert lines[6] == "1"
    assert lines[7] == "1 0 0 0 1 1 1"
    assert lines[8] == "2 0 8"


@pytest.mark.octave
@skipif_no_octave
def test_octave_reference_cache_trusts_checksum_until_explicit_validation(tmp_path):
    cache_path = str(tmp_path / "ref_cache.mat")
    bad_cache_path = str(tmp_path / "ref_cache_bad_payload.mat")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"statgen.save_reference_cache(ref, '{cache_path}'); "
        f"L = load('{cache_path}'); "
        "metadata = L.metadata; bp = L.bp; snp_text_by_shard = L.snp_text_by_shard; a1_text_by_shard = L.a1_text_by_shard; a2_text_by_shard = L.a2_text_by_shard; a1_hash64 = L.a1_hash64; a2_hash64 = L.a2_hash64; "
        "a1_text_by_shard{1} = strrep(a1_text_by_shard{1}, 'A', 'T'); "
        f"save('{bad_cache_path}', 'metadata', 'bp', 'snp_text_by_shard', 'a1_text_by_shard', 'a2_text_by_shard', 'a1_hash64', 'a2_hash64'); "
        f"ref2 = statgen.load_reference_cache('{bad_cache_path}'); "
        "ok1 = ref2.num_snp == 8; "
        "ok2 = 0; try; ref2.validate_checksums(); catch; ok2 = 1; end; "
        "fprintf('%d\\n', ok1); "
        "fprintf('%d\\n', ok2); "
        "fprintf('%s\\n', ref2.shards{1}.checksum);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
    assert lines[1] == "1"
    assert lines[2] == CHR1_CHECKSUM


@pytest.mark.octave
@skipif_no_octave
def test_octave_reference_cache_accepts_format_option(tmp_path):
    cache_path = str(tmp_path / "ref_cache_v5.mat")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"statgen.save_reference_cache(ref, '{cache_path}', 'format', 'v5'); "
        f"ref2 = statgen.load_reference_cache('{cache_path}'); "
        f"fprintf('%d\\n', ref2.num_snp);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "8"


@pytest.mark.octave
@skipif_no_octave
def test_octave_reference_cache_rejects_zero_shards(tmp_path):
    cache_path = str(tmp_path / "ref_cache.mat")
    bad_cache_path = str(tmp_path / "ref_cache_zero_shards.mat")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"statgen.save_reference_cache(ref, '{cache_path}'); "
        f"L = load('{cache_path}'); "
        "L.metadata.n_shards = 0; "
        "L.metadata.shard_labels = {}; "
        "L.metadata.shard_checksums = {}; "
        "L.metadata.shard_start0 = []; "
        "L.metadata.shard_stop0 = []; "
        f"save('{bad_cache_path}', '-struct', 'L'); "
        f"ok = 0; try; statgen.load_reference_cache('{bad_cache_path}'); catch ME; ok = ~isempty(strfind(ME.message, 'n_shards must be positive')); end; "
        "fprintf('%d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_is_object_compatible():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "ref2 = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "fprintf('%d\\n', ref.is_object_compatible(ref2));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_select_shards():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "sub = ref.select_shards({'X'}); "
        "fprintf('%d\\n', sub.num_snp); "
        "fprintf('%s\\n', sub.shards{1}.label);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "3"
    assert lines[1] == "X"
