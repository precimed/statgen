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
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, run_octave, skipif_no_octave

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


# ---------------------------------------------------------------------------
# Python: is_object_compatible
# ---------------------------------------------------------------------------

class _FakeShard:
    def __init__(self, num_snp, checksum=None):
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
    fake = _FakePanel([_FakeShard(5)])  # only one shard instead of two
    assert p.is_object_compatible(fake) is False


def test_compat_row_count_mismatch():
    p = load_reference(SHARDED)
    # chr1 has 5 SNPs, chrX has 3 — swap the counts
    fake = _FakePanel([_FakeShard(3), _FakeShard(5)])
    assert p.is_object_compatible(fake) is False


def test_compat_checksum_mismatch():
    p = load_reference(SHARDED)
    fake = _FakePanel([
        _FakeShard(5, checksum="deadbeef" * 4),
        _FakeShard(3, checksum=CHRX_CHECKSUM),
    ])
    assert p.is_object_compatible(fake) is False


def test_compat_no_checksum_passes():
    p = load_reference(SHARDED)
    # Shards with correct num_snp but no checksum attribute → compatible
    fake = _FakePanel([_FakeShard(5), _FakeShard(3)])
    assert p.is_object_compatible(fake) is True


def test_compat_does_not_raise(caplog):
    import logging
    p = load_reference(SHARDED)
    fake = _FakePanel([_FakeShard(99), _FakeShard(3)])
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


# ---------------------------------------------------------------------------
# Python: BIM validation
# ---------------------------------------------------------------------------

def test_bim_wrong_column_count(tmp_path):
    bad = tmp_path / "bad.bim"
    bad.write_text("1\trs1\t0\t100\tA\n")  # 5 columns
    with pytest.raises(ValueError, match="6 tab-separated"):
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
    with pytest.raises(ValueError, match="sorted by canonical chromosome"):
        load_reference(bad)


def test_bim_unsorted_chr_fails(tmp_path):
    bad = tmp_path / "unsorted_chr.bim"
    bad.write_text(
        "X\trsX\t0\t100\tA\tG\n"
        "1\trs1\t0\t200\tC\tT\n"
    )
    with pytest.raises(ValueError, match="sorted by canonical chromosome"):
        load_reference(bad)


def test_bim_noncanonical_chr_fails(tmp_path):
    bad = tmp_path / "noncanonical_chr.bim"
    bad.write_text("chr1\trs1\t0\t100\tA\tG\n")
    with pytest.raises(ValueError, match="canonical labels"):
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


def test_cache_load_skips_source_sort_validation(monkeypatch, tmp_path):
    import statgen.reference as ref_mod

    panel = load_reference(SHARDED)
    cache = tmp_path / "ref.npz"
    save_reference_cache(panel, cache)

    def _fail(*_args, **_kwargs):
        raise AssertionError("source sort validation should not run on cache load")

    monkeypatch.setattr(ref_mod, "_validate_reference_sort_order", _fail)
    loaded = load_reference_cache(cache)
    assert loaded.num_snp == panel.num_snp


def test_bim_sharded_no_match(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_reference(str(tmp_path / "@.bim"))


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


@skipif_no_octave
def test_octave_sharded_num_snp():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "printf('%d\\n', ref.num_snp);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "8"


@skipif_no_octave
def test_octave_sharded_shard_labels():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "for i = 1:numel(ref.shards); printf('%s\\n', ref.shards{i}.label); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines == ["1", "X"]


@skipif_no_octave
def test_octave_checksums_match_python():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "for i = 1:numel(ref.shards); printf('%s\\n', ref.shards{i}.checksum); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == CHR1_CHECKSUM
    assert lines[1] == CHRX_CHECKSUM


@skipif_no_octave
def test_octave_nonsharded_split():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/nonsharded/all.bim']); "
        "printf('%d\\n', numel(ref.shards)); "
        "for i = 1:numel(ref.shards); printf('%s %d\\n', ref.shards{i}.label, ref.shards{i}.num_snp); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "2"
    assert lines[1] == "1 5"
    assert lines[2] == "X 3"


@skipif_no_octave
def test_octave_shard_offsets():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "off = ref.shard_offsets; "
        "for i = 1:numel(off); "
        "  printf('%s %d %d\\n', off(i).shard_label, off(i).start0, off(i).stop0); "
        "end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1 0 5"
    assert lines[1] == "X 5 8"


@skipif_no_octave
def test_octave_bp_vector():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "printf('%d\\n', ref.bp);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    bp_vals = [int(x) for x in result.stdout.strip().splitlines()]
    assert bp_vals == [100, 200, 300, 400, 500, 100, 200, 300]


@skipif_no_octave
def test_octave_cache_roundtrip(tmp_path):
    cache_path = str(tmp_path / "ref_cache.mat")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"statgen.save_reference_cache(ref, '{cache_path}'); "
        f"ref2 = statgen.load_reference_cache('{cache_path}'); "
        f"printf('%d\\n', ref2.num_snp); "
        f"for i = 1:numel(ref2.shards); "
        f"  printf('%s %s\\n', ref2.shards{{i}}.label, ref2.shards{{i}}.checksum); "
        f"end; "
        f"ok = 1; "
        f"for i = 1:numel(ref.shards); "
        f"  s1 = ref.shards{{i}}; s2 = ref2.shards{{i}}; "
        f"  ok = ok && isequal(s1.chr, s2.chr) && isequal(s1.snp, s2.snp) "
        f"           && isequal(s1.bp, s2.bp) && isequal(s1.a1, s2.a1) && isequal(s1.a2, s2.a2); "
        f"end; "
        f"printf('%d\\n', ok); "
        f"s = load('{cache_path}'); "
        f"printf('%d\\n', isfield(s.cache_shards, 'cm'));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "8"
    assert lines[1] == f"1 {CHR1_CHECKSUM}"
    assert lines[2] == f"X {CHRX_CHECKSUM}"
    assert lines[3] == "1"
    assert lines[4] == "0"


@skipif_no_octave
def test_octave_is_object_compatible():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "ref2 = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "printf('%d\\n', ref.is_object_compatible(ref2));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
