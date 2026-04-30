"""Phase 0 acceptance tests: skeleton, fixtures, and Octave harness smoke test."""

import gzip
from pathlib import Path

import pytest

from tests.conftest import FIXTURES_DIR, skipif_no_octave, run_octave


# ---------------------------------------------------------------------------
# Fixture presence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel", [
    "reference/sharded/1.bim",
    "reference/sharded/X.bim",
    "reference/nonsharded/all.bim",
    "annotations/anno1.bed",
    "annotations/anno2.bed",
    "sumstats/traits.tsv.gz",
    "ld/1/metadata.json",
    "ld/1/ld_idx1.i32",
    "ld/1/ld_idx2.i32",
    "ld/1/ld_r.f32",
    "ld/1/mafvec.f32",
    "ld/X/metadata.json",
    "ld/X/combined/metadata.json",
    "ld/X/combined/ld_idx1.i32",
    "ld/X/male/metadata.json",
    "ld/X/female/metadata.json",
    "genotype/sharded/1.bim",
    "genotype/sharded/1.fam",
    "genotype/sharded/1.bed",
    "genotype/sharded/X.bim",
    "genotype/sharded/X.fam",
    "genotype/sharded/X.bed",
])
def test_fixture_exists(rel):
    assert (FIXTURES_DIR / rel).is_file(), f"missing fixture: {rel}"


# ---------------------------------------------------------------------------
# Fixture sanity
# ---------------------------------------------------------------------------

def test_bim_chr1_row_count():
    lines = (FIXTURES_DIR / "reference/sharded/1.bim").read_text().splitlines()
    assert len(lines) == 5


def test_bim_chrX_row_count():
    lines = (FIXTURES_DIR / "reference/sharded/X.bim").read_text().splitlines()
    assert len(lines) == 3


def test_nonsharded_bim_contains_both_chromosomes():
    lines = (FIXTURES_DIR / "reference/nonsharded/all.bim").read_text().splitlines()
    chrs = {ln.split("\t")[0] for ln in lines}
    assert chrs == {"1", "X"}


def test_sumstats_readable():
    with gzip.open(FIXTURES_DIR / "sumstats/traits.tsv.gz", "rt") as f:
        header = f.readline().strip().split("\t")
    assert "chr" in header and "z" in header and "p" in header


def test_ld_chr1_binary_lengths():
    import numpy as np
    base = FIXTURES_DIR / "ld/1"
    idx1 = np.fromfile(base / "ld_idx1.i32", dtype="<i4")
    idx2 = np.fromfile(base / "ld_idx2.i32", dtype="<i4")
    r    = np.fromfile(base / "ld_r.f32",    dtype="<f4")
    maf  = np.fromfile(base / "mafvec.f32",  dtype="<f4")
    assert len(idx1) == len(idx2) == len(r) == 4   # num_ld
    assert len(maf) == 5                            # num_snp
    assert (idx1 < idx2).all()


def test_ld_chr1_metadata():
    import json
    meta = json.loads((FIXTURES_DIR / "ld/1/metadata.json").read_text())
    assert meta["object_type"] == "ld_shard"
    assert meta["chr"] == "1"
    assert meta["num_snp"] == 5
    assert meta["num_ld"] == 4
    assert meta["byte_order"] == "little_endian"
    assert len(meta["reference_checksum"]) == 32   # MD5 hex


def test_ld_chrX_shard_group_metadata():
    import json
    meta = json.loads((FIXTURES_DIR / "ld/X/metadata.json").read_text())
    assert meta["object_type"] == "ld_shard_group"
    assert set(meta["sex_shards"]) == {"combined", "male", "female"}


def test_plink_bed_magic():
    bed = (FIXTURES_DIR / "genotype/sharded/1.bed").read_bytes()
    assert bed[:3] == b"\x6c\x1b\x01"


# ---------------------------------------------------------------------------
# Python package importable
# ---------------------------------------------------------------------------

def test_statgen_importable():
    import statgen
    assert statgen.__version__ == "0.1.0"


# ---------------------------------------------------------------------------
# Octave smoke test
# ---------------------------------------------------------------------------

@pytest.mark.octave
@skipif_no_octave
def test_octave_version_smoke():
    result = run_octave("v = statgen.version(); disp(v);")
    assert result.returncode == 0, f"Octave stderr:\n{result.stderr}"
    assert "0.1" in result.stdout
