"""Phase 0 acceptance tests: skeleton, fixtures, and Octave harness smoke test."""

import gzip
import json
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import numpy as np
import pytest
from scipy.io import loadmat

from tests.conftest import FIXTURES_DIR, matlab_data_lines, skipif_no_octave, run_octave


def _python_project_version() -> str:
    pyproject = Path(__file__).parents[1] / "python" / "pyproject.toml"
    return tomllib.loads(pyproject.read_text(encoding="utf-8"))["project"]["version"]


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
    "sumstats/traits_complete.tsv.gz",
    "ld/python/ld_manifest.json",
    "ld/python/reference_chr1.bim",
    "ld/python/reference_chrX.bim",
    "ld/python/ld_chr1.npz",
    "ld/python/ld_chrX_female.npz",
    "ld/python/ld_chrX_male.npz",
    "ld/python/ld_chrX_combined.npz",
    "ld/matlab/ld_manifest.json",
    "ld/matlab/reference_chr1.bim",
    "ld/matlab/reference_chrX.bim",
    "ld/matlab/ld_chr1.mat",
    "ld/matlab/ld_chrX_female.mat",
    "ld/matlab/ld_chrX_male.mat",
    "ld/matlab/ld_chrX_combined.mat",
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


def test_ld_python_manifest():
    meta = json.loads((FIXTURES_DIR / "ld/python/ld_manifest.json").read_text())
    assert meta["object_type"] == "ld_panel_manifest"
    assert meta["runtime_format"] == "python_npz_csc32"
    assert [(s["chr"], s["sex"]) for s in meta["shards"]] == [
        ("1", None),
        ("X", "female"),
        ("X", "male"),
        ("X", "combined"),
    ]
    assert all(len(s["reference_checksum"]) == 32 for s in meta["shards"])
    assert all(len(s["file_md5"]) == 32 for s in meta["shards"])
    assert [s["reference_bim"] for s in meta["shards"]] == [
        "reference_chr1.bim",
        "reference_chrX.bim",
        "reference_chrX.bim",
        "reference_chrX.bim",
    ]


def test_ld_chr1_npz_payload():
    with np.load(FIXTURES_DIR / "ld/python/ld_chr1.npz", allow_pickle=False) as z:
        assert set(z.files) == {"data", "indices", "indptr", "shape", "a1freq", "metadata"}
        assert z["data"].dtype == np.float32
        assert z["indices"].dtype == np.int32
        assert z["indptr"].dtype == np.int32
        assert z["shape"].tolist() == [5, 5]
        assert len(z["a1freq"]) == 5
        assert z["a1freq"].dtype == np.float32
        assert z["indptr"][-1] == len(z["data"]) == len(z["indices"]) == 13
        meta = json.loads(z["metadata"].tobytes().decode())
    assert meta["object_type"] == "ld_shard"
    assert meta["schema_version"] == "1.0"
    assert meta["format"] == "statgen_ld_npz_csc32"
    assert meta["chr"] == "1"
    assert meta["sex"] is None
    assert meta["num_snp"] == 5
    assert meta["nnz"] == 13
    assert meta["num_monomorphic_snps"] == 0
    assert meta["reference_bim"] == "reference_chr1.bim"


def test_ld_chrX_npz_metadata():
    with np.load(FIXTURES_DIR / "ld/python/ld_chrX_female.npz", allow_pickle=False) as z:
        meta = json.loads(z["metadata"].tobytes().decode())
        assert z["shape"].tolist() == [3, 3]
        assert z["indptr"][-1] == len(z["data"]) == len(z["indices"]) == 7
    assert meta["chr"] == "X"
    assert meta["sex"] == "female"
    assert meta["diagonal"] == "explicit_unit"
    assert meta["reference_bim"] == "reference_chrX.bim"


def test_ld_matlab_fixture_payload():
    mat = loadmat(FIXTURES_DIR / "ld/matlab/ld_chr1.mat", squeeze_me=True)
    assert set(["ld_r", "a1freq", "metadata"]).issubset(mat)
    assert mat["ld_r"].shape == (5, 5)
    assert mat["ld_r"].nnz == 13
    assert mat["a1freq"].dtype == np.float64
    assert len(np.atleast_1d(mat["a1freq"])) == 5
    assert str(mat["metadata"]["reference_bim"]) == "reference_chr1.bim"


def test_plink_bed_magic():
    bed = (FIXTURES_DIR / "genotype/sharded/1.bed").read_bytes()
    assert bed[:3] == b"\x6c\x1b\x01"


# ---------------------------------------------------------------------------
# Python package importable
# ---------------------------------------------------------------------------

def test_statgen_importable():
    import statgen
    assert statgen.__version__ == _python_project_version()


def test_matlab_data_lines_skips_warning_stack_frames():
    stdout = (
        "[\bWarning: load_sumstats: zvec is absent]\b \n"
        "[> In <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback"
        "('warn_one_optional_', '/tmp/file.m', 10)\" style=\"font-weight:bold\">"
        "warn_one_optional_</a> (<a href=\"matlab: opentoline('/tmp/file.m',10,0)\">"
        "line 10</a>)\n"
        "8\n"
        "0\n"
    )
    assert matlab_data_lines(stdout) == ["8", "0"]


# ---------------------------------------------------------------------------
# Octave smoke test
# ---------------------------------------------------------------------------

@pytest.mark.octave
@skipif_no_octave
def test_octave_version_smoke():
    import statgen
    result = run_octave("v = statgen.version(); disp(v);")
    assert result.returncode == 0, f"Octave stderr:\n{result.stderr}"
    assert statgen.__version__ in result.stdout
