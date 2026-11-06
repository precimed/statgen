import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
REPO_ROOT = Path(__file__).parent.parent
MATLAB_DIR = REPO_ROOT / "matlab"
R_PACKAGE_DIR = REPO_ROOT / "R-package"
R_EXTDATA_COPIES = {
    "allele_swaps/reference.bim": "allele_swaps/reference.bim",
    "allele_swaps/sumstats.tsv.gz": "allele_swaps/sumstats.tsv.gz",
    "allele_swaps/expected.json": "allele_swaps/expected.json",
    "allele_swaps/genotype/1.bed": "allele_swaps/genotype/1.bed",
    "allele_swaps/genotype/1.bim": "allele_swaps/genotype/1.bim",
    "allele_swaps/genotype/1.fam": "allele_swaps/genotype/1.fam",
    "allele_swaps/genotype/X.bed": "allele_swaps/genotype/X.bed",
    "allele_swaps/genotype/X.bim": "allele_swaps/genotype/X.bim",
    "allele_swaps/genotype/X.fam": "allele_swaps/genotype/X.fam",
    "allele_swaps/genotype/X.ploidy": "allele_swaps/genotype/X.ploidy",
    "reference_chr1.bim": "reference/sharded/1.bim",
    "reference_chrX.bim": "reference/sharded/X.bim",
    "traits.tsv.gz": "sumstats/traits.tsv.gz",
    "traits_complete.tsv.gz": "sumstats/traits_complete.tsv.gz",
    "anno1.bed": "annotations/anno1.bed",
    "anno2.bed": "annotations/anno2.bed",
    "grouped.annot": "annotations/grouped.annot",
    "genotype_1.bed": "genotype/sharded/1.bed",
    "genotype_1.bim": "genotype/sharded/1.bim",
    "genotype_1.fam": "genotype/sharded/1.fam",
    "genotype_X.bed": "genotype/sharded/X.bed",
    "genotype_X.bim": "genotype/sharded/X.bim",
    "genotype_X.fam": "genotype/sharded/X.fam",
    "genotype_X.ploidy": "genotype/sharded/X.ploidy",
    "ld/reference_cache.npz": "ld/python/reference_cache.npz",
    "ld/reference_chr1.bim": "ld/python/reference_chr1.bim",
    "ld/reference_chrX.bim": "ld/python/reference_chrX.bim",
    "ld/ld_chr1.npz": "ld/python/ld_chr1.npz",
    "ld/ld_chrX_female.npz": "ld/python/ld_chrX_female.npz",
    "ld/ld_chrX_male.npz": "ld/python/ld_chrX_male.npz",
    "ld/ld_chrX_combined.npz": "ld/python/ld_chrX_combined.npz",
}

# Set STATGEN_MATLAB=1 to run octave-marked tests via native MATLAB instead.
_USE_MATLAB = os.environ.get("STATGEN_MATLAB", "0") == "1"
_ENGINE = "matlab" if _USE_MATLAB else "octave"
_ENGINE_AVAILABLE = shutil.which(_ENGINE) is not None
_RSCRIPT_AVAILABLE = shutil.which("Rscript") is not None
_MATLAB_ENGINE = None
_MATLAB_ENGINE_ERROR = None


@pytest.fixture(scope="session")
def octave_available():
    return _ENGINE_AVAILABLE


@pytest.fixture(scope="session")
def rscript_available():
    return _RSCRIPT_AVAILABLE


skipif_no_octave = pytest.mark.skipif(
    not _ENGINE_AVAILABLE,
    reason=f"{_ENGINE} not installed",
)

skipif_no_rscript = pytest.mark.skipif(
    not _RSCRIPT_AVAILABLE,
    reason="Rscript not installed",
)

skipif_no_matlab = pytest.mark.skipif(
    shutil.which("matlab") is None,
    reason="MATLAB not installed",
)

skipif_matlab_engine = pytest.mark.skipif(
    _USE_MATLAB,
    reason="Octave-specific test",
)


def _get_matlab_engine():
    global _MATLAB_ENGINE, _MATLAB_ENGINE_ERROR
    if _MATLAB_ENGINE is not None or _MATLAB_ENGINE_ERROR is not None:
        return _MATLAB_ENGINE

    try:
        import matlab.engine
    except Exception as exc:  # pragma: no cover - environment dependent
        _MATLAB_ENGINE_ERROR = (
            "MATLAB engine for Python is required when STATGEN_MATLAB=1. "
            f"Import failed: {exc}"
        )
        return None

    try:
        eng = matlab.engine.start_matlab("-nosplash -nodesktop")
        eng.addpath(str(MATLAB_DIR), nargout=0)
        _MATLAB_ENGINE = eng
        return _MATLAB_ENGINE
    except Exception as exc:  # pragma: no cover - environment dependent
        _MATLAB_ENGINE_ERROR = f"Failed to start MATLAB engine: {exc}"
        return None


def _r_ld_fixture_is_prepared() -> bool:
    r_manifest_path = R_PACKAGE_DIR / "inst/extdata/ld/ld_manifest.json"
    r_reference_cache = R_PACKAGE_DIR / "inst/extdata/ld/reference_cache.rds"
    canonical_manifest_path = FIXTURES_DIR / "ld/python/ld_manifest.json"
    if not r_manifest_path.is_file() or not r_reference_cache.is_file():
        return False
    try:
        r_manifest = json.loads(r_manifest_path.read_text())
        canonical_manifest = json.loads(canonical_manifest_path.read_text())
    except json.JSONDecodeError:
        return False
    return (
        r_manifest.pop("r_reference_cache", None) == "reference_cache.rds"
        and isinstance(r_manifest.pop("r_reference_cache_md5", None), str)
        and r_manifest == canonical_manifest
    )


@pytest.fixture(scope="session", autouse=True)
def _close_matlab_engine_at_end():
    yield
    global _MATLAB_ENGINE
    if _MATLAB_ENGINE is not None:
        try:
            _MATLAB_ENGINE.quit()
        except Exception:
            pass
        _MATLAB_ENGINE = None


@pytest.fixture(scope="session", autouse=True)
def _prepare_r_extdata_fixtures():
    needs_prepare = False
    for r_name, canonical_rel in R_EXTDATA_COPIES.items():
        r_fixture = R_PACKAGE_DIR / f"inst/extdata/{r_name}"
        canonical = FIXTURES_DIR / canonical_rel
        if not r_fixture.is_file() or r_fixture.read_bytes() != canonical.read_bytes():
            needs_prepare = True
            break
    if not _r_ld_fixture_is_prepared():
        needs_prepare = True
    if needs_prepare:
        result = subprocess.run(
            ["make", "prepare-r-fixtures"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0, result.stderr


def run_octave(expr: str, timeout: int = 30) -> subprocess.CompletedProcess:
    """Run an expression via Octave or MATLAB (controlled by STATGEN_MATLAB=1)."""
    if _USE_MATLAB:
        eng = _get_matlab_engine()
        if eng is None:
            return subprocess.CompletedProcess(
                args=["matlab-engine", expr],
                returncode=1,
                stdout="",
                stderr=_MATLAB_ENGINE_ERROR or "MATLAB engine initialization failed",
            )
        wrapped = (
            "try; "
            + expr
            + "; "
            + "catch ME; fprintf(2, '%s\\n', getReport(ME, 'extended')); rethrow(ME); end"
        )
        try:
            out = eng.evalc(wrapped, nargout=1)
            return subprocess.CompletedProcess(
                args=["matlab-engine", expr],
                returncode=0,
                stdout=out,
                stderr="",
            )
        except Exception as exc:
            return subprocess.CompletedProcess(
                args=["matlab-engine", expr],
                returncode=1,
                stdout="",
                stderr=str(exc),
            )
    else:
        cmd = ["octave", "--no-gui", "--quiet", "--eval",
               f"addpath('{MATLAB_DIR}'); {expr}"]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def run_rscript(expr: str, timeout: int = 30, env=None) -> subprocess.CompletedProcess:
    """Run an expression via Rscript."""
    cmd = ["Rscript", "-e", expr]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)


def matlab_data_lines(stdout: str) -> list[str]:
    """Return non-warning stdout lines from Octave/MATLAB test snippets."""
    out = []
    in_warning_block = False
    for raw in stdout.splitlines():
        line = raw.replace("\b", "").strip()
        if not line:
            continue
        unwrapped = line.lstrip("[")
        if line.startswith("[Warning:") or line.startswith("Warning:"):
            in_warning_block = line.startswith("[Warning:") and not line.endswith("]")
            continue
        if in_warning_block:
            if line.endswith("]"):
                in_warning_block = False
            continue
        if unwrapped.startswith("> In ") or unwrapped.startswith("In "):
            continue
        if line == "]":
            continue
        out.append(line)
    return out


GENOTYPE_CHR1_CALLS = np.array(
    [
        [2, -1, 1, 0],
        [0, 1, 2, -1],
        [1, 1, 0, 2],
        [2, 0, -1, 1],
        [-1, 2, 0, 1],
    ],
    dtype=np.int8,
)

GENOTYPE_SAMPLES = [
    ("FAM1", "IND1", 0, 0, 1, -9),
    ("FAM1", "IND2", 0, 0, 2, -9),
    ("FAM2", "IND3", 0, 0, 1, -9),
    ("FAM2", "IND4", 0, 0, 2, -9),
]

GENOTYPE_CHR1_BIM = [
    ("1", "rs1001", 0, 100, "A", "G"),
    ("1", "rs1002", 0, 200, "C", "T"),
    ("1", "rs1003", 0, 300, "A", "C"),
    ("1", "rs1004", 0, 400, "G", "A"),
    ("1", "rs1005", 0, 500, "T", "C"),
]

GENOTYPE_CHR2_BIM = [
    ("2", "rs2001", 0, 100, "A", "G"),
    ("2", "rs2002", 0, 200, "C", "T"),
]


def write_plink_bed(path: Path, *, num_snp: int, num_sample: int) -> None:
    bytes_per_snp = (num_sample + 3) // 4
    path.write_bytes(b"\x6c\x1b\x01" + b"\x00" * (num_snp * bytes_per_snp))


def write_plink_bed_calls(path: Path, calls) -> None:
    calls = np.asarray(calls, dtype=np.int8)
    code = {
        2: 0b00,
        -1: 0b01,
        1: 0b10,
        0: 0b11,
    }
    bytes_per_snp = (calls.shape[1] + 3) // 4
    payload = bytearray()
    for row in calls:
        row_bytes = [0] * bytes_per_snp
        for j, value in enumerate(row.tolist()):
            row_bytes[j // 4] |= code[int(value)] << (2 * (j % 4))
        payload.extend(row_bytes)
    path.write_bytes(b"\x6c\x1b\x01" + bytes(payload))


def write_plink_bim(path: Path, rows) -> None:
    path.write_text("".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n" for r in rows))


def write_plink_fam(path: Path, rows) -> None:
    path.write_text("".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n" for r in rows))


def copy_genotype_shard_files(dst: Path, label: str, *, ploidy: bool = True) -> None:
    for suffix in (".bim", ".fam", ".bed"):
        (dst / f"{label}{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/{label}{suffix}").read_bytes()
        )
    src_ploidy = FIXTURES_DIR / f"genotype/sharded/{label}.ploidy"
    if ploidy and src_ploidy.exists():
        (dst / f"{label}.ploidy").write_bytes(src_ploidy.read_bytes())


def copy_sharded_genotype(dst: Path) -> None:
    for label in ("1", "X"):
        copy_genotype_shard_files(dst, label)
