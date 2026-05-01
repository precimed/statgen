import os
import shutil
import subprocess
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
REPO_ROOT = Path(__file__).parent.parent
MATLAB_DIR = REPO_ROOT / "matlab"

# Set STATGEN_MATLAB=1 to run octave-marked tests via native MATLAB instead.
_USE_MATLAB = os.environ.get("STATGEN_MATLAB", "0") == "1"
_ENGINE = "matlab" if _USE_MATLAB else "octave"
_ENGINE_AVAILABLE = shutil.which(_ENGINE) is not None


@pytest.fixture(scope="session")
def octave_available():
    return _ENGINE_AVAILABLE


skipif_no_octave = pytest.mark.skipif(
    not _ENGINE_AVAILABLE,
    reason=f"{_ENGINE} not installed",
)

skipif_no_matlab = pytest.mark.skipif(
    shutil.which("matlab") is None,
    reason="MATLAB not installed",
)


def run_octave(expr: str, timeout: int = 30) -> subprocess.CompletedProcess:
    """Run an expression via Octave or MATLAB (controlled by STATGEN_MATLAB=1)."""
    if _USE_MATLAB:
        cmd = ["matlab", "-batch", f"addpath('{MATLAB_DIR}'); {expr}"]
    else:
        cmd = ["octave", "--no-gui", "--quiet", "--eval",
               f"addpath('{MATLAB_DIR}'); {expr}"]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
