import shutil
import subprocess
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
REPO_ROOT = Path(__file__).parent.parent
MATLAB_DIR = REPO_ROOT / "matlab"


@pytest.fixture(scope="session")
def octave_available():
    return shutil.which("octave") is not None


skipif_no_octave = pytest.mark.skipif(
    shutil.which("octave") is None,
    reason="Octave not installed",
)


def run_octave(expr: str, timeout: int = 30) -> subprocess.CompletedProcess:
    """Run an Octave expression with matlab/ on the path. Returns CompletedProcess."""
    return subprocess.run(
        ["octave", "--no-gui", "--quiet", "--eval",
         f"addpath('{MATLAB_DIR}'); {expr}"],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
