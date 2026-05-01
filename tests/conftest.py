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
_MATLAB_ENGINE = None
_MATLAB_ENGINE_ERROR = None


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
