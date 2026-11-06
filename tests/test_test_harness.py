"""Regression checks for runtime-specific pytest setup."""

from pathlib import Path

import pytest

pytest_plugins = ["pytester"]


@pytest.mark.parametrize("selection,prepares_r", [("octave", False), ("not octave and not r", False), ("r", True)])
def test_r_fixture_preparation_only_runs_when_requested(pytester, selection, prepares_r):
    # A fresh checkout with no derived R files must still run unrelated tests.
    harness = Path(__file__).with_name("conftest.py").read_text()
    pytester.makeconftest(
        harness
        + '''

def pytest_configure(config):
    original_run = subprocess.run
    def run(args, **kwargs):
        if args == ["make", "prepare-r-fixtures"]:
            print("PREPARING_R_FIXTURES")
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        return original_run(args, **kwargs)
    subprocess.run = run
'''
    )
    pytester.makepyfile(
        '''
import pytest

@pytest.mark.octave
def test_matlab():
    pass

def test_python():
    pass

@pytest.mark.r
def test_r(prepare_r_extdata_fixtures):
    pass
'''
    )
    result = pytester.runpytest_subprocess("-m", selection, "-s", "-q")
    result.assert_outcomes(passed=1, deselected=2)
    assert ("PREPARING_R_FIXTURES" in result.stdout.str()) == prepares_r
