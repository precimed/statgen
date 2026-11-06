# Development Notes

## Setup

Create and activate the conda environment:

```sh
conda create -n statgen python=3.11 numpy scipy pytest -y
conda activate statgen
pip install -e python/
```

Run the full test suite. Octave tests are skipped automatically when Octave is
not installed:

```sh
make test
```

Regenerate committed test fixtures after a format change:

```sh
make fixtures
```

## MATLAB tests from a Python 3.9 environment

For a machine whose installed MATLAB Engine uses Python 3.9, use an existing
Python 3.9 environment with `pytest`, `numpy`, `scipy`, `pandas`, `tomli`, and
MATLAB Engine installed. Run from the repository root with `matlab` on `PATH`:

```sh
STATGEN_MATLAB=1 python -m pytest tests/ -m octave
```

This is a development setup for the MATLAB test harness, not a Python 3.9
support commitment for the package. The package requirement remains
`requires-python = ">=3.10"`; the editable-install setup above uses a supported
Python version.

Pytest imports test modules before applying the `octave` marker filter, so
collection errors can come from Python imports before MATLAB starts. The
current compatibility adjustments allow those imports and the Python helpers
used by MATLAB tests to run in this environment. `tomli` supplies the TOML
parser on Python versions without `tomllib`.

## R package documentation and checks

R package fixture storage and synchronization are defined in
[spec/testing.md](spec/testing.md); CRAN-facing runtime and test requirements
are defined in [spec/R.md](spec/R.md#cran-facing-tests).

Commit bundled fixture copies under `R-package/inst/extdata/` so native R
package builds, vignettes, and CRAN checks work from a fresh checkout without
a preparation step. When canonical fixtures change, run
`make prepare-r-fixtures` and commit the refreshed copies with the source
changes. The repository's R test/build and documentation targets also refresh
these copies automatically.

## Releases (brief)

- Manual validation and optional TestPyPI publish are run via GitHub Actions
  `Release` workflow dispatch.
- PyPI publish and GitHub release are triggered by pushing a new `v*` tag.

## Version alignment checklist

Before release, ensure version values are aligned in:

- `python/pyproject.toml` (`project.version`)
- `python/statgen/__init__.py` (`__version__`)
- `matlab/+statgen/version.m` (MATLAB/Octave `statgen.version()` return value)

Note: MATLAB currently returns a short form (for example `0.1`). Keep it
consistent with the intended release series when bumping versions.
