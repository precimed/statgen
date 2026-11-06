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

## R package documentation and checks

Run `make prepare-r-fixtures` before directly building the R package or its
vignette from a checkout. This prepares bundled copies of the canonical test
fixtures. The repository's R test/build and documentation targets already run
this preparation step.

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
