# Development Notes

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
