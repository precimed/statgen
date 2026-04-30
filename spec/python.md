# Python implementation conventions

## Package layout

The Python package lives under `python/statgen/` with a `pyproject.toml` at
`python/pyproject.toml`. It is installable via `pip install -e python/` from
the repository root. There is no `sys.path` manipulation in test helpers or
conftest; the editable install is the only supported import mechanism.

Module boundaries mirror the object families:

```
python/
  pyproject.toml
  statgen/
    __init__.py
    reference.py
    sumstats.py
    annotations.py
    ld.py
    genotype.py
    _utils.py        # shared internals, not public API
```

## Naming conventions

- Modules, functions, and variables: `snake_case`.
- Classes: `PascalCase` matching the spec object names (`ReferenceShard`,
  `ReferencePanel`, `LDPanel`, etc.).
- Public loader functions are module-level free functions: `load_reference`,
  `load_ld`, `load_sumstats`, `load_annotations`.
- Private helpers are prefixed with a single underscore.

## Types and return values

- Per-SNP vectors are `numpy.ndarray` with shape `(num_snp,)`.
- Matrices are `numpy.ndarray` with shape `(num_snp, num_col)`.
- Missing numeric values use `numpy.nan` (float arrays) or masked arrays where
  the spec requires distinguishing "absent" from "zero".
- Optional arguments default to `None`; implementations test with
  `if arg is None`.
- String arguments (paths, labels) accept `str` or `pathlib.Path`; internally
  convert to `pathlib.Path` at the loader boundary.

## Shard representation

Panels hold an ordered list of shard objects. Genome-wide accessors
(e.g. `ReferencePanel.snp`) concatenate shard arrays lazily or eagerly as
appropriate; callers receive a single numpy array regardless of shard count.

## Error handling

Raise `ValueError` for malformed inputs (missing required columns, length
mismatches, out-of-range indices). Raise `FileNotFoundError` for missing files.
Do not use assertions for user-facing validation; assertions are for internal
invariants only.

Compatibility checks (`is_object_compatible`) log via the standard `logging`
module at `WARNING` level for non-strict mismatches and `ERROR` level for
strict mismatches. They do not raise.

## Optional dependencies

- `numpy` and `scipy.sparse` are required for all numeric objects.
- No other third-party dependencies are required for the core package.
- Optional dependencies (e.g. for cache formats) are declared as extras in
  `pyproject.toml`.
