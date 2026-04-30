# MATLAB/Octave implementation conventions

## Compatibility target

All MATLAB/Octave code must run under Octave (current stable release). MATLAB
compatibility is desirable but secondary; when the two differ, Octave behavior
takes precedence.

## Package layout

Internal implementation uses the `+statgen` package directory, which Octave
resolves as the `statgen` namespace. Flat `statgen_*` wrapper scripts are added
only for the public loader entry points and top-level primitives — not for
internal classes or helpers.

```
matlab/
  +statgen/
    load_reference.m
    load_ld.m
    load_sumstats.m
    load_annotations.m
    fast_prune.m
    ReferenceShard.m
    ReferencePanel.m
    LDShard.m
    LDPanel.m
    AnnotationShard.m
    AnnotationPanel.m
    SumstatsShard.m
    Sumstats.m
    private/           % internal helpers not visible outside the package
  statgen_load_reference.m
  statgen_load_ld.m
  statgen_load_sumstats.m
  statgen_load_annotations.m
  statgen_fast_prune.m
```

The `matlab/` directory itself must be on the Octave/MATLAB path. The
`+statgen` subdirectory is discovered automatically from there; it must not
be added to the path directly.

## Wrapper pattern

Each flat entry point is a thin forwarding wrapper using the standard
`varargout`/`varargin` pattern:

```matlab
function varargout = statgen_load_reference(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_reference(varargin{:});
end
```

The `max(nargout, 1)` allocation avoids an empty index expression when the
caller discards all outputs. The internal `statgen.*` function is responsible
for all argument validation and logic; wrappers contain no other code.

Wrappers exist only for: `statgen_load_reference`, `statgen_load_ld`,
`statgen_load_sumstats`, `statgen_load_annotations`, `statgen_fast_prune`.
No wrappers are created for classes or internal helpers.

## Naming conventions

- Functions and variables: `snake_case`.
- Classes: `PascalCase` matching the spec object names.
- Private package helpers: placed under `+statgen/private/`; not accessible
  outside the package.

## Types and return values

- Per-SNP vectors are column vectors (shape `n × 1`).
- Matrices have SNPs along rows (shape `n × k`).
- Missing numeric values use `NaN`.
- Optional arguments use `[]` as the absent sentinel; implementations test
  with `isempty(arg)`.
- String arguments are character arrays or string scalars; convert to
  `char` at the loader boundary for compatibility with Octave string handling.

Cache format requirements are defined normatively in
[performance-contract.md](performance-contract.md). In particular, MATLAB/Octave
caches must be `.mat` files with SNP-axis vectors/matrices stored as native
numeric/logical/sparse arrays.

## Shard representation

Panels hold an ordered cell array or struct array of shard objects. Genome-wide
accessors concatenate shard arrays and return a single matrix or vector.

## Error handling

Use `error('statgen:errorID', 'message ...')` with a namespaced error ID for
all user-facing validation failures. Do not use `assert` for input validation.

Compatibility checks log via `warning('statgen:compat', ...)` for mismatches.
They return a scalar logical and do not throw.
