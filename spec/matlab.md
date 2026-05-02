# MATLAB/Octave implementation conventions

## Compatibility target

The MATLAB implementation is the primary runtime contract. Octave is used as a
unit-test harness and compatibility proxy, not as the normative behavior
source. When MATLAB and Octave differ, implementation behavior is defined by
MATLAB semantics unless explicitly documented otherwise.

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
accessors return a single matrix or vector in panel order.

For panel-like objects (`ReferencePanel`, `Sumstats`, `AnnotationPanel`),
genome-wide accessors should be implemented as **lazy, non-cached**
concatenations from shard payloads:

- constructors should not eagerly materialize/store full panel-wide SNP-axis
  vectors or matrices in addition to per-shard payloads;
- accessors should build the concatenated output on demand from `shards`;
- "non-cached" means the panel object should not persist a memoized full
  panel-wide SNP-axis payload (for example a hidden cached `annomat`/`zvec`)
  across accessor calls; each accessor call may recompute concatenation;
- callers that need repeated use should store one local copy, e.g.
  `A = panel.annomat`.

This keeps MATLAB/Octave behavior aligned with the memory model in the object
specs and with Python panel-accessor semantics.

## Error handling

Use `error('statgen:errorID', 'message ...')` with a namespaced error ID for
all user-facing validation failures. Do not use `assert` for input validation.

Compatibility checks log via `warning('statgen:compat', ...)` for mismatches.
They return a scalar logical and do not throw.

## Runtime compatibility guardrails

- For BIM/TSV-like source inputs, native table readers must use explicit text
  schema for mixed-label columns (for example `chr`) rather than auto-inference.
  This avoids label coercion (for example `X` becoming `NaN`) and preserves
  input-contract semantics.
- Internal struct metadata keys must be valid MATLAB identifiers (for example
  `statgen_var_names__`), not names that rely on permissive dynamic-field
  behavior (for example leading-underscore keys such as `_var_names`), because
  this can fail at runtime under MATLAB even when other environments appear to
  tolerate it.
- Use `fprintf` for formatted output in shared runtime/test snippets intended
  to execute under MATLAB. Octave-only output helpers (for example `printf`)
  must be wrapped or normalized at the harness boundary.
- Cross-runtime differences must be isolated in boundary helper functions
  (I/O parsing, hashing, output formatting, engine/session integration), not
  duplicated across core loader/object logic.
- Test execution may use one persistent MATLAB engine session per pytest run.
  MATLAB-facing code and test harness paths must not assume process-per-call
  isolation.
