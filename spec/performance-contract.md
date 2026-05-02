# Performance and cache contract

## Purpose

Cache files exist to improve performance. The primary objective of caching is
to reduce repeated load/parse/transform cost while preserving object semantics.

Caches are accelerators, not canonical data interchange formats.

## Normative requirements

- Cache reads and writes MUST preserve the logical object contract defined in
  object specs.
- Cache validity MUST be tied to compatibility checks (schema version and
  reference checksums where applicable).
- A cache hit MUST be observationally equivalent to loading from canonical
  source inputs, up to documented numeric tolerance.
- Cache formats MAY be language-specific and are not required to be portable
  across Python and MATLAB/Octave. Cross-language cache sharing is explicitly a
  non-goal: a cache written by one runtime is not required to be loadable by the
  other. Workflows that need to share painted annotation data across runtimes
  should regenerate from canonical BED inputs in each runtime.

Cache loaders SHOULD assume cache payloads are already valid from cache-build
time and SHOULD avoid full source-style revalidation on the default load path.
At minimum, cache loaders MUST enforce lightweight compatibility gates (for
example schema/version checks and reference checksum checks where applicable).

For this contract, "native binary storage" means all SNP-axis numeric payloads
are persisted as language-native array variables in binary container formats,
loaded without text parsing of array values.

Cache round-trips MUST preserve array shape, numeric/logical dtype/class, and
sparse-vs-dense representation unless an object spec explicitly documents a
different normalization.

For panel-like objects, shard subsetting (`select_shards`) SHOULD avoid deep
copying shard payload arrays by default. Implementations should construct subset
panels from existing shard payloads whenever possible, while preserving object
immutability and observable API semantics.

For performance-sensitive tabular source inputs (for example BIM/TSV-like
files), implementations MUST use language-native tabular readers instead of
line-by-line manual parsing.

- MATLAB/Octave: `readtable` (or a semantically equivalent native table reader)
  is the required default path.
- Python: dataframe-style readers are the required default path.

Documented exceptions are allowed only when a native table reader cannot
correctly represent required input semantics for a specific input shape or
encoding. Any fallback parser MUST be documented in code and tests for that
path.

On the default path, SNP-axis conversion and validation MUST be implemented as
column-wise/vectorized operations. Per-row loops over SNP-axis fields MUST NOT
be used for routine type conversion or validation on the default path.

Row-wise loops are allowed only under documented exceptions where vectorized
operations cannot correctly express required semantics for a specific input.

Suggested comment format for retained loops in hot paths:

- `PERF: loop retained because ...; vectorization not used because ...`

When using native table loaders for throughput, implementations MAY apply
lighter per-field validation than fully manual parsers, as long as object-level
contract checks and compatibility checks remain enforced.

## Python cache format requirements

Python caches MUST use Python-native binary storage suitable for direct numeric
array loading.

For SNP-axis vectors and matrices, caches MUST store data in native numeric,
logical, or sparse array form (for example NumPy/SciPy-native array payloads).
Implementations MUST NOT encode these arrays as JSON blobs, string payloads, or
other text-serialized representations inside cache files.

Small metadata fields (for example schema/version strings, shard labels, and
checksums) MAY use plain JSON-compatible scalars/maps, but numeric array
payloads remain native binary arrays.

Metadata-only caches MAY use text-oriented serialization when they contain no
SNP-axis numeric/logical/sparse vectors or matrices.

## MATLAB/Octave cache format requirements

MATLAB/Octave cache files MUST use MATLAB-native `.mat` storage.

For SNP-axis vectors and matrices, caches MUST store data as native numeric,
logical, or sparse arrays in `.mat` variables. Implementations MUST NOT encode
these arrays as JSON blobs, string payloads, or other text-serialized
representations inside cache files.

Small metadata fields (for example schema/version strings, shard labels, and
checksums) MAY use native MATLAB structs/cells/chars/strings, but numeric array
payloads remain native arrays.

Metadata-only caches MAY use text-oriented serialization when they contain no
SNP-axis numeric/logical/sparse vectors or matrices.
