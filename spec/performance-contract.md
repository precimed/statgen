# Performance and cache contract

## Purpose

Cache files exist to improve performance. The primary objective of caching is
to reduce repeated load/parse/transform cost while preserving object semantics.

Caches are accelerators, not canonical data interchange formats.

LD distribution artifacts are outside this cache contract. They are
runtime-native sparse matrix files by design and are specified in
[ld.md](ld.md). They are not considered cache files, and LD does not expose
public cache save/load APIs.

## Normative requirements

Allele validation and hash generation follow the consolidated requirements in
[allele-handling.md](allele-handling.md), including its bounded-processing
exceptions to the SNP-axis loop restrictions below.

- Cache reads and writes MUST preserve the logical object contract defined in
  object specs.
- Cache validity MUST be tied to compatibility checks (schema version and
  reference checksums where applicable).
- A cache hit MUST be observationally equivalent to loading from canonical
  source inputs, up to documented numeric tolerance. Legacy caches preserve
  the values, mappings, and missingness stored by their writer, subject to
  schema migrations defined in object specs; they need not reproduce additional
  matches available when reloading sources with a newer implementation.
- Cache formats MAY be language-specific and are not required to be portable
  across Python, MATLAB/Octave, and R. Cross-language cache sharing is
  explicitly a non-goal: a cache written by one runtime is not required to be
  loadable by another runtime. Workflows that need to share painted annotation
  data across runtimes should regenerate from canonical BED inputs in each
  runtime.

Cache loading and internal hot paths follow the validation-boundary contract in
[SPEC.md](SPEC.md): cache loaders SHOULD avoid full source-style revalidation
on the default load path, and internal constructors/downstream functions SHOULD
trust already-validated payload values while enforcing only structural
invariants needed for well-formed objects.

For this contract, "native binary storage" means all SNP-axis numeric payloads
are persisted as language-native array variables in binary container formats,
loaded without text parsing of array values.

Cache round-trips MUST preserve array shape, numeric/logical dtype/class, and
sparse-vs-dense representation unless an object spec explicitly documents a
different normalization.

The native-array requirement applies to numeric, logical, and sparse SNP-axis
payloads. Python string-typed cache fields should use Python-native array
storage when present. MATLAB/Octave and R object specs may define
cache-specific encodings for string-typed fields when that improves load
performance while preserving the logical object contract; for example,
MATLAB/Octave reference caches store per-shard string vectors as
newline-delimited character payloads.

For panel-like objects, shard subsetting (`select_shards`) SHOULD avoid deep
copying shard payload arrays by default. Implementations should construct subset
panels from existing shard payloads whenever possible, while preserving object
immutability and observable API semantics.

Loaders SHOULD avoid defensive deep copies of loaded payload arrays or sparse
matrix components when the source container has already materialized normal
runtime-owned arrays and the object contract does not require independent
mutable buffers. Copies SHOULD be reserved for correctness requirements such as
dtype/layout conversion, lifetime management, or explicit mutation isolation.

For performance-sensitive tabular source inputs (for example BIM/TSV-like
files), implementations MUST use language-native tabular readers instead of
line-by-line manual parsing.

- MATLAB/Octave: `textscan` with an explicit schema, as documented in
  [matlab.md](matlab.md), is the required default path.
- Python: dataframe-style readers are the required default path.
- R: language-native tabular readers with explicit column schemas, as documented
  in [R.md](R.md), are the required default path.

Runtime table readers are not required to expose identical permissive parsing
behavior for non-canonical inputs. Implementations should choose the simplest
native-reader path that preserves the canonical source contract and should not
add whole-file prefilters or row-by-row compatibility shims merely to mimic
another runtime's treatment of comments, blank lines, malformed rows, or other
reader-specific recovery behavior. Small probes used to identify a header,
skip leading metadata/comment lines, infer a column count, or choose an
explicit schema are allowed because they support the native-reader path rather
than replacing it.

Documented exceptions are allowed only when a native table reader cannot
correctly represent required input semantics for a specific input shape or
encoding. Any fallback parser MUST be documented in code and tests for that
path.

On the default path, SNP-axis and annotation-axis conversion, validation, and
matrix operations MUST be implemented as column-wise/vectorized operations.
Per-row loops over SNP-axis fields and per-column loops over annotation fields
MUST NOT be used for routine type conversion, validation, annotation painting,
cache serialization, cache loading, or analytical matrix operations on the
default path.

Loops over a small list of source file paths are allowed for file-level
dispatch only; each file's SNP-axis and annotation-axis payload processing must
still use vectorized/native tabular operations. Row-wise or annotation-wise
loops are allowed only under documented exceptions where vectorized operations
cannot correctly express required semantics for a specific input.

Suggested comment format for retained loops in hot paths:

- `PERF: loop retained because ...; vectorization not used because ...`

When using native table loaders for throughput, implementations MAY apply
lighter per-field validation than fully manual parsers, as long as object-level
contract checks and compatibility checks remain enforced.

System temporary directory use must remain small:
- `tempfile` / system temporary-directory use is allowed only for a few small helper objects, approximately 1 MB total per process.
- Large intermediate files must be placed under a user-controlled output or scratch prefix because HPC environments may have small or quota-limited `/tmp` partitions.
- Operations with a natural output path, such as builders and converters, MUST
  place large intermediates under that output path by default, with an explicit
  user-controlled scratch override when appropriate.
- Operations without a natural output path, such as source loaders that need to
  expand compressed inputs before parsing, SHOULD use `STATGEN_SCRATCH` when it
  is set. If `STATGEN_SCRATCH` is unset, they may fall back to a source-adjacent
  scratch directory when writable, but they MUST NOT silently place large
  intermediates in the system temporary directory.

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

## R cache format requirements

R caches MUST use R-native binary storage suitable for direct numeric array
loading. The default R cache format is RDS unless an object-specific R spec
defines another R-native format.

For SNP-axis vectors and matrices, caches MUST store data as native R numeric,
logical, integer64, or sparse matrix objects. Implementations MUST NOT encode
these arrays as JSON blobs, string payloads, or other text-serialized
representations inside cache files.

Small metadata fields (for example schema/version strings, shard labels, and
checksums) MAY use R lists, data frames, and character vectors, but numeric
array payloads remain native binary objects.

Metadata-only caches MAY use text-oriented serialization when they contain no
SNP-axis numeric/logical/sparse vectors or matrices.

## MATLAB/Octave cache format requirements

MATLAB/Octave cache files MUST use MATLAB-native `.mat` storage.
MATLAB/Octave cache writers accept an optional `format` setting using the same
language-specific syntax as other `.mat` producers. Supported values are `v7`,
`v7.3`, and `v5`. The default is `v7.3` under MATLAB and `v7` under Octave.
An explicit `format` overrides the runtime default. Cache writers must raise an
error if the selected format cannot store a variable because of its size,
rather than report success with an incomplete cache. The caller's warning
settings must be restored after saving, including on failure.

For SNP-axis vectors and matrices, caches MUST store numeric, logical, and
sparse data as native arrays in `.mat` variables. Implementations MUST NOT
encode numeric, logical, or sparse arrays as JSON blobs, string payloads, or
other text-serialized representations inside cache files. This prohibition does
not apply to string-typed SNP-axis fields; those fields follow the
object-specific MATLAB/Octave cache encoding, such as the per-shard
newline-delimited character payloads defined for reference caches.

Small metadata fields (for example schema/version strings, shard labels, and
checksums) MAY use native MATLAB structs/cells/chars/strings, but numeric array
payloads remain native arrays.

Metadata-only caches MAY use text-oriented serialization when they contain no
SNP-axis numeric/logical/sparse vectors or matrices.
