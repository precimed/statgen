# R implementation conventions

## Compatibility target

The R package is a CRAN-compatible runtime for the language-agnostic `statgen`
object specs. Unless a behavior is explicitly Python-only, MATLAB/Octave-only,
or tied to a runtime-native storage format, R should implement the same
user-facing object contracts as Python and MATLAB/Octave.

Python remains the authoritative implementation for constructing LD panels
from genotype data. R does not reimplement `statgen_build_ld.py`.

The R package must install and load with standard R tooling on Linux, macOS, and
Windows. CRAN checks must not require internet access, large bundled data,
external binaries, or runtime code generation.

## Package layout

The R package should live under `R-package/` rather than the repository root.
This keeps the mixed Python/MATLAB/R repository root from becoming the CRAN
package root and avoids confusion with the package's required `R/` source
directory.

## Naming conventions

- Package, functions, variables, and file names use `snake_case`.
- User-visible S3 classes use the spec object names, for example
  `ReferencePanel`, `AnnotationPanel`, `GenotypePanel`, `Sumstats`, and
  `LDPanel`.
- R uses ordinary exported functions, not MATLAB package syntax. For example,
  R LD distribution preparation is exposed as
  `statgen::prepare_ld_npz_for_r()`, with an optional `extract_npz` argument
  for creating extracted LD shard caches.

When MATLAB and Python naming conventions conflict, R should follow the logical
API names from object specs and express them in idiomatic R `snake_case`. The R
API should not rename existing logical APIs unless the relevant object spec is
changed.

## Dependencies

Required CRAN dependencies should be minimal:

- `Matrix` for sparse matrices;
- `jsonlite` for manifest and metadata JSON;
- `digest` for MD5 checksums;
- `bit64` for exact allele-hash vectors;
- `data.table` for explicit-schema tabular readers;
- `testthat` in `Suggests` for package tests.

The package must not depend on `reticulate`, Python, MATLAB/Octave, PLINK, or
network-backed resources for installation, loading, examples, vignettes, or
CRAN tests. Optional developer-only integrations may live in repository-level
pytest or GitHub Actions workflows outside CRAN checks.

## Exported API

R should export the same logical public API families as the object specs, with
R-native argument syntax and S3 methods where useful.

Panel-level accessors named in object specs may be implemented as S3 methods,
plain functions, or both. R should avoid ambiguous `$`-only public APIs for
important operations because `$` bypasses validation and makes help text less
discoverable.

### S3 generic mapping

Object-spec property and method notation maps to R accessor and operation
generics with the dispatch object as the first argument. For example,
`ReferencePanel.bp` is exposed as `bp(reference_panel)`, and
`ReferencePanel.select_shards(shards)` is exposed as
`select_shards(reference_panel, shards)`.

When a logical method takes another argument in the object specs, R still puts
the dispatch object first. For example, the LD operation specified as
`LDPanel.multiply_r2(M, optional chrX_sex)` is exposed in R as
`multiply_r2(ld_panel, M, chrX_sex = NULL)`.

Compatibility checks intentionally dispatch on the reference panel:
`is_object_compatible(reference, object)` uses the
`is_object_compatible.ReferencePanel` method. Other panel classes do not need
their own `is_object_compatible` methods for the shared reference-alignment
check. To participate in this check, each reference-aligned panel class must
provide a `shards.<Class>` method returning its ordered shard objects, and those
shard objects must expose the reference-alignment metadata required by
[reference.md](reference.md).

The `$` operator is not the public API for important operations or accessors.
Implementations may store ordinary list fields internally, but user-facing code
should use exported accessors such as `bp(x)`, `shards(x)`, and
`reference(ld_panel)`.

## Runtime verbosity

R exports `set_verbosity(level)` and `get_verbosity()` with levels `quiet` and
`info`, matching [SPEC.md](SPEC.md). The current verbosity level may be stored
in a package-private environment. Invalid levels are user-facing validation
errors.

## Shard Representation

R panel objects should be S3 objects backed by ordinary lists. Panels store an
ordered list of shard objects; they must not store an additional panel-wide copy
of concatenated SNP-axis vectors or matrices.

Panel-level accessors should be lazy, non-cached concatenations from shard
payloads. They build the returned vector or matrix on demand using base R or
`Matrix` combiners as appropriate. The panel object must not persist a memoized
full panel-wide SNP-axis payload across accessor calls. Callers that need
repeated use should store a local copy.

`select_shards` should construct a new panel object by subsetting the shard
list, without forcing panel-wide accessor materialization. Under R
copy-on-modify semantics, selected shard payloads should not be deep-copied as
long as shard objects remain immutable after construction.

R S3 classes may define `print.<Class>` methods for interactive summaries.
Print methods should be metadata-only and must not materialize large panel-wide
accessors or lazily decoded reference string payloads.

## Reference string storage

The in-memory/accessor API for `ReferencePanel` and `ReferenceShard` exposes
ordinary R character vectors for `snp`, `a1`, and `a2`.

R reference cache storage should use one newline-delimited character payload per
shard for each of `snp`, `a1`, and `a2`, with no trailing newline. Cache-loaded
references should decode those payloads lazily when the corresponding accessor
is called. This is not required for correctness, but keeps large reference
caches compact, avoids eager materialization of high-cardinality SNP IDs, and
mirrors the MATLAB/Octave reference-cache strategy.

## Type conventions

- Per-SNP vectors are ordinary R vectors with length `num_snp`.
- Matrices have SNPs along rows. Dense inputs to `multiply_r2` may be numeric
  vectors or numeric matrices.
- Missing numeric values use `NaN` for aligned numeric vectors where parity
  with other runtimes is required. Optional absent fields use `NULL`.
- `AnnotationPanel.annomat` is exposed as a `Matrix::dgCMatrix` with
  column names set to `annonames(annotation_panel)`. This is an R-native
  enrichment for ergonomic subsetting such as `annomat(panel)[, "coding"]`;
  annotation order and the authoritative annotation-name accessor remain the
  object contract defined in [annotations.md](annotations.md).
- Paths accept character scalars and are resolved with portable R path
  functions.
- Allele hashes use `bit64::integer64`, not double. This relies on the
  normative `< 2^63` allele-hash bound in [reference.md](reference.md), so
  signed `integer64` storage is exact for the current contract.
- External coordinate metadata such as `start0`, `stop0`, and `index_base`
  remains zero-based exactly as documented. Ordinary R vector indexing in user
  APIs remains one-based.

## Tabular input

Performance-sensitive tabular readers must use explicit column schemas, not
auto-inference. Chromosome labels, SNP IDs, alleles, and FAM identifier columns
must be read as character data; base-pair positions and numeric fields must be
parsed with explicit numeric/integer conversions after column selection. In
particular, readers must preserve chromosome `X` as the string `"X"` and must
not allow mixed numeric/string chromosome columns to be coerced or guessed.

## Error handling

User-facing validation failures should call `stop(..., call. = FALSE)` with a
clear message. Implementations may attach a package-specific condition class
such as `statgen_error`, but they must not depend on non-base condition
packages for ordinary errors.

Use `warning(..., call. = FALSE)` for compatibility warnings. Do not use
`stopifnot` or assertions for user-facing validation.

Compatibility checks such as `is_object_compatible` warn for mismatches and
return scalar logical values. They do not throw.

## CRAN-facing tests

CRAN tests live under `R-package/tests/testthat/`, use only tiny fixtures
included in the built package under `R-package/inst/extdata/`, and cover
ordinary R package behavior: package load, exported API availability, tiny
source loading, non-LD cache round-trips, tiny `.npz` LD loading and
operations, and R LD reference-cache sidecar preparation.

Fixture storage and synchronization follow [testing.md](testing.md).

CRAN tests must not invoke Python, MATLAB/Octave, PLINK, shell scripts, or
internet access. Broader checks belong outside CRAN-facing tests.
