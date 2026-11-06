Tiny CRAN-safe fixtures for R package tests are included here in built package
tarballs.

Portable-format fixtures copied into this directory must be small and derived
from the repository-level canonical fixtures under `tests/fixtures/` when
possible. They are generated working-tree copies, not independent committed
sources. Broader cross-runtime parity tests belong in repository-level pytest
tests.

R-native fixtures with no portable counterpart may be committed directly here.

Current copied fixtures:

- `reference_chr1.bim` is a byte-for-byte copy of
  `tests/fixtures/reference/sharded/1.bim`.
- `traits.tsv.gz` is a byte-for-byte copy of
  `tests/fixtures/sumstats/traits.tsv.gz`; it intentionally exercises
  warning paths for missing optional summary-statistics values.
- `traits_complete.tsv.gz` is a byte-for-byte copy of
  `tests/fixtures/sumstats/traits_complete.tsv.gz`; it is the warning-free
  summary-statistics fixture used by examples and vignettes.
- `anno1.bed`, `anno2.bed`, and `grouped.annot` are byte-for-byte copies of
  tiny annotation fixtures under `tests/fixtures/annotations/`.
- `allele_swaps/` is a byte-for-byte copy of
  `tests/fixtures/allele_swaps/`, used by native R allele-swap tests. This
  directory is ignored by git and must be prepared after a fresh clone.
- `ld/` is an R-ready copy of `tests/fixtures/ld/python/` with an added
  `reference_cache.rds` and manifest `r_reference_cache` fields. It is used by
  examples and vignettes so installed-package examples can call `load_ld()`
  directly without mutating files under `inst/extdata`.

Prepare copied R fixtures before running native package tests or building the
package from a fresh clone. Repository Makefile test/build targets do this
automatically:

```sh
make prepare-r-fixtures
```

Repository pytest checks that copied fixtures remain synchronized with their
canonical source files.
