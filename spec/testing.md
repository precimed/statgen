# Testing strategy

`statgen/tests` should use `pytest`. MATLAB/Octave validation should be done
with Octave; installed Octave is a prerequisite for tests that exercise MATLAB
functionality.

## Fixture data

Tests should generate tiny portable object fixtures with:

- 2 chromosomes, including an autosome and chrX where possible;
- a few SNPs per chromosome;
- one PLINK bfile genotype shard per chromosome;
- one LD shard per chromosome;
- one annotation panel with a few annotations;
- one `Sumstats` object with finite and missing entries.
- both sharded and one-shard/non-sharded reference layouts where supported.

## Cross-language consistency

Python and Octave loaders should read the same portable source files and report
matching:

- dimensions;
- vector values;
- reference shard/panel ordering;
- sparse LD triplets and materialized sparse matrices;
- genotype metadata; genotype slices only after the deferred genotype access
  contract is specified;
- annotation masks and LD-weighted annotation arrays within numeric tolerance.

Cache conversion tests should verify that native cache outputs match the
portable source files, not that caches match each other directly.

## Edge cases and failures

Tests should also cover:

- reference checksum computation and `ReferencePanel.is_object_compatible`
  mismatch reports;
- cache invalidation when a supplied reference has different shard checksums;
- little-endian LD triplet loading and binary array length validation;
- chrX LD shard groups, including sex-specific default selection;
- sharded versus one-shard/non-sharded loading behavior;
- malformed source files with missing required columns;
- summary-statistic `p` handling for missing values, `0`, negative values, and
  values greater than `1`;
- exact BED/reference chromosome matching and interval boundary behavior;
- mismatch fixtures proving that `statgen` does not normalize aliases such as
  `chr1` to `1` or `23` to `X`;
- genotype tests only after the deferred genotype access contract is specified.

## Success criteria

The library is complete when univariate LDSC-style preparation and bivariate
genetic-correlation preparation can be expressed using only public loaders,
accessors, `LDPanel.multiply_r2`, `fast_prune`, and plain arrays. Tests should
exercise at least loading reference, LD, annotations, and one or two sumstats
objects; computing an LD-weighted annotation matrix; applying one random
`fast_prune` mask; and deriving bivariate `z1 * z2` arrays. Python and
MATLAB/Octave should produce matching fixture outputs for those workflows.

## Test boundaries

Tests should cover `statgen` behavior, not upstream harmonization. Input test
fixtures should already satisfy the reference-genome and allele contract.
