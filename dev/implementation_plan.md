# statgen phased implementation plan

This plan turns the object-centered specification in `spec/` into a
working Python and MATLAB/Octave package. The guiding constraints are:

- loaders create immutable, load-only objects from portable disk formats;
- all per-SNP vectors and matrix rows align to `ReferencePanel` order;
- portable storage stays language-agnostic;
- MATLAB/Octave and Python behavior is validated from the same fixtures;
- genotype matrix access remains deferred until its contract is specified.

## Phase 0: package skeleton and shared test fixtures

Goal: establish the repository layout, dependency assumptions, and fixture data
used by every later phase.

Implementation tasks:

- create `python/pyproject.toml` and `python/statgen/` with empty module stubs
  for `reference`, `sumstats`, `annotations`, `ld`, `genotype`, and `_utils`;
  package is installable via `pip install -e python/` from the repo root;
- create `matlab/+statgen/` package with a minimal `version.m` function as the
  first callable entry point, plus the corresponding `statgen_version.m` flat
  wrapper using the `max(nargout, 1)` varargout pattern; no other stubs needed
  yet;
- create `script/` with placeholder scripts for LD conversion/build tooling;
- create `tests/conftest.py` with:
  - a session-scoped `octave_available` fixture that probes `shutil.which("octave")`;
  - a `pytest.mark.octave` marker (registered in `pyproject.toml`);
  - a `skipif_no_octave` skip decorator;
- create `tests/fixtures/generate.py`: standalone script that writes all
  committed fixture files deterministically using explicit little-endian dtypes;
  run with `python tests/fixtures/generate.py` to regenerate;
- generate and commit baseline fixtures under `tests/fixtures/` (this is the
  authoritative concrete baseline fixture composition):
  - one autosome and chrX;
  - sharded `.bim` files and one non-sharded `.bim`;
  - inputs using genomatch NCBI contig naming;
  - small BED annotation files with overlapping, adjacent, and boundary
    intervals;
  - one `.tsv.gz` sumstats file with finite, missing, `p == 0`, invalid, and
    absent-reference variants;
  - minimal LD triplet directories, including a chrX shard group with
    `combined`, `male`, and `female` shards;
  - tiny PLINK bfile metadata fixtures for genotype metadata loading.

Tests and acceptance criteria:

- pytest runs without Octave and skips `@pytest.mark.octave` tests cleanly;
- at least one `@pytest.mark.octave` smoke test: invoke `statgen.version()` via
  subprocess, verify the returned value matches the expected string;
- fixture-generation script is deterministic and its output is committed;
- `Makefile` with three thin targets: `install`, `fixtures`, `test`;
- CI/local test commands are documented in `README.md`.

## Phase 1: shared conventions and reference panels

Goal: implement the base coordinate system that every other object depends on,
while preserving the post-genomatch variant contract exactly.

Implementation tasks:

- implement BIM parsing and validation for required columns:
  `chr`, `snp`, `cm`, `bp`, `a1`, `a2`;
- implement reference checksum computation from exact loaded `chr:bp:a1:a2\n`
  rows with lowercase MD5 output;
- implement Python `ReferenceShard` and `ReferencePanel` classes with
  read-only accessors, shard offsets, and index conversion helpers;
- implement MATLAB/Octave reference loading, checksum, and panel accessors;
- implement `load_reference(path, no_shard=False)` for:
  - sharded path templates containing `@`;
  - single BIM split by chromosome;
  - single BIM as one `"all"` shard when `no_shard=True`;
- implement reference cache save/load in both languages;
- implement `ReferencePanel.is_object_compatible(object, strict=False)` as a
  boolean compatibility check that may log warnings/errors but does not raise
  on mismatches.

Tests and acceptance criteria:

- Python and Octave return the same shard labels, row order, offsets, and
  checksums for all fixture reference layouts;
- single-file default loading splits by chromosome, while `no_shard=True`
  preserves one `"all"` shard;
- compatibility checks return `false` without raising for row-count or checksum
  mismatches;
- chrX is supported without chromosome 1-22 hard-coding.

## Fixture sourcing policy for Phase 2+

Follow `spec/testing.md` as the source of truth for fixture growth policy.
From Phase 2 onward, tests should use:

- Phase 2 (`Sumstats`): checked-in canonical TSV fixtures for cross-language
  alignment and `logpvec` contract; on-the-fly malformed/missing-column files
  for failure-path tests.
- Phase 3 (`Annotations`): checked-in canonical BED fixtures for overlap,
  adjacency, boundary, and exact chromosome-label behavior; on-the-fly
  permutation/stress boundary inputs.
- Phase 4 (`LD`): checked-in canonical LD shards/shard-group fixtures for
  parity and compatibility; on-the-fly binary/metadata corruption cases.
- Phase 5 (LD build): on-the-fly synthetic PLINK-like tables for unit tests
  and on-the-fly generated outputs for integration validation.
- Phase 6 (`Genotype` metadata): checked-in canonical bfile metadata fixtures;
  on-the-fly inconsistent/missing shard-file cases.
- Phase 7 (workflow validation): default to checked-in canonical workflow
  inputs; only add new checked-in fixtures when a new workflow contract cannot
  be represented by temporary derived inputs.

## Phase 2: summary statistics alignment

Goal: load one portable GWAS TSV per trait/source and project it into reference
order.

Implementation tasks:

- implement required-column validation for `chr`, `bp`, `a1`, `a2`, `z`, `n`;
- implement optional field handling for `p`, `beta`, `se`, `eaf`, and `info`;
- join rows to the reference by exact `chr:bp:a1:a2`;
- split aligned vectors into `SumstatsShard`s matching the reference shards;
- represent absent variants and missing numeric values as `NaN`;
- derive `logpvec` from `p` only:
  - missing `p` -> `NaN`;
  - `p == 0` -> `Inf`;
  - `p < 0` or `p > 1` -> `NaN`;
- implement Python and MATLAB/Octave accessors for genome-wide concatenated
  vectors;
- implement aligned cache save/load tied to reference checksums.

Tests and acceptance criteria:

- Python and Octave aligned vectors match the fixture values and missingness;
- absent reference variants remain rows with `NaN`, not dropped rows;
- optional fields have consistent missing-field sentinels;
- malformed TSV files with missing required columns fail with clear errors;
- cache loading detects incompatible references.

## Phase 3: annotation painting

Goal: paint BED interval annotations onto a reference panel and expose binary
annotation matrices.

Implementation tasks:

- implement BED loading using only columns 1-3;
- sort and merge overlapping or adjacent intervals per chromosome;
- convert BIM 1-based `bp` to 0-based positions for interval membership;
- implement binary-search interval painting in reference row order;
- derive annotation names from BED filenames;
- implement Python `AnnotationShard` and `AnnotationPanel` classes with
  `annomat` and `annonames` accessors;
- implement MATLAB/Octave annotation painting with matching semantics;
- implement annotation cache save/load tied to reference checksums.

Tests and acceptance criteria:

- boundary tests cover `[start, end)` membership exactly;
- adjacent BED intervals merge and produce the same mask as separate intervals;
- chromosomes absent from a BED file produce all-zero masks;
- Python and Octave `annomat` and `annonames` match for sharded and one-shard
  reference layouts;
- complement masks and LD-weighted matrices are left to user-space arrays.

## Phase 4: LD portable loader and numerical primitives

Goal: load portable LD triplet directories, validate metadata, and expose the
two public LD operations.

Implementation tasks:

- implement metadata parsing for plain `ld_shard` directories and chrX
  `ld_shard_group` directories;
- load little-endian `ld_idx1.i32`, `ld_idx2.i32`, `ld_r.f32`, and
  `mafvec.f32`, with explicit byte-order handling;
- implement `load_ld(ld_dir, reference, chrX_default_sex=None)` discovery from
  either a single shard directory or a panel root; for panel roots, derive
  expected shard directories from `reference` shard labels as
  `ld_dir/<shard_label>/`, using `ld_dir/all/` for one-shard references;
- validate `num_snp`, `num_ld`, triplet lengths, `idx1 < idx2`, bounds, value
  arrays, diagonal policy, and reference checksums;
- implement sparse LD materialization in Python and MATLAB/Octave;
- implement `LDPanel.mafvec(sex=None)` with chrX default-sex selection;
- implement `LDPanel.multiply_r2(M, sex=None)` without building a dense
  genome-wide LD matrix;
- implement `fast_prune(logpvec, ld_panel, r2_threshold=0.2, sex=None)` with
  stable significance ordering and per-shard pruning;
- implement LD cache save/load in both languages, preserving chrX sex shards
  and following the cache path-template convention: `@` is replaced by the
  reference shard label, and `#` is replaced by the chrX sex label when
  sex-specific cache files are used; autosomal cache paths omit the sex label.

Tests and acceptance criteria:

- Python and Octave sparse matrices, `mafvec`, `multiply_r2`, and `fast_prune`
  outputs match on fixtures within numeric tolerance;
- chrX default sex is restored from cache and can be overridden by callers;
- invalid metadata or binary array lengths fail clearly;
- `multiply_r2` accepts both vectors and matrices and returns the same shape;
- non-sharded `"all"` references and LD panels are supported.

## Phase 5: LD build and conversion script

Goal: create portable LD panels from PLINK bfiles using reproducible command
metadata.

Implementation tasks:

- implement `script/statgen_build_ld.py`;
- support sharded bfile input with `@` and non-sharded input;
- for non-sharded input, split by chromosome by default and support
  `--no-shard` for one `"all"` LD output;
- run PLINK2 `--freq` and `--r` with `--keep-allele-order`;
- default to a 10,000 kb LD window and `r2 >= 0.05` storage threshold, with
  command-line overrides recorded in metadata;
- convert PLINK output into portable triplets and `mafvec`;
- require unambiguous SNP row mapping for `--no-shard`;
- for chrX, build `combined`, `male`, and `female` shard-group outputs by
  default using FAM sex codes, with `--no-sex-split` for combined-only output;
- record build command, PLINK version if available, sample count, window,
  threshold, and reference checksum in metadata.

Tests and acceptance criteria:

- unit-test conversion from synthetic PLINK-like LD/frequency tables to portable
  triplets without requiring PLINK2;
- integration tests that require PLINK2 are optional and skipped when PLINK2 is
  unavailable;
- generated portable LD loads successfully through Phase 4 loaders;
- chrX sex splitting rejects missing FAM sex values unless `--no-sex-split` is
  used.

## Phase 6: genotype metadata objects

Goal: implement the deferred, safe subset of genotype support: metadata loading
without dense genotype access.

Implementation tasks:

- implement `GenotypeShard` and `GenotypePanel` loaders for PLINK bfile
  prefixes;
- support single-prefix input and ordered lists/templates of sharded prefixes;
- validate `.bed`, `.bim`, and `.fam` presence;
- load BIM and FAM metadata, retaining file paths and sample order;
- for sharded genotype panels, require identical FAM rows in the same order
  across all shards, reflecting standard per-chromosome PLINK bfile practice;
- expose read-only `num_snp`, `num_sample`, `bim`, and `fam` accessors;
- do not expose genotype matrix readers until slicing, missingness, and
  reference-alignment behavior are specified.

Tests and acceptance criteria:

- Python and Octave genotype metadata accessors match fixture files;
- sharded panels reject inconsistent FAM rows or sample order across shards;
- tests do not depend on dense genotype decoding.

## Phase 7: workflow-level validation

Goal: prove that the implemented object API can support intended downstream
statistical-genetics workflows without internal-field access.

Implementation tasks:

- implement small Python and Octave examples for:
  - univariate partitioned LDSC-style array preparation;
  - bivariate/cross-trait array preparation;
  - random-prune weight estimation;
- keep regression and model fitting outside `statgen`, using plain arrays from
  accessors;
- document common usage in `README.md`.

Tests and acceptance criteria:

- examples use only public loaders, accessors, `LDPanel.multiply_r2`, and
  `fast_prune`;
- Python and Octave fixture workflows produce matching arrays within tolerance;
- no test reaches into private shard internals to reimplement loading or LD
  logic.

Workflow test sketches:

```text
ref    = load_reference(ref_path)
ld     = load_ld(ld_path, ref)
annot  = load_annotations(bed_paths, ref)
ss     = load_sumstats(sumstats_path, ref)

mafvec = ld.mafvec
sig2_i = 2 * mafvec * (1 - mafvec)
annomat_ld = ld.multiply_r2(annot.annomat * sig2_i[:, None])

mask = isfinite(ss.zvec) & isfinite(ss.nvec) & (mafvec > maf_threshold)
randvec = uniform(0, 1, ref.num_snp)
randvec[~mask] = NaN
mask = mask & isfinite(fast_prune(randvec, ld, r2_threshold=0.1))
weights = 1.0 / ld.multiply_r2(mask.astype(float))[mask]
```

```text
ss1 = load_sumstats(sumstats_path_1, ref)
ss2 = load_sumstats(sumstats_path_2, ref)
mask = isfinite(ss1.zvec) & isfinite(ss2.zvec) & isfinite(ss1.nvec) & isfinite(ss2.nvec)
zprod = ss1.zvec[mask] * ss2.zvec[mask]
```

Random-prune weight estimation should be documented and tested as repeated
calls to `fast_prune` with excluded SNPs pre-set to `NaN`:

```text
weights = zeros(ref.num_snp)
for i in 1..n_iter:
    randvec = uniform(0, 1, ref.num_snp)
    randvec[~mask] = NaN
    weights += isfinite(fast_prune(randvec, ld, r2_threshold))
weights /= n_iter
```

## Phase 8: hardening, packaging, and documentation

Goal: make the package usable by project pipelines and stable enough for
external analysis scripts.

Implementation tasks:

- document installation and optional dependencies:
  Python numeric stack, Octave, PLINK2, and any cache libraries;
- finalize public API names and import paths in Python;
- finalize MATLAB/Octave path setup and class/function naming;
- add concise API reference pages or docstrings for each object family;
- add command-level documentation for `statgen_build_ld.py`;
- add schema-version checks for portable metadata and caches;
- define backwards-compatibility policy for schema `0.1`.

Tests and acceptance criteria:

- full pytest suite passes locally;
- Octave tests pass when Octave is installed and skip cleanly otherwise;
- documentation examples run against fixture data;
- project-specific pipelines can load independent object paths without a
  bundle-root manifest.

## Dependency ordering

The critical path is:

1. shared conventions and references;
2. objects that align to references without sparse LD: sumstats and annotations;
3. LD loading and numerical primitives;
4. LD building;
5. genotype metadata;
6. workflow examples and hardening.

Genotype dense matrix access should remain out of scope until the genotype spec
is expanded. Collection manifests should also remain optional until object
loaders and caches are stable.

## Language implementation sequencing

Python and MATLAB/Octave should implement the same public object contracts, but
they should not be developed in lockstep while low-level details are still
moving. Python should lead each implementation phase, and MATLAB/Octave should
catch up at stable phase gates using the same portable fixtures and expected
outputs.

Recommended sequencing:

1. Implement Phase 0 and Phase 1 in Python first, including fixture generation,
   reference checksums, reference loading, and compatibility checks.
2. Record fixture expectations in language-neutral files or pytest helpers that
   compare plain arrays and metadata, not Python object internals.
3. Implement MATLAB/Octave Phase 1 against those same fixtures before advancing
   too far into higher-level object behavior.
4. Keep Python roughly one phase ahead for sumstats and annotations, then add
   MATLAB/Octave parity before treating each phase as complete.
5. Defer MATLAB/Octave LD implementation until the Python LD portable loader,
   triplet validation, sparse matrix construction, and cache layout are stable.
   LD has the highest risk of format churn and should not be implemented twice
   while those details are changing.
6. Use cross-language tests as phase gates: a phase is not complete until
   Python and MATLAB/Octave produce matching values from the same portable
   fixtures, except for explicitly deferred genotype matrix access.

This is Python-first, not Python-only. MATLAB/Octave parity remains part of the
acceptance criteria, but it should follow stable contracts instead of driving
early implementation choices before the storage and cache details are proven.

## Open implementation decisions

- Python cache format for reference, sumstats, and annotations.
- MATLAB/Octave cache field layout: it should optimize simple load/save and
  checksum validation, not mirror Python internals.
- Whether the Python package should be installable as part of the repository or
  via a nested `pyproject.toml` under `python`.
- How much of existing `ld/*` code should be mined for tests or conversion
  details after the clean LD implementation is underway.
