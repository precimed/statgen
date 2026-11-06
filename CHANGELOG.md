# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Python, MATLAB/Octave, and R source loaders automatically match reversed
  SNVs, including A/T and C/G under the same-strand assumption. Both loaders
  accept `allow_indel_swaps` to also match reversed indels with exactly one
  single-base allele; pairs with two multi-base alleles remain exact-only.
- Genotype panels expose `is_allele_swapped` on the full reference SNP axis.
  Genotype fetching via `fetch_genotypes_int8` and `fetch_genotypes` orients
  hardcalls to the reference panel's `a1` (ALT) allele before ploidy scaling.

### Changed
- When source alleles are reversed, Z scores and beta change sign and EAF
  becomes `1 - EAF`, so they describe the reference panel's `a1` (ALT) allele.
  Other statistics and missing values are preserved.
- `load_reference` now rejects repeated variants and SNPs listed with both
  allele orders, such as A/G and G/A at the same position.
- Genotype caches now use `genotype_cache/0.2` to record allele swaps.

### Fixed
- R variant matching now sorts allele hashes at or above
  `0x7FF0000000000000` correctly. Previously, two such alleles at one position
  could lose exact matches, and duplicate keys could load silently.

## [0.3.4] - 2026-05-27

### Added
- Added grouped binary annotation loading in Python, MATLAB/Octave, and R.
  `load_annotation(...)` now accepts `group_column` to create one binary
  annotation per distinct group value in a BED-like source column.

## [0.3.3] - 2026-05-26

### Added
- Added continuous numeric annotation support in Python, MATLAB/Octave, and R.
  New `load_annotation(...)` APIs load one BED-like source file, support
  selecting one or more value columns, optional headers, annotation name
  overrides, and metadata supplied directly or through sidecar files.
- `AnnotationPanel` now exposes per-column binary status and metadata through
  `is_binary` / `annotation_metadata` accessors or properties. Binary
  `load_annotations(...)` batch loaders can also attach per-file metadata.

### Changed
- Annotation matrices are now sparse numeric matrices instead of binary-only
  masks.
- Annotation caches now use schema `annotations_cache/0.2` to store
  `is_binary` and `annotation_metadata`. Existing `annotations_cache/0.1`
  caches remain loadable and are upgraded with binary status set to true and
  empty metadata.

### Fixed
- Updated R package metadata to satisfy CRAN-style quoting of software names in
  `DESCRIPTION`.

## [0.3.2] - 2026-05-19

### Changed
- Refined R examples, vignette, tutorial install notes, and package metadata for
  clearer installed-package workflows.

## [0.3.1] - 2026-05-16

### Changed
- The R package now imports `data.table` and `R.utils` for faster source-file
  parsing, including `.gz` summary-statistics inputs.

### Fixed
- Improved `load_sumstats(...)` performance in Python, MATLAB/Octave, and R,
  especially for large summary-statistics files with mostly single-base
  alleles.

## [0.3.0] - 2026-05-14

### Added
- Added the R package runtime with reference, summary statistics, annotation,
  genotype, and LD support. Release automation now checks the R package and
  attaches a CRAN-style source tarball.

### Changed
- Python and MATLAB/Octave `create_sumstats(...)` now use argument names
  `p`, `z`, `n`, `beta`, `se`, `eaf`, and `info`; `create_annotations(...)`
  now uses `annotation_matrix` and `annotation_names`.
- Python cache save helpers and panel `.save_cache(...)` methods no longer
  accept the unused `format` argument.

### Fixed
- MATLAB/Octave summary-statistics loading now parses numeric TSV columns with
  native numeric `textscan` formats instead of reading all fields as strings
  first.

 ## [0.2.6] - 2026-05-11

### Added
- Python and MATLAB/Octave now provide `load_ld_reference(...)` to load the
  manifest-declared `ReferencePanel` from an LD distribution root without
  loading LD matrices or manually parsing `ld_manifest.json`.
- Python and MATLAB/Octave now provide global runtime verbosity controls with
  `set_verbosity(...)` and `get_verbosity(...)`. `load_ld` reports per-shard
  loading progress at the default `info` level and suppresses it at `quiet`.
- MATLAB/Octave panel and shard objects now have concise display summaries for
  interactive use.

### Changed
- LD distributions now use the manifest-declared reference cache for
  `load_ld(...)`; runtime loaders no longer rely on hard-coded reference cache
  filenames. Manifest validation also checks the manifest-declared reference
  cache against bundled reference BIM shards.
- MATLAB/Octave `create_ld_mat_manifest(...)` now writes the manifest-declared
  reference cache into the LD distribution.

### Fixed
- MATLAB/Octave LD manifest creation and validation now stream the internal MD5
  fallback instead of reading entire large `.mat` LD shards into Java heap
  memory. A warning is emitted when the system MD5 command is unavailable and
  MATLAB-side hashing is used.

 ## [0.2.5] - 2026-05-10

### Fixed
- MATLAB reference checksums now use the same unpadded base-pair formatting and
  MD5 implementation path as Python/Octave, avoiding checksum mismatches for
  LD reference shards.

 ## [0.2.4] - 2026-05-09

### Added
- `GenotypePanel.fetch_genotypes(...)` now support
  `haploid_mode="ploidy_scaled"` to scale decoded allele counts by per-variant
  male/female ploidy.
- `ReferencePanel` and `ReferenceShard` expose `is_single_nucleotide_variant` and
  `is_strand_ambiguous` masks.
- Source-to-reference loading for sumstats and genotype data now warns when an
  unmatched source-side variant would match the reference if `a1` and `a2` were
  swapped. The variant remains unmatched; loaders still never flip alleles.

### Changed
- Updated documentation and tutorials for genotype loading, genotype fetching,
  LD distribution workflows, and MATLAB/Octave API usage.
- Missing chrX genotype `.ploidy` sidecars now default chrX rows to male
  haploid and female diploid and warn.

### Fixed
- Reference loading now fails clearly when a reference variant has identical
  `a1` and `a2` alleles.

## [0.2.3] - 2026-05-08

### Added
- Python and MATLAB/Octave can load reference-aligned PLINK 1 genotype datasets
  with `load_genotype(bfile_prefix, reference)`, including non-sharded and
  `@`-sharded bfile layouts.
- `GenotypePanel.fetch_genotypes(...)` fetches selected genotype calls without
  loading the full genotype dataset.
- Genotype metadata caches can be saved and reloaded with
  `save_genotype_cache(...)` and `load_genotype_cache(...)` to speed up repeated
  genotype loading.
- chrX genotype loading supports `.ploidy` sidecars and FAM files with chrX-only
  sample subsets.
- `GenotypePanel.is_present` and `Sumstats.is_present` provide a boolean mask in
  reference coordinates for variants present after aligning input data to the
  reference.

### Changed
- Summary statistics now require `p` values and allow `z` and `n` to be omitted

## [0.2.2] - 2026-05-06

### Added
- Sumstats loaders now accept genomatch-style summary-statistic headers in
  Python and MATLAB/Octave.
- MATLAB/Octave flat cache layouts now expose user-inspectable top-level
  variables for reference, sumstats, and annotation caches.
- Reference objects now cache deterministic `a1_hash64` and `a2_hash64` vectors
  for fast variant matching across Python and MATLAB/Octave.
- MATLAB/Octave reference caches store SNP identifiers and alleles as
  shard-local text payloads while keeping base-pair positions and allele hashes
  as native arrays, improving `load_reference_cache` performance.

### Fixed
- Improved MATLAB/Octave reference and sumstats loading performance.

## [0.2.1] - 2026-05-05

### Added
- LD reference builds now follow an explicit shard-then-finalize workflow. Build
  one Python LD shard at a time, then finalize the Python manifest:

  ```bash
  for shard in 1 2 X; do
    python script/statgen_build_ld.py \
      --bfile /path/to/ref/chr@ --out /path/to/ld_npz \
      --shard "$shard" --scratch /path/to/scratch/shard_"$shard"
  done

  python script/statgen_create_ld_manifest.py --ld /path/to/ld_npz
  ```

- MATLAB/Octave LD distributions use the same shard-then-finalize pattern from
  the Python `.npz` handoff:

  ```matlab
  shards = {'1', '2', 'X'};

  for i = 1:numel(shards)
      statgen.convert_ld_npz_to_mat(ld_npz, ld_mat, shards{i});
  end
  statgen.create_ld_mat_manifest(ld_npz, ld_mat, shards);
  ```

- LD distributions now bundle the matching reference BIM shards. `load_ld(...)`
  can load the paired `ReferencePanel` directly when a reference is not supplied,
  and loaded `LDPanel` objects carry that reference for compatibility checks.
- LD builders now detect monomorphic SNPs and report them in LD metadata. By
  default, builds fail when monomorphic SNPs are present; use
  `--allow-monomorphic-snps` to force output and write a monomorphic SNP report.

### Changed
- `statgen_build_ld.py` no longer exposes PLINK sample-QC controls such as
  `--mind`. Apply sample and variant QC upstream before building the LD
  reference supplied to `statgen`.

### Fixed
- Improved MATLAB/Octave annotation loading performance by replacing the slow
  row-by-row BED payload path with vectorized parsing and painting.

## [0.2.0] - 2026-05-04

### Added
- Summary-statistics loading, alignment, creation, shard subsetting, and cache
  support for Python and MATLAB/Octave.
- BED annotation painting into sparse reference-aligned annotation panels,
  including annotation/shard selection and panel union helpers.
- LD distribution loading and validation for Python `.npz` and MATLAB/Octave
  `.mat` artifacts, including manifest and reference-checksum checks.
- LD panel operations: allele-frequency accessors, shard selection,
  element-wise `r²` multiplication, and fast LD pruning.

## [0.1.0] - 2026-04-30

### Added
- Reference loading for Python and MATLAB/Octave.
- Python and MATLAB/Octave reference cache save/load.
