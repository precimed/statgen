# SumstatsShard and Sumstats

## Represents

GWAS summary statistics for one trait or source: Z scores, sample sizes, and
optional effect estimates and quality fields, all oriented to `a1` under the
allele contract in [SPEC.md](SPEC.md). Multiple traits are
independent objects that may share the same reference panel. A future
multi-trait container should use a different name, such as `SumstatsCollection`.

## Disk representation

Portable disk representation is one gzip-compressed TSV file per trait/source:

```text
TRAIT.tsv.gz
```

Canonical sumstats TSV files are produced by external GWAS pipelines, not by
`statgen`; `statgen` consumes these files for loading, alignment, and caching.

Required columns:

- `chr`
- `bp`
- `a1`
- `a2`
- `z`
- `n`

Allele-column contract follows [SPEC.md](SPEC.md).

`z` is required. Files that provide only `beta` and `se` must be converted
upstream before loading into `statgen`.

Optional columns:

- `p`;
- `beta`;
- `se`;
- `eaf`;
- `info`.

The TSV does not have to contain every variant in the reference panel. Loading
against a `ReferencePanel` projects rows into reference order using
`chr:bp:a1:a2` as the exact join key, splits the result into
`SumstatsShard`s matching the reference shards, and represents variants absent
from the TSV as missing values. The loader does not normalize or alias contig
labels; sumstats `chr` values must already match the reference labels.

## In-memory objects

A `SumstatsShard` is aligned to one `ReferenceShard`; a `Sumstats` object is an
ordered collection of `SumstatsShard` objects aligned to a `ReferencePanel`.
Each shard retains the paired reference checksum.

Each `SumstatsShard` contains:

- `zvec`: float vector, length `num_snp`, containing signed Z scores;
- `nvec`: float vector, length `num_snp`, containing effective sample size;
- `logpvec`: float vector, length `num_snp`, containing `-log10(p)` with the
  zero-p convention described below;
- `beta_vec`: optional float vector, length `num_snp`, containing effect sizes
  oriented to `a1`;
- `se_vec`: optional float vector, length `num_snp`, containing standard errors;
- `eaf_vec`: optional float vector, length `num_snp`, containing effect allele
  frequency for `a1` where available;
- `info_vec`: optional float vector, length `num_snp`, containing imputation or
  variant-quality information where available.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: each shard's vector length and row order correspond to the
paired `ReferenceShard`, and any aligned cache is valid only for that
reference.

## Panel-level accessors

`Sumstats` exposes read-only genome-wide accessors that concatenate across
shards in reference panel order. The result is a plain array; downstream code
works with it natively without modifying the `Sumstats` object.

```text
Sumstats.zvec     -> num_snp float vector
Sumstats.nvec     -> num_snp float vector
Sumstats.logpvec  -> num_snp float vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec   -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec  -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
```

SNPs absent from the source file are `NaN`. Missing optional fields use the
language-specific sentinel defined in [SPEC.md](SPEC.md).

## API

```text
load_sumstats(path, reference) -> Sumstats
save_sumstats_cache(sumstats, path)
load_sumstats_cache(path, optional shards) -> Sumstats
create_sumstats(reference, zvec, nvec, optional pvec, optional beta_vec,
                optional se_vec, optional eaf_vec, optional info_vec) -> Sumstats

Sumstats.zvec -> num_snp float vector
Sumstats.nvec -> num_snp float vector
Sumstats.logpvec -> num_snp float vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.select_shards(shards) -> Sumstats
```

Expected behavior:

- `reference` is required; rows are projected into reference order on load.
- `create_sumstats(...)` creates an in-memory `Sumstats` aligned to `reference`
  from full-panel vectors. Canonical sumstats TSV files remain external inputs.
- `chr:bp:a1:a2` joins are exact after basic field parsing; the loader does not
  normalize chromosome labels, swap alleles, or perform strand handling.
- required numeric fields `z` and `n` must parse as finite numeric values for
  all source rows; non-numeric, `NaN`, or infinite values are validation
  errors and must fail load with a clear message.
- cache is a single file (non-sharded); internal per-shard layout is implementation-specific.
- `load_sumstats_cache` performs cache-internal validation only and supports
  optional `shards` subsetting.
- missing variants are represented as `NaN` or masks; row order matches the
  reference panel.
- `logpvec` is derived only from an optional `p` column as `-log10(p)`.
  Missing `p` values yield `NaN`; by convention, `p == 0` yields `Inf`;
  `p < 0` or `p > 1` yields `NaN`.
- for each optional field (`beta_vec`, `se_vec`, `eaf_vec`, `info_vec`):
  when the source column is present, the accessor returns a full aligned vector
  with `NaN` for missing or unmatched rows; when the source column is absent,
  the accessor returns the language-specific missing optional-field sentinel
  from [SPEC.md](SPEC.md) (`None` in Python, `[]` in MATLAB/Octave).
- `create_sumstats(...)` validates vector lengths against `reference.num_snp`.
  Unknown shapes fail clearly; required `zvec` and `nvec` values must satisfy
  the same finite-numeric contract as loaded objects; optional vectors follow
  the same absent/present sentinel semantics.
- cache save/load must preserve the same optional-field semantics (field absent
  remains absent; field present remains a vector).
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native vectors.
- `Sumstats.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md).
- cache payloads store per-shard reference checksums for compatibility checks
  under the general reference-compatibility contract in
  [reference.md](reference.md).
