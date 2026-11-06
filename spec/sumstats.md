# SumstatsShard and Sumstats

## Represents

GWAS summary statistics for one trait or source: p-values, optional Z scores,
sample sizes, effect estimates, and quality fields, aligned to reference `a1` under
the allele contract in [SPEC.md](SPEC.md). Multiple traits are independent
objects that may share the same reference panel. A future multi-trait container
should use a different name, such as `SumstatsCollection`.

## Disk representation

Portable disk representation is one gzip-compressed TSV file per trait/source:

```text
TRAIT.tsv.gz
```

Canonical sumstats TSV files are produced by `genomatch` or equivalent
post-harmonization pipelines, not by `statgen`; `statgen` consumes these files
for loading, alignment, and caching.

Required columns:

- `chr` or `CHR`
- `bp` or `POS`
- `a1` or `EffectAllele`
- `a2` or `OtherAllele`
- `p` or `P`

Allele roles follow [SPEC.md](SPEC.md). An explicit genomic `REF` column is
not currently supported for allele interpretation or matching; source orientation is
determined by the effect and other allele columns.

The `p` column is required, but individual values may be missing. Non-missing
values must be finite and in `[0, 1]`; `p == 0` is allowed, and out-of-range
values are validation errors.

Optional columns:

- `z`;
- `n`;
- `beta`;
- `se`;
- `eaf`;
- `info`.

Column recognition is case-insensitive. Apart from case, the recognized
genomatch cleaned-sumstats names `POS`, `EffectAllele`, and `OtherAllele` map
to internal fields `bp`, `a1`, and `a2`. Genomatch vmap-style `bp`, `a1`, and
`a2` are accepted directly. `SNP` may be present but is not used for matching.
Unrecognized extra columns are ignored.

Logical column types:

| Logical column | Required | Type |
| --- | --- | --- |
| `chr` | yes | string |
| `bp` | yes | float-parsed integer |
| `a1` | yes | string |
| `a2` | yes | string |
| `p` | yes | float |
| `snp` | no | string |
| `z` | no | float |
| `n` | no | float |
| `beta` | no | float |
| `se` | no | float |
| `eaf` | no | float |
| `info` | no | float |

`bp` is parsed as a numeric value and then validated to be integral. It is not
parsed as a string on performance-sensitive paths. Empty fields, `NA`, and
`NaN` are missing values in every numeric column and are represented as `NaN`
after parsing. Implementations may accept additional missing-value spellings
and may read infinities as missing; any other non-numeric value is a
validation error. Missing `bp` is invalid.

Implementations must parse numeric logical columns with native numeric TSV
reader support where the runtime provides it. MATLAB/Octave readers therefore
use `%f` textscan specifiers for `bp`, `p`, `z`, `n`, `beta`, `se`, `eaf`, and
`info`, and `%s` for string columns. This avoids a default path that reads large
numeric columns as strings and converts them later row by row. Python and R may
use their language-native tabular readers and vectorized numeric conversion
paths, provided they preserve the same logical column types and validation
behavior.

The TSV does not have to contain every variant in the reference panel.
`load_sumstats` is a reference-driven projection from a possibly larger source.
Loading against a `ReferencePanel` projects rows into reference order under
[the shared matching contract](reference.md#source-to-reference-matching).
Source rows with missing `p` are validated like other rows, then dropped before
matching; matching and the warning below use only the remaining rows.
Validated source rows on canonical contigs outside the supplied reference
shard set are ignored for alignment, such as `X` rows loaded against a
chr1-only reference. Recognized non-supported contigs, invalid contig labels,
and allele syntax follow the shared rules in
[contigs-and-shards.md](contigs-and-shards.md) and [SPEC.md](SPEC.md). The
result is split into `SumstatsShard`s matching the reference shards, and
variants absent from the TSV are represented as missing values. If a supplied
reference shard has no matching source rows at all, the loader warns and
represents that shard as all missing. The loader does not normalize or alias
contig labels; matched sumstats `chr` values must already match the reference
labels.

For a swapped match, negate source `z` and `beta`, and replace source `eaf`
with `1 - eaf`. Leave `p`, `n`, `se`, and `info` unchanged. Preserve missing
optional values and absent optional fields. Exact matches are unchanged.
`create_sumstats` accepts already reference-oriented vectors and does not
perform matching or swaps. Sumstats objects store aligned values without a
per-variant swap mask.

## In-memory objects

A `SumstatsShard` is aligned to one `ReferenceShard`; a `Sumstats` object is an
ordered collection of `SumstatsShard` objects aligned to a `ReferencePanel`.
Each shard retains the paired reference checksum as `reference_checksum`.

Each `SumstatsShard` contains:

- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`;
- `logpvec`: float vector, length `num_snp`, containing `-log10(p)` with the
  zero-p convention described below;
- `zvec`: optional float vector, length `num_snp`, containing signed Z scores;
- `nvec`: optional float vector, length `num_snp`, containing effective sample
  size;
- `is_present`: derived logical vector, length `num_snp`, equal to
  `~isnan(logpvec)`;
- `beta_vec`: optional float vector, length `num_snp`, containing effect sizes
  oriented to reference `a1`;
- `se_vec`: optional float vector, length `num_snp`, containing standard errors;
- `eaf_vec`: optional float vector, length `num_snp`, containing the frequency
  of reference `a1` (ALT) where available;
- `info_vec`: optional float vector, length `num_snp`, containing imputation or
  variant-quality information where available.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: each shard's vector length and row order correspond to the
paired `ReferenceShard`, and any aligned cache is valid only for that
reference.

## Cache layout

MATLAB/Octave sumstats caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
logpvec
zvec
nvec
beta_vec
se_vec
eaf_vec
info_vec
```

`metadata` is a struct with:

```text
schema = "sumstats_cache/0.1"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
has_z
has_n
has_beta
has_se
has_eaf
has_info
```

`logpvec` is a required panel-wide numeric vector. Optional vectors are also
panel-wide numeric vectors when present; absent optional fields are saved as
`[]` and indicated by the corresponding `has_*` flag. There is no `has_p`
metadata field because `p` is represented by required `logpvec`. `is_present`
is derived from `logpvec` and is not stored. Shard offsets are zero-based
half-open intervals into the panel-wide vectors and are sufficient to
reconstruct `SumstatsShard` objects. `shard_checksums` stores the per-shard
reference checksums used to restore `SumstatsShard.reference_checksum`. Cache
metadata validation should be cheap, depending on shard count and vector
dimensions rather than scanning all SNP values. After structural validation
passes, cache loaders must perform the warning checks for optional
`zvec`/`nvec` completeness; those warning checks inspect vector values as
needed and are not part of the structural compatibility gate.

R sumstats caches are RDS files containing one named list with the same logical
top-level fields and metadata fields as the MATLAB/Octave cache layout above.
R `metadata` is a named list. Present vectors are panel-wide numeric vectors;
absent optional fields are stored as `NULL` and indicated by the corresponding
`has_*` flag. Optional top-level field names remain present in the RDS list even
when their value is `NULL`.

## Panel-level accessors

`Sumstats` exposes read-only genome-wide accessors that concatenate across
shards in reference panel order. The result is a plain array; downstream code
works with it natively without modifying the `Sumstats` object.

```text
Sumstats.logpvec  -> num_snp float vector
Sumstats.zvec     -> num_snp float vector, or missing optional field sentinel
Sumstats.nvec     -> num_snp float vector, or missing optional field sentinel
Sumstats.is_present -> num_snp logical vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec   -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec  -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
```

SNPs without a matched source row have `NaN` in `logpvec` and any present
optional vectors, with `is_present == false`. Where `p` passed to
`create_sumstats` is missing (`NaN`, or `NA` in R), `logpvec` and present
optional vectors are `NaN`.
Missing optional statistics alone do not make a variant with valid `p` absent.
Missing optional fields use the language-specific sentinel defined in
[SPEC.md](SPEC.md).

## API

```text
load_sumstats(path, reference, optional allow_indel_swaps=false) -> Sumstats
save_sumstats_cache(sumstats, path)
load_sumstats_cache(path, optional shards) -> Sumstats
create_sumstats(reference, p, optional z, optional n, optional beta,
                optional se, optional eaf, optional info) -> Sumstats

Sumstats.num_snp -> int
Sumstats.logpvec -> num_snp float vector
Sumstats.zvec -> num_snp float vector, or missing optional field sentinel
Sumstats.nvec -> num_snp float vector, or missing optional field sentinel
Sumstats.is_present -> num_snp logical vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.save_cache(path) -> void
Sumstats.select_shards(shards) -> Sumstats
```

Expected behavior:

- `reference` is required; rows are projected into reference order on load.
- `create_sumstats(...)` creates an in-memory `Sumstats` aligned to `reference`
  from full-panel vectors. Canonical sumstats TSV files remain external inputs.
- sumstats TSV column recognition is case-insensitive and only accepts the
  internal field names plus `POS`, `EffectAllele`, and `OtherAllele` from the
  genomatch cleaned-sumstats schema.
- sumstats-to-reference matching follows the shared source-to-reference
  matching contract in [reference.md](reference.md#source-to-reference-matching).
  The loader does not normalize chromosome labels or perform strand handling.
- `p` follows the validation and missing-value rules above. Invalid values
  must fail load with a clear message. At `info` verbosity, `load_sumstats`
  reports the number of rows dropped for missing `p` when nonzero; `quiet`
  suppresses this report.
- optional numeric fields `z` and `n` may be absent. `load_sumstats` warns when
  either column is absent. When present, each returns a full aligned vector with
  `NaN` for unmatched rows and missing source values, and `load_sumstats` warns
  if any value is missing among aligned variants with `is_present == true`.
- cache is a single file (non-sharded). The MATLAB/Octave cache layout is
  specified in the "Cache layout" section above; absent optional fields are
  serialized as empty arrays (`[]`) so all field variables are always present
  in the `.mat` file, distinguished by the `has_*` flags in metadata.
- `load_sumstats_cache` performs cache-internal validation only and supports
  optional `shards` subsetting; per-shard checksums are trusted from cache
  metadata. It warns when cached `zvec` or `nvec` is absent, and warns when a
  present `zvec` or `nvec` has missing values among cached variants with
  `is_present == true`.
- missing variants are represented as `NaN` or masks; row order matches the
  reference panel.
- `logpvec` is derived from the required `p` column as `-log10(p)`. By
  convention, `p == 0` yields `Inf`.
- `is_present` is a derived accessor defined as `~isnan(logpvec)`, including
  `true` for `p == 0`.
- for each optional field (`zvec`, `nvec`, `beta_vec`, `se_vec`, `eaf_vec`,
  `info_vec`):
  when the source column is present, the accessor returns a full aligned vector
  with `NaN` for missing or unmatched rows; when the source column is absent,
  the accessor returns the language-specific missing optional-field sentinel
  from [SPEC.md](SPEC.md) (`None` in Python, `[]` in MATLAB/Octave, `NULL` in
  R).
- `create_sumstats(...)` validates vector lengths against `reference.num_snp`.
  Unknown shapes fail clearly. `p` accepts `NaN` (and `NA` in R) as missing;
  other values must be finite and in `[0, 1]`. Optional vectors follow the same
  absent/present sentinel semantics.
- cache save/load must preserve the same optional-field semantics (field absent
  remains absent; field present remains a vector).
- `Sumstats.save_cache(...)` is a thin convenience method equivalent to
  `save_sumstats_cache(sumstats, ...)`.
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native vectors.
- `Sumstats.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md).
- cache payloads store per-shard reference checksums for compatibility checks
  under the general reference-compatibility contract in
  [reference.md](reference.md).
