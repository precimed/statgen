# SumstatsShard and Sumstats

## Represents

GWAS summary statistics for one trait or source: Z scores, sample sizes, and
optional effect estimates and quality fields, all oriented to `a1` under the
allele contract in [conventions.md](conventions.md). Multiple traits are
independent objects that may share the same reference panel. A future
multi-trait container should use a different name, such as `SumstatsCollection`.

## Disk representation

Portable disk representation is one gzip-compressed TSV file per trait/source:

```text
TRAIT.tsv.gz
```

Required columns:

- `chr`
- `bp`
- `a1`
- `a2`
- `z`
- `n`

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

The input contract requires no duplicate `chr:bp:a1:a2` keys. Loaders are not
required to validate or resolve duplicate keys; behavior is undefined when this
contract is violated.

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
load_sumstats_cache(path, reference) -> Sumstats

Sumstats.zvec -> num_snp float vector
Sumstats.nvec -> num_snp float vector
Sumstats.logpvec -> num_snp float vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
```

Expected behavior:

- `reference` is required; rows are projected into reference order on load.
- `chr:bp:a1:a2` joins are exact after basic field parsing; the loader does not
  normalize chromosome labels, swap alleles, or perform strand handling.
- cache is a single file (non-sharded); internal per-shard layout is implementation-specific.
- missing variants are represented as `NaN` or masks; row order matches the
  reference panel.
- `logpvec` is derived only from an optional `p` column as `-log10(p)`.
  Missing `p` values yield `NaN`; by convention, `p == 0` yields `Inf`;
  `p < 0` or `p > 1` yields `NaN`.
- optional vectors are loaded when the corresponding source columns are
  present and otherwise omitted or filled with `NaN` in aligned caches.
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native vectors.
