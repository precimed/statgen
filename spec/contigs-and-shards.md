# Contigs, Shards, and Panels

## Canonical contig set

`statgen` supports contigs `1`–`22` and `X`, in that order. This is the
canonical contig order used for all shard sequences, row-order keys, and
shard-label validation.

- `Y` and `MT` are recognized but silently ignored on load.
- Labels such as `chr1`/`chrX` indicate non-preprocessed input and are rejected.
- Contig labels are never normalized: `chr1`, `1`, and `NC_000001.11` are
  distinct labels; mismatches are not resolved silently.

`X` is a first-class chromosome. Implementations must not hard-code
chromosomes 1–22 in APIs, manifests, or test fixtures.

## Row-order contract

Source-like SNP tables must already be sorted by:

```text
(chr_rank, bp_numeric, a1_lexicographic, a2_lexicographic)
```

where `chr_rank` follows canonical contig order. Duplicate `(chr, bp, a1, a2)`
tuples are not allowed.

This contract applies to both non-sharded files (enforced across the full file)
and sharded files (shard sequence must follow canonical contig order; rows
within each shard must follow row order).

Loaders validate and fail clearly on violations. They must not reorder or
deduplicate rows.

## Shard discovery

Sharded paths use `@` as the shard-label placeholder.

- **Reference and genotype loaders**: substitute `@` with canonical contig
  labels in order (no glob-based discovery) and load existing matches.
- **LD loaders**: derive expected shard paths from the supplied `reference`
  object rather than performing canonical substitution independently.

Non-sharded single-file inputs are split by the `chr` column into
per-chromosome shards in canonical order.

## Shard subsetting

All panel objects expose `select_shards(shards)`.
Cache loaders that operate without a required reference
(`load_reference_cache`, `load_ld_cache`, `load_annotations_cache`,
`load_sumstats_cache`) accept an optional `shards` parameter.
Source loaders that are defined around an explicit `reference` argument
(`load_ld`, `load_annotations`, `load_sumstats`) use the supplied reference
shard structure; subsetting is done via `select_shards`.
The rules are uniform:

- Omitting `shards` loads or returns all available shards present in the input
  or cache.
- Provided `shards` must be a non-empty list of unique contig labels in
  canonical subsequence order.
- Requesting a shard not present in the input or cache is an error.
