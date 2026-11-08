# Contigs, Shards, and Panels

## Canonical contig set

`statgen` supports contigs `1`–`22` and `X`, in that order. This is the
canonical contig order used for all shard sequences, row-order keys, and
shard-label validation.

- `Y` and `MT` are recognized non-supported contigs. They are ignored when
  outside a loader's supported output/projection rather than normalized into a
  supported shard. Reference validation excludes `Y`/`MT` rows.
- Labels such as `chr1`/`chrX` indicate non-preprocessed input and are rejected.
- Contig labels are never normalized: `chr1`, `1`, and `NC_000001.11` are
  distinct labels; mismatches are not resolved silently.

`X` is a first-class chromosome. Implementations must not hard-code
chromosomes 1–22 in APIs, manifests, or test fixtures.

## Row-order contract

Reference inputs must already be sorted by:

```text
(chr_rank, bp_numeric)
```

where `chr_rank` follows canonical contig order. Base-pair positions may
repeat; distinct allele pairs within the same `(chr, bp)` may appear in any
order.

This contract applies to both non-sharded files (enforced across the full file)
and sharded files (shard sequence must follow canonical contig order; rows
within each shard must follow row order).

Reference loaders must preserve loaded row order; it defines the aligned SNP
axis. They must not reorder or deduplicate reference rows.

For genotype and sumstats matching, source order is not relevant. Matching
follows [reference.md](reference.md#source-to-reference-matching). Sumstats
loaders do not separately validate source row order or duplicate keys.
Genotype BIM row order and ordered-key duplicate checks are
implementation-dependent.

Reference duplicate-key validation and source ambiguity are defined by the
[matching policy](reference.md#source-to-reference-matching), rather than a
separate requirement to reject all source duplicate keys.

## Shard discovery

Sharded paths use `@` as the shard-label placeholder.

- **Reference loaders**: substitute `@` with canonical contig labels in order
  (no glob-based discovery) and load existing matches.
- **Genotype loaders**: substitute `@` with reference shard labels in order
  when a reference is supplied. Every requested reference shard must have a
  matching bfile shard; users who want a shard subset should first subset the
  reference and load genotype against that smaller reference.
- **LD loaders**: always resolve shard files through `ld_manifest.json`.
  `load_ld` and `load_ld_reference` load the manifest-declared reference cache.
  The optional `shards` parameter further subsets which shards are loaded.
  `load_ld` errors if a requested shard is absent from the reference cache or
  LD panel; `load_ld_reference` errors if a requested shard is absent from the
  reference cache.

Except for sumstats ID mode as defined in [sumstats.md](sumstats.md),
non-sharded single-file inputs are split by the `chr` column into
per-chromosome shards in canonical order. For source loaders with an explicit
`reference` argument, the split is a reference-driven projection: validated
source rows outside the supplied reference shard set are ignored for alignment.
Invalid or ambiguous contig labels remain errors.

## Shard subsetting

All panel objects expose `select_shards(shards)`.
Cache loaders that operate without a required reference
(`load_reference_cache`, `load_annotations_cache`, `load_sumstats_cache`,
`load_genotype_cache`) accept an optional `shards` parameter.
Source loaders with an explicit `reference` argument
(`load_annotations`, `load_annotation`, `load_sumstats`, `load_genotype`) use
the supplied reference shard structure; subsetting is done via `select_shards`
on the reference before passing it in. Such loaders do not accept an additional
`shards` parameter.
`load_ld` follows the cache-loader pattern because the LD distribution includes
its own manifest-declared reference cache. The `shards` argument is honored as
a subset of both that reference cache and the LD panel. `load_ld_reference`
honors the same `shards` argument for the manifest-declared reference cache.
The rules are uniform:

- Omitting `shards` loads or returns all available shards present in the input
  or cache.
- Provided `shards` must be a non-empty list of unique contig labels in
  canonical subsequence order.
- Requesting a shard not present in the input or cache is an error.
