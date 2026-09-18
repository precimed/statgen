# Contigs, Shards, and Panels

## Canonical contig set

`statgen` supports contigs `1`–`22` and `X`, in that order. This is the
canonical contig order used for all shard sequences, row-order keys, and
shard-label validation.

- `Y` and `MT` are recognized non-supported contigs. They are ignored when
  outside a loader's supported output/projection rather than normalized into a
  supported shard. Source loaders drop `Y`/`MT` rows before row-order and
  duplicate-key validation for supported output shards.
- Labels such as `chr1`/`chrX` indicate non-preprocessed input and are rejected.
- Contig labels are never normalized: `chr1`, `1`, and `NC_000001.11` are
  distinct labels; mismatches are not resolved silently.

`X` is a first-class chromosome. Implementations must not hard-code
chromosomes 1–22 in APIs, manifests, or test fixtures.

## Row-order contract

Reference and genotype BIM inputs must already be sorted by:

```text
(chr_rank, bp_numeric)
```

where `chr_rank` follows canonical contig order. Base-pair positions may
repeat; distinct allele pairs within the same `(chr, bp)` may appear in any
order.

This contract applies to both non-sharded files (enforced across the full file)
and sharded files (shard sequence must follow canonical contig order; rows
within each shard must follow row order).

Loaders must not reorder or deduplicate BIM rows. Loaded reference order defines
the aligned SNP axis, and genotype BIM/BED/PLOIDY row correspondence must remain
intact.

Sumstats source rows may appear in any order. Loading projects them onto the
reference axis using the ordered allele identity keys in
[reference.md](reference.md), independently of source row order.

## Duplicate-key validation

Duplicate variant matching keys, as defined in [reference.md](reference.md),
must fail within the following scope:

- **BIM inputs**: all supported-contig rows in each loaded file, including
  chromosomes outside the selected reference shard set.
- **Sumstats inputs**: all rows on chromosomes selected by the supplied
  reference, including variants absent from that reference. Duplicate keys on
  unselected chromosomes are ignored.

Loaders must not silently deduplicate these rows or choose one occurrence.

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

Non-sharded single-file inputs are split by the `chr` column into
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
