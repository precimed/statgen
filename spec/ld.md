# LDShard and LDPanel

## Represents

LD describes the pairwise signed correlation structure between variants in the
same processed shard. Values are signed correlation `r`, not `r²`; consumers
that need `r²` compute it after loading. Cross-chromosome LD is not stored.
Signs are oriented to `a1` dosages, consistent with the allele contract in
[SPEC.md](SPEC.md).

LD is paired with a reference: row and column order, SNP count, and the
accompanying `mafvec` all correspond to a specific `ReferenceShard`.

## Disk representation

Portable LD shard storage is raw binary arrays plus JSON metadata. An LD shard
is tied to a specific reference shard: its row order, SNP count, and `mafvec`
entries correspond to the paired `ReferenceShard`, identified by checksum in
`metadata.json`.

### Autosomal shard (single sex-agnostic shard per chromosome)

```text
ld/chrN/
  metadata.json
  ld_idx1.i32
  ld_idx2.i32
  ld_r.f32
  mafvec.f32
```

### chrX shard group (multiple sex-specific shards in one directory)

When sex-specific chrX shards are present, the chromosome directory acts as a
**shard group**: the top-level `metadata.json` has `object_type: "ld_shard_group"`
and lists the available sex labels. Each sex-specific shard lives in a named
subdirectory:

```text
ld/chrX/
  metadata.json          ← object_type: "ld_shard_group", lists sex labels
  combined/
    metadata.json        ← sex: "combined"
    ld_idx1.i32
    ld_idx2.i32
    ld_r.f32
    mafvec.f32
  male/
    metadata.json        ← sex: "male"
    ld_idx1.i32
    ld_idx2.i32
    ld_r.f32
    mafvec.f32
  female/
    metadata.json        ← sex: "female"
    ld_idx1.i32
    ld_idx2.i32
    ld_r.f32
    mafvec.f32
```

Not all sex labels need to be present; the group `metadata.json` is the
authoritative list of which subdirectories exist.

### Required payload (per shard directory)

- `ld_idx1.i32`: little-endian signed int32, zero-based within-shard row
  indices for upper-triangle entries;
- `ld_idx2.i32`: little-endian signed int32, zero-based within-shard column
  indices for upper-triangle entries;
- `ld_r.f32`: little-endian IEEE-754 float32, signed LD correlations relative
  to `a1` dosage;
- `mafvec.f32`: little-endian IEEE-754 float32, allele frequencies aligned to
  the paired reference shard.

`posvec` and `snpvec` must not be stored in the LD shard. SNP identity,
position, and row order come from the paired `ReferenceShard` `.bim`.

### `metadata.json` for a plain shard

```json
{
  "object_type": "ld_shard",
  "schema_version": "0.1",
  "chr": "N",
  "sex": null,
  "num_snp": 123,
  "num_ld": 456,
  "index_base": 0,
  "byte_order": "little_endian",
  "triangle": "upper",
  "diagonal": "implicit_unit",
  "value": "r",
  "reference_checksum": "...",
  "build_tool": "plink2",
  "build_command": "plink2 ...",
  "ld_window_kb": 10000,
  "ld_r2_threshold": 0.05,
  "num_sample": 1000
}
```

`sex` is `null` for autosomal shards and one of `"combined"`, `"male"`, or
`"female"` for sex-specific shards.

### `metadata.json` for a shard group

```json
{
  "object_type": "ld_shard_group",
  "schema_version": "0.1",
  "chr": "X",
  "sex_shards": ["combined", "male", "female"]
}
```

Little-endian encoding is required for cross-platform safety. Loaders must
byteswap on big-endian platforms rather than silently reading native-endian raw
arrays.

An `LDPanel` is stored as an ordered collection of LD shard directories
(plain or shard-group).

## Building an LD panel

`statgen_build_ld.py` builds an LD panel from a PLINK bfile. It accepts both
input layouts:

- **Sharded bfile input** (`@` in path, e.g. `chr@`): one bfile per
  chromosome; discovery is by canonical substitution (`1`-`22`, `X`) and each
  discovered shard is processed independently to produce one shard directory.
- **Non-sharded bfile input** (single path, no `@`): the script splits by
  chromosome and produces one shard directory per chromosome.

For each processed unit (chromosome shard, or chromosome extracted from a
non-sharded bfile) the steps are:

1. **Compute allele frequencies** using PLINK2 `--freq`. The output provides
   `a1f` aligned to the BIM row order. Any installed PLINK2 version is
   acceptable; the exact command line should be recorded in metadata.

2. **Compute pairwise LD** using PLINK2 with `--r` and `--keep-allele-order`
   with a window
   of 10,000 kb and an `r²` threshold of 0.05 (pairs with `|r| < sqrt(0.05)`
   are omitted). `--keep-allele-order` ensures signs are consistent with `a1`.

3. **Convert to portable format**: write binary triplets (`ld_idx1.i32`,
   `ld_idx2.i32`, `ld_r.f32`, `mafvec.f32`) and `metadata.json` into the
   output shard directory.

The r² threshold, window size, PLINK2 command, and sample count are build
metadata and must be recorded in `metadata.json`. Do not construct the binary
files by hand.

The sharding of the resulting LD panel must match the sharding of the reference
used with it.

### chrX sex-specific build

By default `statgen_build_ld.py` produces sex-specific shards for chrX: one
combined shard (all samples) and one male-only and one female-only shard. These
are written as a shard group:

```text
out/chrX/
  metadata.json          ← object_type: "ld_shard_group"
  combined/  male/  female/
```

Sex is determined from the FAM file column 5 (1 = male, 2 = female). The script
requires a FAM file with non-missing sex for chrX builds. To suppress sex
splitting and write only a single combined shard for chrX, pass `--no-sex-split`.

For chr1–22, the script always writes a single sex-agnostic (combined) shard
regardless of `--no-sex-split`.

## In-memory objects

An `LDShard` holds the upper-triangle triplets for one specific LD matrix and
the associated `mafvec`. Each `LDShard` is unambiguously identified by its
`(chr, sex)` pair.

An `LDPanel` is an ordered collection of shard lists matching the paired
reference panel. For autosomal chromosome shards, each chromosome has exactly
one `LDShard` (sex=`None`). For chrX, the panel may hold any non-empty subset
of the available sex-specific shards (`"combined"`, `"male"`, `"female"`).
Loading always loads the full set of shards present on disk for each shard.

`LDPanel` carries a `chrX_default_sex` field (`"combined"`, `"male"`, or
`"female"`) that selects which chrX shard downstream operations use when they
need a single LD matrix for chrX. It must name a shard actually present in the
panel. For panels built by `statgen_build_ld.py`, the default is `"combined"`.

`LDShard` fields:

- `chr`: contig/shard label (`1`-`22` or `X`);
- `sex`: `None` for autosomal/sex-agnostic shards; `"combined"`, `"male"`, or
  `"female"` for sex-specific chrX shards;
- `num_snp`: number of SNPs in the shard;
- `ld_idx1`: int32 vector of zero-based row indices;
- `ld_idx2`: int32 vector of zero-based column indices;
- `ld_r`: float32 vector of signed LD correlations;
- `mafvec`: float32 vector of allele frequencies aligned to the paired
  `ReferenceShard`;
- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`.

The triplets store the upper triangle only and exclude the diagonal. Required
invariant: `ld_idx1 < ld_idx2`. Loaders materialize a symmetric sparse matrix
and set the diagonal to `1.0`. Missing off-diagonal entries are interpreted as
zero.

## Representation

Implementations may use language-native sparse matrix containers. The object is
tied to one reference panel: sparse matrix row/column order and `mafvec`
entries correspond to the paired reference shard, and any cache is valid only
for that reference.

## API

```text
load_ld(ld_dir, reference, optional chrX_default_sex) -> LDPanel
save_ld_cache(panel, path)
load_ld_cache(path, optional shards, optional chrX_default_sex) -> LDPanel

LDPanel.mafvec(optional sex) -> num_snp float vector
LDPanel.chrX_default_sex -> "combined" | "male" | "female"
LDPanel.select_shards(shards) -> LDPanel
LDPanel.multiply_r2(M, optional sex) -> vector or matrix with same shape as M
fast_prune(logpvec, ld_panel, optional r2_threshold, optional sex) -> logpvec
```

`ld_dir` may identify a panel root containing one or more shard directories, or
a single shard directory. For a panel root, `load_ld` derives expected shard
directories from the supplied reference shard labels and loads
`ld_dir/<shard_label>/`.

`path` for `save_ld_cache` and `load_ld_cache` must contain `@`, replaced by
the shard label. One cache file is written or read per shard or shard group.

Expected behavior:

- `reference` is required; the loader retains on-disk shard checksums and may
  warn when they do not match the supplied reference. It must validate that all
  shards are internally consistent with their metadata, including matching
  `num_snp` and binary array lengths.
- for `load_ld(ld_dir, reference, ...)`, every shard declared by the supplied
  reference must be present under `ld_dir/<shard_label>/`; missing expected
  shards are errors.
- when a chromosome directory contains `metadata.json` with
  `object_type: "ld_shard_group"`, all listed sex-specific subdirectories are
  loaded and stored as a per-chromosome shard list; otherwise the directory is
  loaded as a single sex-agnostic shard.
- `chrX_default_sex` is stored on the returned `LDPanel` and must name a shard
  present in the panel; the loader raises an error if it does not.
- caching saves and restores all shards for the panel (the full shard group for
  chrX, the single shard for autosomes); `chrX_default_sex` is stored in the
  cache and restored on load.
- `load_ld_cache` does not require a reference object and performs cache-internal
  integrity validation only.
- `load_ld_cache` and `LDPanel.select_shards` shard subsetting follow
  [contigs-and-shards.md](contigs-and-shards.md).
- the cache stores per-shard reference checksums so compatibility with a
  `ReferencePanel` can be checked after load via
  `ReferencePanel.is_object_compatible`.
- Panel accessors are read-only and concatenate shard data in reference panel
  order.
- `chrX_default_sex` defaults to `"combined"` when omitted.
- `r2_threshold` defaults to `0.2` when omitted.
- `LDPanel.mafvec(optional sex)` uses `chrX_default_sex` for chrX by default; a
  sex label overrides chrX selection and must name a shard present in the panel.
- `LDPanel.multiply_r2(M, optional sex)` computes `LD_r² * M` independently per
  reference shard, never materializes a dense genome-wide LD matrix, accepts
  both 1-D and 2-D inputs, and returns the same shape as `M`.
- `fast_prune` applies greedy significance-based pruning independently per
  reference shard using the LD selected by `sex`.
- Shard-level matrix multiplication is an internal implementation detail, not
  part of the public API.

## r² Matrix Multiply

The core LD operation multiplies the per-shard squared LD matrix by a
genome-wide matrix or vector. This is the only numerical primitive exposed by
`LDPanel`; all downstream computations (LD scores, annotation LD weighting,
per-SNP variance) are expressed in terms of it.

```text
LDPanel.multiply_r2(M, optional sex) -> vector or matrix with same shape as M
```

`M` is a `num_snp × k` matrix or a `num_snp` vector aligned to the reference
panel row order. The result has the same shape as `M`.

`LDPanel.multiply_r2(M, optional sex)` splits `M` by reference shard, computes
`LD_r² * M_shard` for each shard, and concatenates the results in reference
panel order. Omitting `sex` uses `chrX_default_sex` for the chrX shard; pass
`"male"`, `"female"`, or `"combined"` to override in language-specific syntax.

The operation must not materialize a dense genome-wide LD matrix. Callers
pre-exclude SNPs by setting the corresponding rows of `M` to zero.

## Pruning

```text
fast_prune(logpvec, ld_panel, optional r2_threshold, optional sex) -> logpvec
```

Greedy significance-based LD pruning. Input and output are genome-wide vectors
aligned to the reference panel.

Algorithm (applied independently per reference shard):

1. Stable-sort SNPs by `|logpvec|` descending; ties keep reference order. Skip
   `NaN` entries (treated as pre-excluded).
2. For each SNP in sorted order: if not already pruned, retain the current SNP
   and mark all other SNPs with `r² ≥ r2_threshold` as pruned (output set to
   `NaN`). The current SNP is never pruned by its own diagonal.
3. Retained positions keep their original `logpvec` values.

Omitting `sex` uses `chrX_default_sex` for the chrX shard. Omitting
`r2_threshold` uses `0.2`.
