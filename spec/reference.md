# ReferenceShard and ReferencePanel

## Represents

A reference panel is a genome-wide table of SNP positions and alleles that
defines the coordinate system for all other objects. Each SNP is described by
chromosome, base-pair position, identifier, and an ordered allele pair under
the input contract in [SPEC.md](SPEC.md):

- `chr`: upstream contig label, expected to use the selected `genomatch`
  contig naming mode, normally NCBI labels such as `1`-`22` and `X`;
- `snp`: SNP identifier;
- `bp`: base-pair position;
- `a1`, `a2`: ordered alleles under the allele contract in
  [SPEC.md](SPEC.md); `a1` is the non-reference allele, `a2` is
  the reference allele.

Optional columns include `source` and `variant_id`.

## Disk representation

Canonical disk representation is PLINK-style `.bim`:

```text
chr  snp  cm  bp  a1  a2
```

The `cm` column is a legacy placeholder required only for 6-column schema
compatibility. Loaders MUST validate that the `cm` column is present; `cm`
MUST NOT be stored in in-memory objects and MUST NOT participate in checksums,
compatibility checks, caches, or public APIs.

A reference panel may be sharded or non-sharded on disk:

- **Sharded**: one `.bim` file per chromosome via an `@` placeholder in the
  path (e.g. `chr@.bim`).
- **Non-sharded**: a single `.bim` file spanning all chromosomes.

Shard row order defines shard-local SNP indices. Panel row order defines global
SNP indices.

## In-memory objects

A `ReferenceShard` is the in-memory representation of one `.bim` shard.
A `ReferencePanel` is an ordered collection of `ReferenceShard` objects and
defines the multi-chromosome SNP order used for all aligned in-memory objects.
Each shard carries the MD5 reference checksum defined below.

## Reference checksum definition

Each `ReferenceShard` has a shard-local reference checksum. The checksum input
is the UTF-8 text sequence:

```text
chr:bp:a1:a2\n
```

one line per variant in shard row order. The `chr` value is the contig label as
loaded; no normalization is applied. `bp` is the decimal integer base-pair
coordinate; `a1`/`a2` are the allele strings as loaded. The algorithm is MD5
over that byte sequence, encoded as a lowercase hexadecimal string. The checksum
identifies row order, contig naming, position, and allele orientation for
alignment checks; it is not a security primitive.

Regardless of how the panel is sharded, callers access genome-wide vectors and
matrices through read-only panel-level accessors that concatenate across shards
transparently.

## Representation

Implementations may use language-native containers. They must preserve shard
order, expose the required vectors, retain per-shard checksums, and provide
zero-based shard offsets for compatibility with portable metadata.

## API

```text
load_reference(path, optional shards) -> ReferencePanel
save_reference_cache(panel, path)
load_reference_cache(path, optional shards) -> ReferencePanel

ReferencePanel.num_snp -> int
ReferencePanel.chr -> num_snp chromosome-label vector
ReferencePanel.snp -> num_snp string vector
ReferencePanel.bp  -> num_snp integer vector
ReferencePanel.a1  -> num_snp string vector
ReferencePanel.a2  -> num_snp string vector
ReferencePanel.shard_offsets -> table with shard_label, start0, stop0
ReferencePanel.select_shards(shards) -> ReferencePanel
ReferencePanel.is_object_compatible(object) -> bool
```

Expected behavior:

- Reference panels are external inputs; `statgen` never writes `.bim` files.
- Cache is a single file (non-sharded); internal layout is implementation-specific.
- Shard discovery, contig validation, row-order validation, and shard subsetting
  follow [contigs-and-shards.md](contigs-and-shards.md).
- Compute and retain each shard reference checksum.
- Accessors are read-only, concatenate shard columns in reference panel order,
  and return plain language-native vectors or tables.
- `shard_offsets` uses zero-based half-open intervals into genome-wide arrays.
- `load_reference_cache` skips source-style row validation; `shards` subsetting
  applies against cached shard labels per
  [contigs-and-shards.md](contigs-and-shards.md).
- `is_object_compatible` checks whether a loaded statgen object is aligned to
  this reference panel. Compatibility requires the same ordered shard labels,
  matching shard row counts, and matching shard checksums where available.
  Returns `true` when compatible, `false` otherwise. Implementations may log
  informational messages or warnings. Compatibility mismatches do not raise
  exceptions; callers branch on the returned boolean.
