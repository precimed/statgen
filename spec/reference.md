# ReferenceShard and ReferencePanel

## Represents

A reference panel is a genome-wide table of SNP positions and alleles that
defines the coordinate system for all other objects. Each SNP is described by
chromosome, base-pair position, identifier, genetic position, and an ordered
allele pair under the input contract in [conventions.md](conventions.md):

- `chr`: upstream contig label, expected to use the selected `genomatch`
  contig naming mode, normally NCBI labels such as `1`-`22` and `X`;
- `snp`: SNP identifier;
- `cm`: genetic position, retained for BIM compatibility;
- `bp`: base-pair position;
- `a1`, `a2`: ordered alleles; `a1` is the non-reference allele, `a2` is the
  reference allele.

Optional columns include `source` and `variant_id`.

The `cm` column is loaded and retained in shard metadata for BIM compatibility,
but it is not part of the standard panel-level accessor API.

## Disk representation

Canonical disk representation is PLINK-style `.bim`:

```text
chr  snp  cm  bp  a1  a2
```

A reference panel may be sharded or non-sharded on disk:

- **Sharded**: one `.bim` file per chromosome, addressed via an `@` placeholder
  in the path (e.g. `chr@.bim` expands to `chr1.bim`, `chr2.bim`, …, `chrX.bim`).
  Canonical chromosome order is 1–22, X. A chromosome subset such as chr20–22
  is represented by those shard files only.
- **Non-sharded**: a single `.bim` file spanning all chromosomes.

Shard row order defines shard-local SNP indices. Panel row order defines global
SNP indices.

## In-memory objects

A `ReferenceShard` is the in-memory representation of one `.bim` shard.
A `ReferencePanel` is an ordered collection of `ReferenceShard` objects and
defines the multi-chromosome SNP order used for all aligned in-memory objects.
Each shard carries the MD5 reference checksum defined in
[conventions.md](conventions.md).

Regardless of how the panel is sharded, callers access genome-wide vectors and
matrices through read-only panel-level accessors that concatenate across shards
transparently.

## Representation

Implementations may use language-native containers. They must preserve shard
order, expose the required vectors, retain per-shard checksums, and provide
zero-based shard offsets for compatibility with portable metadata.

## API

```text
load_reference(path, no_shard=False) -> ReferencePanel
save_reference_cache(panel, path)
load_reference_cache(path) -> ReferencePanel

ReferencePanel.num_snp -> int
ReferencePanel.chr -> num_snp chromosome-label vector
ReferencePanel.snp -> num_snp string vector
ReferencePanel.bp  -> num_snp integer vector
ReferencePanel.a1  -> num_snp string vector
ReferencePanel.a2  -> num_snp string vector
ReferencePanel.shard_offsets -> table with shard_label, start0, stop0
ReferencePanel.is_object_compatible(object, optional strict) -> bool
```

Expected behavior:

- Reference panels are external inputs; statgen never writes `.bim` files.
- Cache is a single file (non-sharded); internal layout is implementation-specific.
- If `path` contains `@`: sharded input; one `.bim` per chromosome; `no_shard`
  is ignored.
- If `path` is a single file and `no_shard=False` (default): split by the `chr`
  column into per-chromosome shards.
- If `path` is a single file and `no_shard=True`: load as a single shard
  spanning all chromosomes, with shard label `"all"`.
- Preserve row order and validate required columns.
- Compute and retain each shard reference checksum.
- Assume the upstream allele contract has already been satisfied; loaders
  validate syntax and preserve allele order, but `.bim` alone is not sufficient
  to prove reference/alternate allele truth.
- Accessors are read-only, concatenate shard columns in reference panel order,
  and return plain language-native vectors or tables.
- `shard_offsets` uses zero-based half-open intervals into genome-wide arrays.
- `is_object_compatible` checks one loaded statgen object against this
  reference panel by comparing shard counts, row counts, and stored shard
  checksums where available. It returns `true` when compatible and `false`
  otherwise. Omitting `strict` warns on checksum mismatch; strict validation
  logs mismatches as errors. Compatibility mismatches do not raise exceptions;
  callers branch on the returned boolean.
