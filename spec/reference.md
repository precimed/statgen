# ReferenceShard and ReferencePanel

## Represents

A reference panel is a genome-wide table of SNP positions and alleles that
defines the coordinate system for all other objects. Each SNP is described by
chromosome, base-pair position, identifier, and an ordered allele pair under
the input contract in [conventions.md](conventions.md):

- `chr`: upstream contig label, expected to use the selected `genomatch`
  contig naming mode, normally NCBI labels such as `1`-`22` and `X`;
- `snp`: SNP identifier;
- `bp`: base-pair position;
- `a1`, `a2`: ordered alleles; `a1` is the non-reference allele, `a2` is the
  reference allele.

Optional columns include `source` and `variant_id`.

## Disk representation

Canonical disk representation is PLINK-style `.bim`:

```text
chr  snp  cm  bp  a1  a2
```

Per [conventions.md](conventions.md), `cm` is validated as part of BIM schema
compatibility and then ignored by `statgen`.

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
load_reference(path) -> ReferencePanel
save_reference_cache(panel, path)
load_reference_cache(path) -> ReferencePanel

ReferencePanel.num_snp -> int
ReferencePanel.chr -> num_snp chromosome-label vector
ReferencePanel.snp -> num_snp string vector
ReferencePanel.bp  -> num_snp integer vector
ReferencePanel.a1  -> num_snp string vector
ReferencePanel.a2  -> num_snp string vector
ReferencePanel.shard_offsets -> table with shard_label, start0, stop0
ReferencePanel.is_object_compatible(object) -> bool
```

Expected behavior:

- Reference panels are external inputs; statgen never writes `.bim` files.
- Cache is a single file (non-sharded); internal layout is implementation-specific.
- If `path` contains `@`: sharded input; discover all matching `.bim` files and
  load each as one chromosome shard, sorted in canonical order.
- If `path` is a single file: split by the `chr` column into per-chromosome
  shards.
- Preserve row order and validate required columns.
- Source reference rows MUST already be sorted by canonical chromosome order
  (`1`-`22`, `X`) and then by ascending `bp` within chromosome.
- `load_reference` MUST validate this ordering and fail clearly on violations.
  It MUST NOT silently reorder rows.
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
  otherwise. Compatibility mismatches are logged as warnings, do not raise
  exceptions, and callers branch on the returned boolean.
- `load_reference_cache` assumes cache payloads were validated at cache
  creation time and should not re-run full source-style row validation (for
  example per-row parsing/sort-order checks) on the default load path.
