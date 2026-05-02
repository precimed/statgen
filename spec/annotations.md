# AnnotationShard and AnnotationPanel

## Represents

Annotations are SNP-level feature masks. A feature may represent functional
annotation, gene membership, geneset membership, or other binary SNP category
definitions.

## Disk representation

Canonical annotation inputs are BED interval files, one file per annotation.
Only BED columns 1–3 are used: chromosome, 0-based start, 0-based exclusive end.
Additional BED columns are ignored. The set of annotations for an analysis is
one BED file per annotation:

```text
annotations/
  baseline.bed
  coding.bed
  conserved.bed
  ...
```

Painted matrices are not canonical—they are derived from BED files and a
specific reference and must be reproducible from them.

## In-memory objects

`AnnotationShard` and `AnnotationPanel` are reference-specific in-memory
objects. An `AnnotationShard` is produced by painting one or more BED interval
files onto a `ReferenceShard`; an `AnnotationPanel` is an ordered collection of
`AnnotationShard` objects aligned to a `ReferencePanel`. The per-shard
representation keeps annotation data partitioned by chromosome, matching the
shard structure of the reference. Each shard retains the paired reference
checksum.

Annotations loaded from BED files are binary. Non-binary annotation values are
out of scope for this phase. LD-weighted or otherwise numeric annotation
matrices are user-space arrays derived from `annomat`, not
`AnnotationPanel` fields.

After painting, each annotation occupies one column in `annomat`. The full set
of annotations is:

- `annomat`: binary matrix, shape `num_snp x num_annot`, one column per
  annotation, rows aligned to the reference shard;
- `annonames`: string vector, length `num_annot`.

## Painting logic

Painting maps a BED interval file onto a `ReferenceShard` or `ReferencePanel`
to produce a 0/1 column in `annomat`. The procedure follows
`annot/paint_bed_to_bim.py`, except that `statgen` does not normalize
chromosome labels. BED chromosome labels and reference `chr` labels must already
use the same upstream contig naming mode, normally `genomatch` NCBI naming.

1. **Exact chromosome matching**: compare the BED chromosome field to the BIM
   `chr` column exactly as loaded. Do not strip `chr`, map `23` to `X`, or
   apply any other alias conversion. Inputs with inconsistent contig naming must
   be fixed upstream before annotation loading.

2. **Interval merging**: within each chromosome, sort intervals by start
   position and merge overlapping or adjacent intervals into a non-overlapping
   sorted array of `(start, end)` pairs.

3. **Position conversion**: BIM `BP` values are 1-based. Convert to 0-based by
   subtracting 1 before testing interval membership.

4. **Mask assignment**: for each SNP, use binary search on the merged start
   array to find the candidate interval. A SNP at 0-based position `p` is
   inside interval `[start, end)` when `start <= p < end`. SNPs on chromosomes
   absent from the BED file receive `0`.

The result is a boolean or 0/1 vector in BIM row order, one entry per SNP.
Complement masks are a user-space operation on the returned `annomat`; the
annotation loader does not provide a negation option.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: row order and SNP count correspond to the paired reference,
and any cache is valid only for that reference.
Internal representation should default to sparse storage. `annomat` may be
exposed as dense or sparse depending on implementation heuristics (for example
mask density), but semantics are identical. Implementations should avoid
materializing dense `num_snp × num_annot` arrays unless explicitly requested.

## Panel-level accessors

`AnnotationPanel` is immutable after loading. It exposes read-only genome-wide
accessors that concatenate across shards in reference panel order. Downstream
code works with the resulting plain arrays natively.

```text
AnnotationPanel.annomat    -> num_snp × num_annot binary matrix
AnnotationPanel.annonames  -> num_annot string vector
```

`annonames` is identical across shards and returned once. The `annomat` type
follows the shard representation: `bool`, `uint8`, or sparse binary matrix.

## API

```text
load_annotations(bed_paths, reference) -> AnnotationPanel
save_annotations_cache(panel, path)
load_annotations_cache(path, optional shards) -> AnnotationPanel
create_annotations(reference, annomat, annonames) -> AnnotationPanel
create_annotation(reference, annovec, annoname) -> AnnotationPanel

AnnotationPanel.annomat -> num_snp × num_annot binary matrix
AnnotationPanel.annonames -> num_annot string vector
AnnotationPanel.select_shards(shards) -> AnnotationPanel
AnnotationPanel.select_annotations(names) -> AnnotationPanel
AnnotationPanel.union_annotations(other, optional mode) -> AnnotationPanel
```

Expected behavior:

- `bed_paths` is a list of BED files, one per annotation; annotation names are
  derived deterministically from BED basenames (without extension).
- `reference` is required; the output shard structure matches the reference
  panel; `annomat` rows are aligned to the reference.
- `annonames` must be unique. Duplicate names from BED basenames or
  programmatic inputs are an error.
- caching saves and restores the painted `annomat`; the cache is a single file
  (non-sharded). `load_annotations_cache` performs cache-internal validation
  only and supports optional `shards` subsetting.
- internal per-shard layout within the cache file is implementation-specific.
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native matrices or vectors.
- `annonames` must be identical across shards and is returned once.
- `AnnotationPanel.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md).
- `create_annotations` and `create_annotation` validate strict alignment to the
  supplied `reference` (`num_snp` shape match, binary values only) and return
  immutable `AnnotationPanel` objects.
- `AnnotationPanel.select_annotations(names)` preserves requested name order and
  fails on unknown names (no silent drops).
- `AnnotationPanel.union_annotations(other, optional mode)` requires
  `ReferencePanel.is_object_compatible(other) == true`. `mode` defaults to
  `by_name`; name collisions are errors.
- cache payloads store per-shard reference checksums so compatibility with a
  `ReferencePanel` can be checked after load via
  `ReferencePanel.is_object_compatible`.
