# Conventions

## Scope

`statgen` owns reusable statistical-genetics data structures and local
transformations:

- reference panels and allele semantics;
- GWAS summary-statistic projection into aligned vectors, assuming inputs
  already follow the reference-genome and allele contract below;
- PLINK bfile genotype shards as inputs and, when needed, on-the-fly genotype
  readers;
- LD panel construction, storage, loading, and cache conversion;
- BED-to-BIM annotation painting, including row-order-preserving masks and
  annotation shards/panels;
- shared Python and MATLAB/Octave loaders plus consistency tests.

`statgen` does not own external data acquisition workflows. Downloading,
versioning, and choosing annotation sources belongs to project-specific
pipelines. Once a BED-like annotation source exists locally, painting that
source onto a reference BIM/SNP table is `statgen` functionality.

`statgen` assumes variant-bearing inputs have already been prepared by
`genomatch` or an equivalent upstream pipeline. In particular, upstream
preprocessing is responsible for genome-build compatibility, contig naming,
reference allele assignment, allele ordering, strand handling, liftover,
intersection, union, and payload projection. `statgen` consumes the resulting
tables as-is and does not rewrite contig labels or allele fields.

Existing `ld/*` functionality conceptually belongs in `statgen`, but the current
code is considered prior art rather than implementation to preserve directly.
The redesigned LD path should be implemented cleanly in `statgen`.

## Input contract

- on-disk inputs must already be aligned to a specific reference genome;
- GRCh38 is the default reference genome unless an object explicitly declares a
  different build;
- contig labels must already use `genomatch` NCBI target naming (`1`-`22`,
  `X`, `Y`, `MT` where applicable), unless an object explicitly declares a
  different `genomatch` contig naming mode;
- input variant records must already use `a1` as the non-reference allele and
  `a2` as the reference allele for that reference genome;
- GWAS signed effects and Z scores must already be oriented to `a1`;
- genotype dosage/count data must count `a1`;
- LD correlations must be signed with respect to `a1` dosages.

`statgen` validates and preserves this contract where practical, but it does
not resolve raw data into the contract. Reference liftover, contig
normalization, allele harmonization, REF/ALT assignment, strand resolution, and
source-specific GWAS cleaning belong to upstream preprocessing tools such as
`genomatch` or project-specific pipelines.

## BIM `cm` column policy

For BIM-like inputs, the `cm` column is a legacy placeholder required only for
6-column schema compatibility (`chr`, `snp`, `cm`, `bp`, `a1`, `a2`).

`statgen` MUST treat `cm` as input-format scaffolding, not biological signal:

- loaders MUST validate BIM row shape includes the `cm` column;
- `cm` MUST NOT be stored in in-memory statgen objects;
- `cm` MUST NOT participate in checksums, compatibility checks, caches, or
  public APIs.

Loaders must not silently convert between chromosome naming schemes. For
example, `chr1`, `1`, and `NC_000001.11` are distinct labels to `statgen`.
If a BED file uses `chr1` and the reference uses `1`, annotation painting must
produce no match or fail validation rather than rewriting either side. The
correct fix is to regenerate or project the input with the selected upstream
contig naming mode before loading it into `statgen`.

## Alignment and missingness

The main rule is that every per-SNP vector or matrix row is aligned to the same
reference-panel SNP order in memory. Rows are never silently dropped from
in-memory aligned objects; missing or filtered values are represented by masks
or `NaN`.

On disk, canonical source-like objects do not have to be expanded to the full
reference panel. They may store only observed rows plus enough metadata to
project onto the reference panel. Aligned full-panel files are allowed and
useful for caches, model-ready objects, and cross-language fixtures, but
alignment is primarily an in-memory contract.

## Reference identity and checksums

Each `ReferenceShard` has a shard-local reference checksum computed from the
ordered sequence of variants in that shard. The checksum input is the UTF-8 text
sequence:

```text
chr:bp:a1:a2\n
```

one line per variant, in shard row order. The `chr` value is the contig label
as loaded from the input after basic field parsing; no chromosome normalization
or alias mapping is applied. `bp` is the decimal integer base-pair coordinate,
and `a1`/`a2` are the allele strings as loaded from the input. The checksum
algorithm is MD5 over that byte sequence, encoded as a lowercase hexadecimal
string. The checksum is intended to identify row order, contig naming,
position, and allele orientation for aligned objects; it is not a security
primitive.

Objects tied to a reference shard store the corresponding shard checksum in
memory and in their caches or portable metadata where applicable. `LDShard`
metadata stores the same checksum for its paired reference shard.

Compatibility checks are exposed through the reference panel:

```text
ReferencePanel.is_object_compatible(object) -> bool
```

`object` is one loaded statgen object expected to be aligned to that reference
panel. The method compares stored shard checksums where available and checks
compatible shard counts and row counts. It returns `true` when the object is
compatible and `false` otherwise. Implementations may log informational
messages or warnings to the console. Compatibility mismatches do not raise
exceptions; callers branch on the returned boolean.

## Shards and panels

A shard is a single chromosome or otherwise local partition. A panel is an
ordered collection of shards. Within each shard, rows define zero-based local
SNP indices. Across a panel, rows define zero-based global SNP indices.

In memory, all panel-like objects are represented as ordered shard vectors.
Simple non-sharded usage is represented as a one-shard panel, including a
single shard that spans multiple chromosomes when the caller deliberately
constructs or loads it that way.

Portable reference, genotype, and LD panels may be sharded or non-sharded on
disk. A one-chromosome analysis is a one-shard panel. A chromosome subset such
as chr20-22 is represented by an ordered subset of shard files. Variant
subsetting creates a new reference panel with the chosen shard structure;
downstream genotype and LD panels are interpreted relative to that new panel.

chrX is a first-class chromosome. Implementations may initially support only
diploid chrX coding inherited from PLINK bfiles, but they must not hard-code
chromosomes 1-22 in APIs, manifests, or test fixtures.

## Storage policy

Portable object files should be language-agnostic. Use:

- PLINK `.bim` for reference shards and PLINK `.bed/.bim/.fam` for genotype
  shards;
- BED for canonical annotation interval inputs;
- `.tsv.gz` for portable non-sharded GWAS summary statistics;
- JSON for manifests and small metadata;
- raw little-endian binary arrays plus JSON metadata for LD triplets.

Avoid pickle, RDS, or MATLAB `.mat` as portable storage. Language-specific
formats are allowed only as caches for performance or interoperability, and
must be reproducible from portable source files.
