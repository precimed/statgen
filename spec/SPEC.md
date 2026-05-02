# statgen specification

This directory defines the object-centered design for `statgen`, a shared
Python and MATLAB/Octave package for statistical genetics data structures and
the basic operations on them needed for downstream inference.

The design is centered on individual objects. A project may have multiple
genotype panels, LD panels, GWAS summary-statistic objects, annotation panels,
and reference panels in the same workspace. Object loaders must not require all
data types to live under a single bundle root.

## Package-wide rules

In this file:

1. [API notation](#api-notation)
1. [Contig naming, Genome build and Allele contract](#contig-naming-genome-build-and-allele-contract)
1. [Alignment and indexing](#alignment-and-indexing)
1. [Portable storage policy](#portable-storage-policy)
1. [Object scope and mutability](#object-scope-and-mutability)
1. [Optional collection manifests](#optional-collection-manifests)

In separate files:

1. [contigs-and-shards.md](contigs-and-shards.md): canonical contig set,
   row-order contract, shard discovery, and shard subsetting rules.
1. [performance-contract.md](performance-contract.md): cache purpose and
   performance-oriented storage requirements.
1. [python.md](python.md): Python package layout and implementation conventions.
1. [matlab.md](matlab.md): MATLAB/Octave package layout and implementation conventions.
1. [testing.md](testing.md): pytest and Octave consistency strategy.

Objects:

1. [reference.md](reference.md): `ReferenceShard` / `ReferencePanel` —
   SNP tables close to PLINK BIM files.
1. [genotype.md](genotype.md): `GenotypeShard` / `GenotypePanel` —
   PLINK bfile-backed genotype objects.
1. [ld.md](ld.md): `LDShard` / `LDPanel` —
   per-chromosome LD objects and multi-shard panels.
1. [annotations.md](annotations.md): `AnnotationShard` / `AnnotationPanel` —
   SNP annotation matrices or masks, including BED-to-BIM painting.
1. [sumstats.md](sumstats.md): `SumstatsShard` / `Sumstats` —
   GWAS summary statistics for one trait/source.

Each object spec should define what the object represents; its portable disk representation;
its public API and behavioral invariants.

## API notation

API descriptions in the spec are language-agnostic. Signatures use generic pseudocode
notation and define the logical contract — what arguments are required, what
is returned, and what invariants hold. Function-like names such as
`load_reference` are logical API names, not a mandate to implement module-level
free functions. Python and MATLAB/Octave specifics (naming conventions, types,
method vs. function style, native containers, and cache internals) belong in
implementation notes, not in the contract spec.

API notation uses neutral data-shape terms:

- `vector` means a one-dimensional logical column of values; Python may expose
  it as a 1-D `numpy.ndarray`, while MATLAB/Octave should expose a column
  vector unless local convention requires otherwise.
- `matrix` means a two-dimensional numeric or logical array.
- `table` means a column-oriented record collection with named fields/columns.
- `optional X` means the argument may be omitted. Python may represent omitted
  optional values with `None`; MATLAB/Octave may use `[]` or an omitted
  name-value argument.
- Missing optional return fields use the language's natural empty sentinel:
  `None` in Python and `[]` in MATLAB/Octave.

## Contig naming, Genome build and Allele contract

`statgen` starts after variant harmonization. Variant-bearing inputs are
expected to have already passed through [`genomatch`](https://github.com/precimed/genomatch)
or an equivalent upstream pipeline so they share one genome build,
contig naming for chromosome labels, and the `a1=non-reference`, `a2=reference`
allele contract. `statgen` uses these inputs as-is. It does not normalize
chromosome labels, swap alleles, infer reference alleles, liftover coordinates,
or repair strand issues.

Within the current contract, supported contigs are `1`-`22` and `X`.
`statgen` currently does not depend on inputs to be in a specific reference genome,
but it's recommended to provide inputs in `GRCh38`.

Input variant records use `a1` as the non-reference allele and `a2` as the
reference allele (i.e. match reference fasta sequence).
Both are non-empty strings of uppercase DNA bases (`A`, `C`,
`T`, `G`); either may be multi-base for indels. Loaders validate allele syntax
and preserve values as-is without normalization or case conversion. They do not
validate whether either allele matches the reference genome sequence at that
position.

Signed quantities must already be oriented to `a1`:

- GWAS signed effects and Z scores must be oriented to `a1`;
- genotype dosage/count data must count `a1`;
- LD correlations must be signed with respect to `a1` dosages.

## Alignment and indexing

Every per-SNP vector or matrix row is aligned to the same reference-panel SNP
order in memory. Rows are never silently dropped from in-memory aligned objects;
missing or filtered values are represented by masks or `NaN`.
`ReferenceShard`/`ReferencePanel` preserve the loaded order in memory; all
other aligned objects use that reference-defined order.

On disk, source-like objects may store only observed rows plus enough metadata
to project onto the reference panel. Aligned full-panel files are allowed for
caches and cross-language fixtures, but alignment is primarily an in-memory
contract.

A shard is a per-contig partition; a panel is an ordered collection of shards.
Within each shard rows have shard-local SNP indices; across a panel they have
global SNP indices. In memory, all panel-like objects are represented as ordered
shard vectors.

Reference checksum construction is defined in [reference.md](reference.md).
Reference-aligned objects store the paired `ReferenceShard` checksum in memory
and in caches/metadata where applicable.

## Portable storage policy

Portable object files should be language-agnostic. Use:

- PLINK `.bim` for reference shards and PLINK `.bed/.bim/.fam` for genotype
  shards;
- BED for canonical annotation interval inputs;
- `.tsv.gz` for portable non-sharded GWAS summary statistics;
- JSON for manifests and small metadata;
- raw little-endian binary arrays plus JSON metadata for LD triplets.

Avoid pickle or MATLAB `.mat` as portable storage. Language-specific
formats are allowed only as caches for performance or interoperability, and
must be reproducible from portable source files.

## Object scope and mutability

`statgen` objects are **load-only**: they are created exclusively by loader
functions (`load_reference`, `load_ld`, etc.) or cache loaders. There are no
public constructors. User code does not build objects field by field; it loads
from canonical disk sources.

Objects are **immutable after loading**. No setter methods or in-place
modifications are defined. The object holds exactly what was on disk; any
derived quantity (a filtered vector, a reweighted matrix) lives in a plain
array in the caller's workspace, not in the object.

The operations exposed by objects are limited to what is needed for basic
statistical genetics workflows:

- **Underlying data access**: read-only accessors that return plain
  language-native arrays (numpy `ndarray`, MATLAB matrix/vector). For panel
  objects these concatenate shards transparently into genome-wide arrays.
  Examples: `AnnotationPanel.annomat`, `Sumstats.zvec`, `LDPanel.mafvec`.
- **Object-specific operations**: defined in each object spec (for example
  `LDPanel.multiply_r2` and `fast_prune` in [ld.md](ld.md)).

Everything beyond this — regression, likelihood, enrichment, model fitting — is
user-space code that operates on plain arrays extracted from these objects.

## Optional collection manifests

Collections are optional convenience manifests, not the central storage concept.
They can record paths and checksums for a named analysis input set, but object
loaders must work directly on object paths.

If used, `manifest.json` records `schema_version`, reference build, chromosome
list, object file paths, and checksums for the objects included in that
analysis. Cache files are not required to be listed unless they are distributed
intentionally.
