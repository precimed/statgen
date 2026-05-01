# statgen specification

This directory defines the object-centered design for `statgen`, a shared
Python and MATLAB/Octave package for statistical genetics data structures and
the basic operations on them needed for downstream inference.

The design is centered on individual objects. A project may have multiple
genotype panels, LD panels, GWAS summary-statistic objects, annotation panels,
and reference panels in the same workspace. Object loaders must not require all
data types to live under a single bundle root.

`statgen` starts after variant harmonization. Variant-bearing inputs are
expected to have already passed through `genomatch` or an equivalent upstream
pipeline so they share one genome build, `genomatch` NCBI contig naming unless
explicitly declared otherwise, and the `a1=non-reference`, `a2=reference`
allele contract. `statgen` uses these inputs as-is. It does not normalize
chromosome labels, swap alleles, infer reference alleles, liftover coordinates,
or repair strand issues.
Within the current default contract, supported contigs are `1`-`22` and `X`.

Portable disk representation choices:

- `ReferencePanel`, `GenotypePanel`, and `LDPanel` may be sharded or
  non-sharded on disk. Sharded paths use `@` as the shard-label placeholder.
  Non-sharded is an input layout only; in-memory panels are contig-sharded.
  The sharding of an `LDPanel` must match the sharding of the reference used
  with it.
- `Sumstats` is non-sharded on disk: one `.tsv.gz` per trait/source.
- annotations are canonical BED interval inputs, non-sharded on disk;
  painted annotation objects are in-memory representations (sharded).

## Objects

- `ReferenceShard` / `ReferencePanel`: SNP tables close to PLINK BIM files.
- `GenotypeShard` / `GenotypePanel`: PLINK bfile-backed genotype objects.
- `LDShard` / `LDPanel`: per-chromosome LD objects and multi-shard panels.
- `AnnotationShard` / `AnnotationPanel`: SNP annotation matrices or masks.
- `SumstatsShard` / `Sumstats`: GWAS summary statistics for one trait/source;
  `Sumstats` is an ordered collection of `SumstatsShard`s, not a multi-trait
  panel.

Each object spec should define:

1. what the object represents;
2. its portable disk representation;
3. its public API and behavioral invariants.

API descriptions are language-agnostic. Signatures use generic pseudocode
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
- **LD r² matrix multiply**: `LDPanel.multiply_r2(M)` — the only numerical
  primitive on LD objects, covering LD scores, annotation LD
  weighting, and per-SNP variance accumulation.
- **Pruning**: `fast_prune(logpvec, ld_panel, r2_threshold)` — greedy
  significance-based LD pruning, the other primitive consumers of LD structure.

Everything beyond this — regression, likelihood, enrichment, model fitting — is
user-space code that operates on plain arrays extracted from these objects.

## Object provenance

Objects differ in whether `statgen` creates their disk representation or only
consumes it:

- **Genotype panel**: created by external QC and imputation pipelines, not by
  `statgen`. The bfile triplets are consumed as-is and serve as a de-facto
  reference for sample and variant identity.
- **Reference panel**: created externally (typically derived from genotype BIM
  files or a curated SNP list), not by `statgen`. The `.bim` files are consumed
  as-is by the reference loader.
- **LD panel**: created by `statgen` tooling. The typical sequence is:
  1. run PLINK to compute pairwise LD for a chromosome;
  2. run a `statgen` conversion script to transform PLINK output into the
     portable binary triplet format (`ld_idx1.i32`, `ld_idx2.i32`, `ld_r.f32`,
     `mafvec.f32`, `metadata.json`).
- **Summary statistics**: created by external GWAS pipelines, not by `statgen`.
  The `.tsv.gz` files are consumed by the sumstats loader.
- **Annotations**: BED interval files are created and curated externally.
  Painting BED intervals onto a reference (producing `annomat`) is part of the
  `statgen` API and is supported in both Python and MATLAB/Octave.

Loading objects and projecting them onto a reference is part of the `statgen`
API.

## Reading order

1. [conventions.md](conventions.md): scope, input contract, shard/panel
   conventions, index rules, and storage policy.
2. [contigs-and-shards.md](contigs-and-shards.md): canonical contig set,
   row-order contract, shard discovery, and shard subsetting rules.
3. [performance-contract.md](performance-contract.md): cache purpose and
   performance-oriented storage requirements.
4. [reference.md](reference.md): `ReferenceShard` and `ReferencePanel`.
5. [genotype.md](genotype.md): `GenotypeShard` and `GenotypePanel`.
6. [ld.md](ld.md): `LDShard` and `LDPanel`.
7. [annotations.md](annotations.md): `AnnotationShard`,
   `AnnotationPanel`, and BED-to-BIM painting.
8. [sumstats.md](sumstats.md): `SumstatsShard` and `Sumstats`.
9. [testing.md](testing.md): pytest and Octave consistency strategy.
10. [python.md](python.md): Python package layout and implementation conventions.
11. [matlab.md](matlab.md): MATLAB/Octave package layout and implementation conventions.

## Optional collection manifests

Collections are optional convenience manifests, not the central storage concept.
They can record paths and checksums for a named analysis input set, but object
loaders must work directly on object paths.

If used, `manifest.json` records `schema_version`, reference build, chromosome
list, object file paths, and checksums for the objects included in that
analysis. Cache files are not required to be listed unless they are distributed
intentionally.
