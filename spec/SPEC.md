# statgen specification

This directory defines the object-centered design for `statgen`, a shared
Python, MATLAB/Octave, and R package for statistical genetics data structures
and the basic operations on them needed for downstream inference.

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
1. [Validation boundaries](#validation-boundaries)
1. [Object scope and mutability](#object-scope-and-mutability)
1. [Runtime verbosity](#runtime-verbosity)
1. [Optional collection manifests](#optional-collection-manifests)

In separate files:

1. [contigs-and-shards.md](contigs-and-shards.md): canonical contig set,
   row-order contract, shard discovery, and shard subsetting rules.
1. [performance-contract.md](performance-contract.md): cache purpose and
   performance-oriented storage requirements.
1. [allele-handling.md](allele-handling.md): long-allele validation, hash limits,
   and exact hash generation across runtimes, with design rationale.
1. [python.md](python.md): Python package layout and implementation conventions.
1. [matlab.md](matlab.md): MATLAB/Octave package layout and implementation conventions.
1. [R.md](R.md): R package layout, CRAN constraints, and implementation conventions.
1. [docstrings.md](docstrings.md): interactive help text policy for public APIs.
1. [testing.md](testing.md): pytest, Octave, and R consistency strategy.

Objects:

1. [reference.md](reference.md): `ReferenceShard` / `ReferencePanel` —
   variant lists with genomic REF/ALT semantics, as in sites-only `.vcf` or
   PLINK 2 `.pvar`; currently loaded from PLINK 1 `.bim` only.
1. [genotype.md](genotype.md): `GenotypeShard` / `GenotypePanel` —
   reference-aligned PLINK bfile-backed genotype access objects.
1. [ld.md](ld.md): `LDShard` / `LDPanel` —
   per-chromosome LD objects and multi-shard panels.
1. [annotations.md](annotations.md): `AnnotationShard` / `AnnotationPanel` —
   SNP annotation matrices or masks, including BED-to-BIM painting.
1. [sumstats.md](sumstats.md): `SumstatsShard` / `Sumstats` —
   GWAS summary statistics for one trait/source.

Each object spec should define what the object represents, its disk
representation, its public API, and behavioral invariants.

## API notation

API descriptions in the spec are language-agnostic. Signatures use generic pseudocode
notation and define the logical contract — what arguments are required, what
is returned, and what invariants hold. Function-like names such as
`load_reference` are logical API names, not a mandate to implement module-level
free functions. Python, MATLAB/Octave, and R specifics (naming conventions,
types, method vs. function style, native containers, and cache internals)
belong in implementation notes, not in the contract spec.

API notation uses neutral data-shape terms:

- `vector` means a one-dimensional logical column of values; Python may expose
  it as a 1-D `numpy.ndarray`, MATLAB/Octave should expose a column vector
  unless local convention requires otherwise, and R may expose it as an
  ordinary vector.
- `matrix` means a two-dimensional numeric or logical array.
- `table` means a column-oriented record collection with named fields/columns.
- `optional X` means the argument may be omitted. Python may represent omitted
  optional values with `None`; MATLAB/Octave may use `[]` or an omitted
  name-value argument; R may use `NULL` or an omitted argument.
- Missing optional return fields use the language's natural empty sentinel:
  `None` in Python, `[]` in MATLAB/Octave, and `NULL` in R.

## Contig naming, Genome build and Allele contract

`statgen` starts after variant harmonization. Variant-bearing inputs are
expected to have already passed through [`genomatch`](https://github.com/precimed/genomatch)
or an equivalent upstream pipeline so they share one genome build,
contig naming for chromosome labels, and strand. `statgen` does not normalize
chromosome labels, determine or verify REF alleles against a reference-genome
FASTA sequence, liftover coordinates, or repair strand issues. Sumstats and genotype source loaders accept allele
swaps under the matching contract in [reference.md](reference.md).

Within the current contract, supported contigs are `1`-`22` and `X`.
`statgen` currently does not depend on inputs to be in a specific reference genome,
but it's recommended to provide inputs in `GRCh38`.

A `ReferencePanel` is a user-supplied variant table defining genomic REF/ALT
identity and the output allele orientation. Allele roles depend on the input:

For reference input in `.vcf` or `.pvar` format (not yet implemented), `REF`
is expected to match the reference genome sequence, and `ALT` represents the
alternative allele. For reference `.bim` input, column 5 (`a1`) represents ALT,
and column 6 (`a2`) represents REF and is expected to match the reference genome
sequence.

For genotype `.bim` input, column 5 (`a1`) is the counted allele, and column 6
(`a2`) is the non-counted allele. Preferably, the counted allele matches the
reference panel's ALT, and the non-counted allele matches its REF.

For sumstats, `a1` / `EffectAllele` is the effect allele, and
`a2` / `OtherAllele` is the other allele. Preferably, the effect allele matches
the reference panel's ALT, and the other allele matches its REF.

For alleles in genotype or sumstats that are reversed relative to the reference,
loaders orient values using the shared
[matching policy](reference.md#source-to-reference-matching), under the
same-strand assumption. For indels, reversing allele roles can change the
insertion/deletion interpretation. Users who set `allow_indel_swaps=true` on
`load_genotype` or `load_sumstats` take responsibility for interpreting reversed
indel alleles as the same biological variant.

Users establish genome build, strand, and reference alignment upstream.
`statgen` neither requires nor reads a reference FASTA to verify REF identity.
Both alleles are non-empty strings of uppercase DNA bases (`A`, `C`,
`T`, `G`); either may be multi-base for indels. Loaders validate allele syntax
and preserve values as-is without normalization or case conversion.
Alleles on supported contigs must satisfy this syntax contract. Validation
of alleles on ignored contigs such as `Y` and `MT` is implementation-dependent;
loaders need not add special handling for these rows.

Source signed effects, Z scores, and EAF refer to the source effect allele
(`EffectAllele`); source genotype counts count source `a1`. Aligned effects, EAF, and
genotype counts refer to reference ALT (`a1`). LD correlations are signed with
respect to reference ALT dosages.

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
Panel objects must contain at least one shard; empty panels are invalid.

Unless a field is explicitly documented as portable-coordinate metadata (for
example `start0`, `stop0`, or `index_base`), indexing semantics in code follow
the host language (Python 0-based, MATLAB/Octave and R 1-based).

Reference checksum construction is defined in [reference.md](reference.md).
Reference-aligned objects store the paired `ReferenceShard` checksum in memory
and in caches/metadata where applicable.

## Portable storage policy

Portable object files should be language-agnostic. Use:

- PLINK `.bim` for reference shards and PLINK `.bed/.bim/.fam` for genotype
  shards;
- BED and BED-like tab-separated interval files for canonical annotation inputs;
- `.tsv.gz` for portable non-sharded GWAS summary statistics;
- JSON for manifests and small metadata.

Avoid pickle or MATLAB `.mat` as portable storage. Language-specific
formats are allowed only as caches for performance or interoperability, and
must be reproducible from portable source files.

LD panels are the documented exception to this portable-storage policy: their
distributed shard files are runtime-native sparse matrix artifacts, not caches.
Python LD distributions use NumPy/SciPy `.npz`; MATLAB/Octave LD distributions
use `.mat`; R loads Python `.npz` LD shards with an R reference-cache sidecar.
See [ld.md](ld.md) for the LD-specific distribution contract.

Portable source formats define the canonical data contract, not identical
parser leniency across runtimes. For valid canonical inputs, all runtimes must
produce equivalent objects. For non-canonical details such as extra comment
lines, blank-line placement, malformed rows, or implementation-specific table
reader recovery behavior, runtimes may differ as long as they do not silently
change valid data semantics.

## Validation boundaries

`statgen` validates data at the boundary where it enters the library.
External source loaders validate source formats and domain contracts while
source-aware diagnostics are still available. Validated factory functions
validate user-supplied in-memory arrays and vectors before assembling objects.

Cache loaders validate cache schema, required metadata, dimensions, shard
selection, and compatibility gates such as reference checksums where applicable.
They may trust payload values written by `statgen` cache writers and are not
required to repeat full source-style validation on the default load path.

Object constructors are not public user APIs. Internal constructors and
downstream object operations enforce structural invariants needed for
well-formed objects, but should not duplicate semantic validation already
performed by loaders, cache loaders, or validated factories.

Object specs may define stricter boundary validation or explicit QA APIs, such
as checksum recomputation or distribution validation, when those checks are
intentionally outside the default load path.

## Object scope and mutability

`statgen` objects are created by loader/cache-loader APIs or by validated
factory functions defined in object specs. There are no public constructors.
User code does not build objects by mutating fields directly.

Objects are **immutable after creation** (loading or validated factory
construction). No setter methods or in-place
modifications are defined. The object holds exactly what was loaded from disk
or supplied to a validated factory; any derived quantity (a filtered vector, a
reweighted matrix) lives in a plain array in the caller's workspace, not in
the object.

The operations exposed by objects are limited to what is needed for basic
statistical genetics workflows:

- **Underlying data access**: read-only accessors that return plain
  language-native arrays (numpy `ndarray`, MATLAB matrix/vector, R vector or
  matrix). For panel objects these concatenate shards transparently into
  genome-wide arrays.
  Examples: `AnnotationPanel.annomat`, `Sumstats.zvec`, `LDPanel.a1freq`.
- **Object-specific operations**: defined in each object spec (for example
  `LDPanel.multiply_r2` and `fast_prune` in [ld.md](ld.md)).

Everything beyond this — regression, likelihood, enrichment, model fitting — is
user-space code that operates on plain arrays extracted from these objects.

## Runtime verbosity

`statgen` defines a global runtime verbosity setting with levels `quiet` and
`info`; the default is `info`. The public API is `set_verbosity(level)` and
`get_verbosity()`.

`load_ld` consults this setting: at `info` it reports progress while loading
individual shards, and at `quiet` it suppresses that reporting. Sumstats and
genotype source loaders also use verbosity for reporting swap counts as
specified in [reference.md](reference.md#source-to-reference-matching).

## Optional collection manifests

Collections are optional convenience manifests, not the central storage concept.
They can record paths and checksums for a named analysis input set, but object
loaders must work directly on object paths.

If used, `manifest.json` records `schema_version`, reference build, chromosome
list, object file paths, and checksums for the objects included in that
analysis. Cache files are not required to be listed unless they are distributed
intentionally.
