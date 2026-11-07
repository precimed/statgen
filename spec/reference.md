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
  [SPEC.md](SPEC.md); `a1` is the alternative allele (ALT), and `a2` is
  the reference allele (REF).

Optional columns include `source` and `variant_id`.

## Disk representation

Sites-only `.vcf` and PLINK 2 `.pvar` illustrate the reference panel's REF/ALT
semantics. Reference source loading supports PLINK-style `.bim`, with allele
roles specified in
[SPEC.md](SPEC.md#contig-naming-genome-build-and-allele-contract):

```text
chr  snp  cm  bp  a1  a2
```

The six BIM columns may be separated by tabs or by runs of ASCII whitespace.
Fields must not contain embedded whitespace.

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
Each shard carries the MD5 reference checksum defined below. A
`ReferenceShard` contains exactly one chromosome label: all loaded BIM rows in
the shard must have `chr == ReferenceShard.label`. Shard constructors/loaders
must validate this invariant.

## Reference identity keys

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

Each variant also has fixed-width allele identity hashes for
sumstats-to-reference matching, which is the fast object-alignment path that
needs exact variant identity. Chromosome identity is supplied by the shard
label, and base-pair position is already stored as a numeric vector, so only the
variable-length allele strings need compact hashed representations.

`a1_hash64` and `a2_hash64` are `uint64` vectors, length `num_snp`, computed
from allele strings as loaded, using a capped prefix. The cap, warning policy,
cross-runtime implementation requirements, and rationale are specified in
[allele-handling.md](allele-handling.md).
Hash inputs are assumed to be allele strings that already satisfy
the source allele contract in [SPEC.md](SPEC.md). The hash algorithm is a
deterministic two-lane polynomial hash over the UTF-8 bytes of the allele
prefix, chosen so all intermediate
arithmetic is exactly representable in MATLAB/Octave `double` without relying
on unsigned integer overflow:

```text
p = 2147483647  # 2^31 - 1, prime
base1 = 257
base2 = 263
h1 = 1
h2 = 1
for byte in allele_utf8_bytes(hash_prefix(allele)):
    x = double(byte) + 1
    h1 = mod(h1 * base1 + x, p)
    h2 = mod(h2 * base2 + x, p)
allele_hash64 = bitshift(uint64(h1), 32) + uint64(h2)
```

Because `h1 < p`, `h2 < p`, and both bases are small, all products and sums in
the recurrence are below `2^53` and can be computed exactly in IEEE-754 double
precision. Implementations must not use saturating `uint64` multiplication to
compute allele hashes.

The construction also bounds every emitted hash below `2^63`:

```text
max_allele_hash64 = ((p - 1) << 32) + (p - 1)
                  = 9223372030412324862
                  < 9223372036854775808  # 2^63
```

This bound is part of the cross-runtime contract. It allows runtimes without a
native unsigned 64-bit scalar, such as R, to store allele hashes exactly in a
signed 64-bit representation.

## Source-to-reference matching

An exact match requires the same chromosome and position, with the source
effect allele (sumstats) or counted allele (genotypes) matching reference `a1`
(ALT), and the source other or non-counted allele matching reference `a2` (REF).
`load_sumstats` and `load_genotype` implement this using `(shard label, bp,
a1_hash64, a2_hash64)`; marker IDs are ignored. Reference validation follows
[contigs-and-shards.md](contigs-and-shards.md#duplicate-key-validation).
Compatibility checks use the orientation-sensitive `reference_checksum`, not
allele hashes; hash collisions are assumed negligible for non-adversarial inputs.

Source order is not relevant. Match in two passes:

1. For each reference slot, assign its exact match, regardless of
   allele length. More than one exact source candidate is an error.
2. Only for still-unmatched reference slots, try reversed alleles from unused
   source rows. More than one eligible swapped candidate is an error.
   SNV swaps are always enabled. Indel swaps require the per-load logical option
   `allow_indel_swaps=true` (default `false`; MATLAB/Octave:
   `'allow_indel_swaps', true`). Here, an eligible indel has exactly one
   single-base allele and one longer allele. Pairs with two multi-base alleles
   match exactly but never through swap fallback.

Each source row fills at most one reference slot. Exact matches win even when
their optional statistics or calls are missing. Do not attempt swaps for source
rows already assigned by exact matching. Swapped candidates for a reference slot
already filled by an exact match do not create ambiguity; retain the exact
match. Source duplicates that cannot populate a reference slot need no check.
Reciprocal source rows such as `A/G` and `G/A` are allowed; the exact row wins.
Reciprocal indels such as `A/AG` and `AG/A` may also coexist in the reference:
a source row used by its exact reference slot cannot fill the reversed slot.

Swaps exchange complete alleles at the same coordinate and assume consistent
strand, including for `A/T` and `C/G`. They do not complement strands, normalize
alleles, rewrite source rows or reference data, or change checksums. Values are
oriented to reference `a1` as defined in the sumstats and genotype specs.

At verbosity `info`, report only nonzero counts, identified by source file
and shard:

- Successful swaps, with `A/T` and `C/G` swaps counted separately as a subset.
- Indel rows unmatched because indel swaps are disabled: unused rows whose
  reversed key matches a still-unmatched reference slot.

These counts are informational, not warnings, and are suppressed at `quiet`.
Unused reciprocal rows superseded by exact matches are not reported.

Matching, swap eligibility, and `A/T`/`C/G` classification must use allele
hashes only, without decoding or inspecting allele strings. Single-base
alleles are identified by comparing hashes with the hashes of `A`, `C`, `G`,
and `T`, computed with the hash function above. Matching must preserve all
64 bits of allele hashes, use native numeric operations, and avoid composite
string keys.

Regardless of how the panel is sharded, callers access genome-wide vectors and
matrices through read-only panel-level accessors that concatenate across shards
transparently.

Reference panels expose two derived logical variant masks computed from
`a1_hash64` and `a2_hash64`:

- `is_single_nucleotide_variant`: true when both alleles are single
  nucleotides (`A`, `C`, `G`, or `T`);
- `is_strand_ambiguous`: true for unordered `A/T` and `C/G` allele pairs.

## Representation

Implementations may use language-native containers. They must preserve shard
order, expose the required vectors, retain per-shard checksums, and provide
zero-based shard offsets for compatibility with portable metadata.

## Cache layout

### Shared requirements

Newline-delimited reference string payloads encode a shard-local string vector
by joining complete values with `\n`, with no trailing newline. An empty shard
has an empty payload and decodes to zero values. Because BIM fields must not
contain embedded whitespace and alleles are non-empty, newline is a safe
separator. Decoded vector lengths MUST equal the shard row count. This encoding
is independent of the reference checksum byte stream, which retains its
record-ending newlines.

Cache writers MUST NOT pad SNP-ID or allele columns to their longest value.
Schema migration MUST preserve stored hashes and checksums unchanged.

### Python

Python reference caches are `.npz` files. Writers use `reference_cache/0.2`;
readers support `reference_cache/0.1` and `reference_cache/0.2`.

Both Python schemas store `_meta` as a one-dimensional `uint8` array containing
UTF-8 JSON with `schema`, `shard_labels`, and `shard_checksums`. Shard index `i`
follows metadata order. Both schemas store `s{i}_chr` as a fixed-width Unicode
array, `s{i}_bp` as an `int64` vector, and `s{i}_a1_hash64` and `s{i}_a2_hash64`
as `uint64` vectors. Shard row counts come from `s{i}_bp`; numeric, hash, and
chromosome vector lengths MUST agree with them.

#### Schema `reference_cache/0.1` (legacy reader support)

`s{i}_snp`, `s{i}_a1`, and `s{i}_a2` are fixed-width Unicode arrays of
individual strings, with lengths equal to the shard row count. Reading this
schema eagerly produces NumPy object arrays of individual strings.

#### Schema `reference_cache/0.2` (current writer format)

`s{i}_snp_text`, `s{i}_a1_text`, and `s{i}_a2_text` replace the legacy string
arrays. Each is a one-dimensional `uint8` array containing a UTF-8 encoded
newline-delimited payload as defined above. Readers reject invalid UTF-8 and
eagerly decode selected shards into the same object-array representation as
schema `0.1`.

### MATLAB/Octave — schema `reference_cache/0.1`

MATLAB/Octave reference caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
bp
snp_text_by_shard
a1_text_by_shard
a2_text_by_shard
a1_hash64
a2_hash64
```

`chr` is not stored as a cache variable; cache-loaded references synthesize it
from `metadata.shard_labels` and shard row counts.

`metadata` is a struct with:

```text
schema = "reference_cache/0.1"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
```

`bp` is a panel-wide numeric vector. `a1_hash64` and `a2_hash64` are
panel-wide `uint64` vectors. The BIM `cm` field is not cached. Shard offsets
are zero-based half-open intervals into the panel-wide variables and are
sufficient to reconstruct `ReferenceShard` objects.

`snp_text_by_shard`, `a1_text_by_shard`, and `a2_text_by_shard` are cell arrays
with one character vector per shard in panel order. Each character vector
uses the newline-delimited encoding above. Shard `num_snp` is implied by
`shard_start0`/`shard_stop0`.

MATLAB/Octave implementations should keep these per-shard string payloads in
their encoded form after cache load and instantiate cell arrays lazily when
`snp`, `a1`, or `a2` is accessed. `ReferencePanel.chr` is synthesized lazily by
repeating each shard label `num_snp` times and concatenating shard vectors in
panel order; it is not stored in the cache. `bp`, `a1_hash64`, and `a2_hash64`
are available without decoding string payloads.

### R — schema `reference_cache/0.1`

R reference caches are RDS files containing one named list with the same logical
top-level fields and metadata fields as the MATLAB/Octave cache layout above.
R `metadata` is a named list, `a1_hash64` and `a2_hash64` are panel-wide
`bit64::integer64` vectors, and `snp_text_by_shard`, `a1_text_by_shard`, and
`a2_text_by_shard` are character vectors with one newline-delimited shard
payload per shard in panel order. The RDS payload stores these encoded shard
strings, not decoded per-SNP character vectors; cache-loaded R references should
keep them encoded until the corresponding accessor is called.

### Validation and API

The `a1_hash64`/`a2_hash64` cache fields and sumstats matching semantics are
cross-runtime object contracts. Cache metadata validation should be cheap,
depending on shard count and array dimensions rather than scanning all SNP
values. MATLAB/Octave and R should avoid decoding string payloads on the default
load path; Python performs the eager decoding specified above.
`load_reference_cache` trusts stored shard checksums, matching the cache
behavior of reference-aligned objects that cannot recompute those checksums
themselves. Callers who want to verify reference cache integrity may explicitly
call `ReferencePanel.validate_checksums()`, which forces string payload
materialization.

There is one logical reference cache API. Runtime storage details are
language-specific cache concerns governed by
[performance-contract.md](performance-contract.md).

## API

```text
load_reference(path, optional shards) -> ReferencePanel
save_reference_cache(panel, path) -> void
load_reference_cache(path, optional shards) -> ReferencePanel

ReferencePanel.num_snp -> int
ReferencePanel.chr -> num_snp chromosome-label vector
ReferencePanel.snp -> num_snp string vector
ReferencePanel.bp  -> num_snp integer vector
ReferencePanel.a1  -> num_snp string vector
ReferencePanel.a2  -> num_snp string vector
ReferencePanel.a1_hash64 -> num_snp uint64 vector
ReferencePanel.a2_hash64 -> num_snp uint64 vector
ReferencePanel.is_single_nucleotide_variant -> num_snp logical vector
ReferencePanel.is_strand_ambiguous -> num_snp logical vector
ReferencePanel.shard_offsets -> table with shard_label, start0, stop0
ReferencePanel.select_shards(shards) -> ReferencePanel
ReferencePanel.is_object_compatible(object) -> bool
ReferencePanel.validate_checksums() -> bool
ReferencePanel.save_cache(path) -> void
```

Expected behavior:

- Reference panels are external inputs; reference APIs never write `.bim`
  files. LD distribution builders and runtime-native LD converters are the
  only tools that write `.bim` files, and only as exact copies of input BIM rows
  bundled alongside LD shards.
- Cache is a single file (non-sharded). Runtime cache formats follow
  [performance-contract.md](performance-contract.md); the MATLAB/Octave cache
  layout is specified in the "Cache layout" section above.
- Shard discovery, contig validation, row-order validation, and shard subsetting
  follow [contigs-and-shards.md](contigs-and-shards.md).
- Reference shard construction and source loading fail clearly if any variant
  has `a1 == a2`.
- Each `ReferenceShard` must contain exactly one chromosome label. This is a
  global invariant of the object model, not just a cache-loader assumption.
- Compute and retain each shard reference checksum and each per-row
  `a1_hash64` and `a2_hash64` value.
- Accessors are read-only, concatenate shard columns in reference panel order,
  and return plain language-native vectors or tables.
- Cache-loaded references expose `snp`, `a1`, and `a2`; MATLAB/Octave may
  materialize those vectors lazily from the per-shard newline-delimited cache
  payloads. `chr` is synthesized lazily by repeating each shard label
  `num_snp` times and concatenating shard vectors in panel order, without
  caching.
- `shard_offsets` uses zero-based half-open intervals into genome-wide arrays.
- `shard_offsets.start0`/`stop0` are cross-language coordinate metadata, not
  direct language indices. MATLAB/Octave and R callers convert at use-site
  (`start0 + 1 : stop0`).
- `load_reference_cache` skips source-style row validation; `shards` subsetting
  applies against cached shard labels per
  [contigs-and-shards.md](contigs-and-shards.md).
- `ReferencePanel.select_shards(shards)` preserves cache-loaded shard payloads
  without forcing lazy string materialization. Selected Python shards retain
  their existing decoded object arrays; MATLAB/Octave and R shards retain
  encoded payloads until the corresponding string accessors are used.
- `load_reference_cache` validates runtime-specific shard metadata and vector
  lengths before reconstructing shard objects. For MATLAB/Octave and R,
  per-shard string payload containers must match `n_shards`; decoded vector
  lengths may be validated lazily when materialized. Python validates decoded
  lengths eagerly as specified above. Cache loaders trust stored shard
  checksums and do not recompute them.
- `ReferencePanel.validate_checksums()` recomputes each shard checksum from the
  current `chr`, `bp`, `a1`, and `a2` values, compares it to the stored
  `ReferenceShard.checksum`, fails on mismatch with the shard label, and
  returns `true` on success. For cache-loaded references, this explicit
  validation path forces materialization of the per-shard `a1` and `a2` string
  payloads needed to recompute the checksum.
- `ReferencePanel.save_cache(...)` is a thin convenience method equivalent to
  `save_reference_cache(panel, ...)`.
- `is_object_compatible` checks whether a loaded statgen object is aligned to
  this reference panel. Compatibility requires the same ordered shard labels,
  matching shard row counts, and matching shard reference checksums where
  available. Reference-aligned object shards expose `reference_checksum`;
  `ReferenceShard` itself exposes `checksum`.
  Returns `true` when compatible, `false` otherwise. Implementations may log
  informational messages or warnings. Compatibility mismatches do not raise
  exceptions; callers branch on the returned boolean.
