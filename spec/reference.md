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
from allele strings as loaded. If any allele exceeds 150 characters,
implementations must warn and compute the hash from only the first 150
characters. Hash inputs are assumed to be allele strings that already satisfy
the source allele contract in [SPEC.md](SPEC.md). The hash algorithm is a
deterministic two-lane polynomial hash over the UTF-8 bytes of the allele
string prefix, chosen so all intermediate
arithmetic is exactly representable in MATLAB/Octave `double` without relying
on unsigned integer overflow:

```text
p = 2147483647  # 2^31 - 1, prime
base1 = 257
base2 = 263
h1 = 1
h2 = 1
for byte in allele_utf8_bytes(first_150_characters(allele)):
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

Source-to-reference matching for variant-bearing source files, including
summary statistics and genotype BIM sources, uses the tuple `(shard label, bp,
a1_hash64, a2_hash64)`. The allele hashes are not used for reference
compatibility checks; compatibility continues to use the stored per-shard
`reference_checksum`. The allele hashes are not security primitives. Collision
risk is negligible for non-adversarial allele strings in this use.
Duplicate-key validation follows
[contigs-and-shards.md](contigs-and-shards.md#duplicate-key-validation).

After exact matching, implementations warn if any unmatched source-side variant
would have matched the reference within the same shard had `a1_hash64` and
`a2_hash64` been swapped. The warning names the source file, shard label, and
count of such variants. This is diagnostic only: those variants remain
unmatched, and callers must not silently flip alleles, repair strand issues, or
otherwise rescue them.

Implementations should match within each shard using native hash-preserving
operations over `(bp, a1_hash64, a2_hash64)`. MATLAB/Octave implementations
should use a native numeric sort/merge; equivalent native numeric approaches
are allowed. Implementations must not reconstruct string join keys and must
not cast `uint64` hashes to `double` for matching because values above `2^53`
would lose precision.

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
encodes the corresponding shard-local string vector as `strjoin(values, "\n")`
with no trailing newline. This payload encoding is independent of the reference
checksum byte stream, which continues to serialize each checksum record with a
trailing newline as defined above. Because BIM fields must not contain embedded
whitespace and alleles are non-empty, newline is not a legal field value and is
safe as the cache separator. Decoders must reject decoded string-vector lengths
that do not equal the shard `num_snp` implied by `shard_start0`/`shard_stop0`.

MATLAB/Octave implementations should keep these per-shard string payloads in
their encoded form after cache load and instantiate cell arrays lazily when
`snp`, `a1`, or `a2` is accessed. `ReferencePanel.chr` is synthesized lazily by
repeating each shard label `num_snp` times and concatenating shard vectors in
panel order; it is not stored in the cache. `bp`, `a1_hash64`, and `a2_hash64`
are available without decoding string payloads.

R reference caches are RDS files containing one named list with the same logical
top-level fields and metadata fields as the MATLAB/Octave cache layout above.
R `metadata` is a named list, `a1_hash64` and `a2_hash64` are panel-wide
`bit64::integer64` vectors, and `snp_text_by_shard`, `a1_text_by_shard`, and
`a2_text_by_shard` are character vectors with one newline-delimited shard
payload per shard in panel order. The RDS payload stores these encoded shard
strings, not decoded per-SNP character vectors; cache-loaded R references should
keep them encoded until the corresponding accessor is called.

The `a1_hash64`/`a2_hash64` cache fields and sumstats matching semantics are
cross-runtime object contracts. Cache metadata validation should be cheap,
depending on shard count and array dimensions rather than scanning all SNP
values or decoding all string payloads on the default load path.
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
  without forcing string materialization. When a cache-loaded panel is subset in
  memory, selected shards retain their encoded per-shard `snp`, `a1`, and `a2`
  payloads until the corresponding string accessors are used.
- `load_reference_cache` validates shard-offset metadata and panel-wide vector
  lengths cheaply before reconstructing shard objects. It validates that
  per-shard string payload containers match `n_shards`; decoded string-vector
  lengths may be validated lazily when those payloads are materialized. Cache
  loaders trust stored shard checksums and do not recompute them.
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
