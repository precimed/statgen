# Long allele handling

This document defines shared validation and hash-generation requirements for
Python, MATLAB/Octave, and R. Allele syntax and orientation remain in
[SPEC.md](SPEC.md#contig-naming-genome-build-and-allele-contract); the hash
is defined in [reference.md](reference.md#reference-identity-keys).

## Hash input cap and identity

Hash generation MUST use at most the first `1000000` bases of each allele.
`hash_prefix(allele)` denotes this prefix, or the complete allele when shorter.
The cap MUST be an internal constant, without a public parameter or runtime
override. Under the DNA syntax contract, base, character, and byte counts coincide.

For rows subject to allele validation, source loaders MUST inspect lengths
once the native reader has materialized the columns, before padded conversion
or hashing. `load_reference`, `load_genotype`, and `load_sumstats` MUST warn
when the cap is exceeded and continue using the prefix for hashing. Emit at
most one cap warning per source file per load, covering both allele columns.
Warnings SHOULD identify the affected allele count, cap, and first affected
row, column, and observed length, without printing complete sequences.
Warnings MUST originate in loaders, where source locations are known;
`allele_hash64` helpers MUST cap silently, including on repeated internal calls.

Stored alleles, checksums, and syntax validation MUST use complete sequences.
Length alone MUST NOT cause load failure or variant removal. Empty alleles or
invalid DNA anywhere, including beyond the prefix, remain errors. Internal
operations may trust validated inputs. Native readers may materialize complete
strings; bounded streaming parsing is not required. Alleles MUST NOT be
normalized or case-converted.

Alleles up to the cap have full-sequence identity hashes. Above it, identical
prefixes yield identical hashes despite different suffixes; matching and
hash-key duplicate detection cannot distinguish those alleles. The warning
makes this intentional loss of suffix identity visible.

## Performance requirements

For `N` alleles with `S` total sequence bytes and maximum length `Lmax`,
validation/hash kernels MUST use `O(N + S)` work and MUST NOT allocate a global
`N`-by-`Lmax` representation or repeatedly scan all rows for each byte position.
Bounded length-bucket dispatch may additionally scan length vectors. Hashing
work MUST follow capped-prefix bytes; extract prefixes before byte expansion.
Reading and validation still process complete sequences. Sorting and other
object operations are outside this kernel cost requirement.

Retain a vectorized single-base fast path. Every temporary sequence matrix or
byte expansion, including on short-input paths, MUST use bounded batches that
account for numeric expansions and masks. Output/state vectors may scale with
`N`; power tables may scale with the longest hashed prefix, bounded by the cap.
Loops over batches, power-table positions, or rare long alleles are permitted
as an exception to [performance-contract.md](performance-contract.md)'s SNP-axis
loop restrictions; ordinary processing must remain native/vectorized.

## Design per runtime

- MATLAB/Octave SHOULD use vectorized padded validation within bounded batches
  grouped by length, keeping padding close to actual sequence volume, and
  vectorized reductions for hashing. Any retained short-allele padded kernel
  MUST also run in bounded batches.
- Python and R SHOULD use native variable-length string validation and grouped
  reductions over raw-byte batches for hashing, without padded matrices.

## Hash generation

Equivalent evaluations of the two polynomial lanes MUST produce identical bits.
For prefix length `L`, lane base `b`, modulus `p`, initial state `1`, and
`x_i = byte_i + 1`, the recurrence equals:

```text
h = (b^L + sum(x_i * b^(L-i), i=1..L)) mod p
```

Powers MUST be reduced modulo `p` before multiplication. Accumulation MUST
remain exact, using bounded reductions where necessary to avoid floating-point
rounding or integer overflow. Combining the lanes retains the existing integer
representation and `< 2^63` bound.

## Rationale

The larger cap preserves distinguishing sequence in the intended structural-
variant range while bounding hash work for potentially malformed multi-megabase
codes. Warning instead of failing retains historical behavior for exceptional
lengths, deliberately sacrificing suffix identity above the cap. Increasing the
cap alone cannot make global padding or repeated whole-vector scans scale;
bounded processing is required independently. Grouped padding preserves
MATLAB's efficient vectorized validation, while variable-length processing
suits Python and R.

## Verification

All runtimes MUST test shared hash vectors, differences within the hash prefix,
structural variants, and cap boundaries. Cover warning ownership and
deduplication, above-cap prefix equivalence, complete-string preservation,
invalid characters beyond the cap, and bounded over-cap processing.
Performance checks SHOULD include ordinary inputs and rare long outliers,
assess scaling and allocation shapes, and avoid machine-specific timing
thresholds.
