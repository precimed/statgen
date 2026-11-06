# Historical allele hashes and caches

Captured before kernel/validator changes at source commit
`c36657c35d762bd6a92f1175313bac4a51167088` using `generate.py` from the
repository root in the statgen environment. Python, R, and GNU Octave historical
outputs were checked bit-for-bit before capture completed. The MATLAB-format
cache was written by Octave using MAT v7, as in the existing cache fixtures.
Native MATLAB was unavailable.

`vectors.tsv` represents each sequence by a repeated motif truncated to the
specified length. Historical values came from all three production kernels.
New expected values came solely from the generator's scalar implementation of
the two-lane recurrence in spec/reference.md. Hashes are decimal strings to
preserve all bits across readers. All selected vectors above length 150 have
new expected values distinct from their historical values.

`legacy.bim` contains the 193 kb vector paired with T. Each `legacy_*` cache
stores this complete sequence, its original checksum and former-prefix hash.
These are compatibility artifacts, not files to regenerate using new writers.
They are deliberately separate from `make fixtures` and the existing committed
cache fixtures, which remain unchanged.
