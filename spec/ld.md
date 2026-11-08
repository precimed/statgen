# LDShard and LDPanel

## Represents

LD describes the pairwise signed correlation structure between variants in the
same processed shard. Values are signed correlation `r`, not `r²`; consumers
that need `r²` use runtime operations that square LD values element-wise after
loading. Cross-chromosome LD is not stored.
Signs are oriented to `a1` dosages, consistent with the allele contract in
[SPEC.md](SPEC.md).

LD is paired with a reference: matrix row/column order, SNP count, allele
frequencies, and `reference_checksum` correspond to a specific
`ReferenceShard`.

LD panels are a deliberate exception to the package-wide portable-source plus
optional-cache pattern. LD distribution artifacts are runtime-native sparse
matrix files because the expensive sparse-matrix construction is part of
reference-panel preparation, not routine user loading.

## Logical LD shard contract

Each LD shard represents:

- `chr`: shard label, one of `1`-`22` or `X`;
- `sex`: `null` for autosomes, one of `"female"`, `"male"`, or `"combined"`
  for chrX;
- `num_snp`: number of SNPs in the paired reference shard;
- `ld_r`: symmetric sparse signed-correlation matrix, shape
  `num_snp x num_snp`;
- diagonal entries of `ld_r`: explicit `1.0`;
- off-diagonal entries: signed `r` oriented to `a1` dosage;
- missing off-diagonal entries: interpreted as zero;
- `a1freq`: float vector of `a1` allele frequencies aligned to the paired
  `ReferenceShard`;
- `reference_checksum`: MD5 checksum of the paired reference shard;
- build metadata: PLINK command, PLINK version if available, window size,
  `r²` storage threshold, sample count, and schema version.

The paired `ReferencePanel` and `LDPanel` must originate from the same BIM row
order and allele orientation. The LD build additionally depends on the paired
BED/FAM genotype source used to estimate LD and allele frequencies.

Off-diagonal LD `r` values are Pearson correlations between two `a1` dosage
vectors. Pearson correlation is statistically undefined when either vector has
zero variance. For a full-shard genotype vector this includes SNPs with no
observed variation among the samples PLINK2 uses for that shard (`a1freq`
equal to `0` or `1`), including chrX sex-specific shards. Pairwise
missing-call handling can also make LD undefined: two SNPs may each vary in
the shard overall, but after restricting the estimate to samples with
non-missing calls for both SNPs, one or both pairwise dosage vectors may have
zero variance.

LD shards must preserve the paired reference shard's full SNP axis. Builders
must not drop SNPs just because their LD is undefined.

Undefined off-diagonal LD values are represented by omission from the sparse
matrix and therefore read as zero under the sparse-missing convention. This is
a storage and downstream-computation convention, not a claim that the
statistical correlation estimate is zero. `ld_r` must not contain `NaN` or
infinite values. The diagonal is always stored as explicit `1.0`, even for
monomorphic or otherwise pairwise-undefined SNPs, so aligned matrix operations
retain the full BIM/reference shape.

Because marginally monomorphic SNPs have undefined LD with every other SNP,
`statgen_build_ld.py` treats them as a build-time QC failure by default. The
builder checks for `a1freq == 0` or `a1freq == 1` immediately after frequency
calculation and before the heavier PLINK2 LD computation. Users may force
output with `--allow-monomorphic-snps`; forced shards retain the full SNP axis
and use the sparse omission convention above for undefined off-diagonal LD.

Pairwise zero variance caused only by missing-call intersection is documented
but is not a default build-time failure. In ordinary sparse PLINK2 `.vcor`
output with a nonnegative `--ld-window-r2` threshold, such undefined pairs are
omitted and are not distinguishable from below-threshold omitted pairs in the
resulting sparse artifact. Users who need stronger guarantees should apply
upstream missingness QC or run separate diagnostics before building LD.

## Distribution formats

`load_ld` reads the distribution artifact native to the runtime:

- Python reads NumPy/SciPy `.npz` files with CSC sparse components.
- MATLAB/Octave reads `.mat` files containing native sparse matrices.
- R reads the Python `.npz` files directly and constructs `Matrix::dgCMatrix`
  matrices in memory; R package layout and CRAN constraints are specified in
  [R.md](R.md).

The Python `.npz` format is also the documented build handoff format from
`statgen_build_ld.py` into runtime-native converters for MATLAB/Octave.
MATLAB/Octave user loading reads `.mat`, not `.npz`. R user loading reads the
Python `.npz` LD shard files and uses an R-native reference-cache sidecar
recorded in the Python `.npz` manifest. Runtime
validation follows this boundary: Python validates Python `.npz`
distributions, MATLAB/Octave validates MATLAB/Octave `.mat` distributions, and
R validates Python `.npz` distributions that carry the required R
reference-cache sidecar. A runtime validator must otherwise reject another
runtime's manifest format rather than acting as the authority for it.

### Python `.npz` shard format

Each `.npz` file represents one LD shard. Required arrays:

```text
data       float32, length nnz
indices    int32, length nnz, zero-based row indices
indptr     int32, length num_snp + 1, zero-based CSC column pointers
shape      int64, length 2, equals [num_snp, num_snp]
a1freq     float32, length num_snp
metadata   uint8, UTF-8 JSON bytes
```

`data`, `indices`, and `indptr` encode a symmetric sparse signed-`r` matrix in
CSC layout with explicit unit diagonal. Sparse indices are zero-based because
the `.npz` format is Python/SciPy-oriented. The archive must not contain object
arrays or pickled payloads.

CSC32 shards must have `nnz < 2^31`. Builders must reject larger shards rather
than overflowing sparse indices; a future CSC64 schema would be a separate
format/version.

The metadata JSON must include:

```json
{
  "object_type": "ld_shard",
  "schema_version": "1.0",
  "format": "statgen_ld_npz_csc32",
  "chr": "1",
  "sex": null,
  "num_snp": 123,
  "nnz": 4567,
  "sparse_layout": "csc",
  "index_base": 0,
  "matrix": "symmetric",
  "diagonal": "explicit_unit",
  "value": "r",
  "num_monomorphic_snps": 0,
  "reference_checksum": "...",
  "reference_bim": "reference_chr1.bim",
  "build_tool": "plink2",
  "build_command": "plink2 ...",
  "plink_version": "...",
  "ld_window_kb": 10000,
  "ld_r2_threshold": 0.05,
  "num_sample": 1000
}
```

`plink_version` may be `null` if unavailable.

### MATLAB/Octave `.mat` shard format

Each `.mat` file represents one LD shard. Required variables:

```text
ld_r       sparse double, num_snp x num_snp
a1freq     numeric vector, length num_snp
metadata   struct
```

Published production `.mat` LD distribution files must be MAT-file version 7.3
(HDF5-based). The legacy v5 MAT-file format is not acceptable for production
distribution because realistic LD shards can exceed its variable-size limits.
Loaders may accept v5 MAT-files for small synthetic fixtures and local tests
when they contain the required sparse `ld_r`, `a1freq`, and `metadata`
variables.

`ld_r` must be native sparse, symmetric, signed `r`, explicitly unit-diagonal,
and aligned to the paired reference shard. The loader must load `ld_r`
directly; it must not reconstruct the sparse matrix from COO or CSC arrays on
the default user path.

`metadata` must contain the same logical fields as the `.npz` metadata, with
MATLAB-valid struct field names and `format: "statgen_ld_mat_sparse_double"`.
Storage-layout fields that are specific to `.npz` CSC arrays, such as
`sparse_layout` and `index_base`, are not required in `.mat` metadata.
Unavailable optional metadata values, such as an unknown `plink_version`, may
use the runtime's natural empty value.

## Panel layout and manifest

A panel root is a directory containing one or more LD shard files, one bundled
reference `.bim` file per chromosome, runtime reference cache files, and
`ld_manifest.json`. `path` passed to `load_ld` must always identify a panel
root directory; single shard files are not a supported load target.

Recommended filenames:

```text
reference_chr1.bim
reference_chr2.bim
...
reference_chr22.bim
reference_chrX.bim
reference_cache.npz
reference_cache.rds
ld_chr1.npz
ld_chr2.npz
...
ld_chr22.npz
ld_chrX_female.npz
ld_chrX_male.npz
ld_chrX_combined.npz
ld_chr1.npz.d/
```

MATLAB/Octave distributions use the same names with `.mat` for LD shard files
and `reference_cache.mat` for the reference cache. Python `.npz` distributions
that are intended to be loaded by R keep the Python `reference_cache.npz` and
additionally include an R `reference_cache.rds` sidecar named by
`r_reference_cache`. Bundled `.bim` files are identical across runtimes.
Autosomal LD shards use no sex suffix. chrX LD shards always use an explicit
sex suffix. Not all chrX sex labels need to be present. `female` and `male`
are the default sex-specific build outputs; `combined` is optional and
assumption-dependent. Each chromosome has exactly one bundled reference `.bim`
regardless of how many chrX sex labels are present.

`ld_manifest.json` is authoritative for file discovery. Loaders must not infer
panels by globbing arbitrary filenames. Files present in the directory but
absent from the manifest are ignored.

For R, a manifest-declared `.npz` shard may have an optional sibling extracted
cache directory named by appending `.d` to the full shard filename, for example
`ld_chr1.npz.d/`. The extracted directory is a cache of the corresponding
`.npz` file, not a manifest target and not an independent LD shard format.
When present and current, R may read the extracted `.npy` members directly;
otherwise it must fall back to the manifest-declared `.npz` file.

Required manifest fields for Python `.npz` distributions:

```json
{
  "object_type": "ld_panel_manifest",
  "schema_version": "1.0",
  "runtime_format": "python_npz_csc32",
  "reference_cache": "reference_cache.npz",
  "reference_cache_md5": "...",
  "shards": [
    {
      "chr": "1",
      "sex": null,
      "file": "ld_chr1.npz",
      "file_md5": "...",
      "num_snp": 123,
      "nnz": 4567,
      "reference_checksum": "...",
      "reference_bim": "reference_chr1.bim"
    },
    {
      "chr": "X",
      "sex": "female",
      "file": "ld_chrX_female.npz",
      "file_md5": "...",
      "num_snp": 456,
      "nnz": 7890,
      "reference_checksum": "...",
      "reference_bim": "reference_chrX.bim"
    },
    {
      "chr": "X",
      "sex": "male",
      "file": "ld_chrX_male.npz",
      "file_md5": "...",
      "num_snp": 456,
      "nnz": 3456,
      "reference_checksum": "...",
      "reference_bim": "reference_chrX.bim"
    }
  ]
}
```

Additional manifest fields required for R-loadable Python `.npz`
distributions:

```json
{
  "r_reference_cache": "reference_cache.rds",
  "r_reference_cache_md5": "..."
}
```

`runtime_format` is `"python_npz_csc32"` for Python/R `.npz` distributions and
`"matlab_mat_sparse_double"` for MATLAB/Octave distributions.
`runtime_format` is a manifest-level distribution identifier: it lets a loader
reject a panel root that belongs to another runtime before opening shard files.
The per-shard metadata `format` field is a shard payload identifier and must
agree with the runtime distribution format. Validators use both fields so a
manifest cannot silently point to shard files with the wrong internal schema.
`reference_cache` is a plain relative filename with no path separators, no
`..` components, and no absolute path prefix. For Python/R `.npz`
distributions, it names the Python reference cache and must not be repurposed
for R. `reference_cache_md5` is manifest-only. Python/R `.npz` distributions
that are freshly produced by the Python builder omit `r_reference_cache` and
`r_reference_cache_md5` until R preparation is run. Python/R `.npz`
distributions that are intended to be loaded by R must additionally contain
those fields. `r_reference_cache` is a plain relative filename naming an R
`ReferencePanel` cache in the same directory as `ld_manifest.json`, and
`r_reference_cache_md5` is its manifest-only MD5. R `load_ld` must reject a
Python `.npz` LD distribution that lacks these R-specific reference-cache
fields rather than rebuilding the reference from bundled BIM files on the
default load path. For MATLAB/Octave distributions, `reference_cache` names
the MATLAB/Octave reference cache.
Manifest MD5 fields (`reference_cache_md5` and per-shard `file_md5`) must be
32-character lowercase hexadecimal strings. Producers must write lowercase
hex, and validators should reject non-lowercase variants rather than accepting
case-insensitive equivalents.
Runtime loaders are not required to compute `file_md5`, `reference_cache_md5`,
or `r_reference_cache_md5` on the default load path.

Manifest entries and per-file metadata must agree on `chr`, `sex`, `num_snp`,
`nnz`, `reference_checksum`, `reference_bim`, and runtime/storage format.
`file` and `file_md5` are manifest-only fields.

`reference_bim` must be a plain relative filename with no path separators, no
`..` components, and no absolute path prefix (e.g. `"reference_chr1.bim"`, not
`"../foo.bim"` or `"/abs/path.bim"`). Loaders must reject manifests that
violate this constraint. chrX sex variants (`female`, `male`, `combined`) all
carry the same `reference_bim` value because they share one paired reference
shard.

`load_ld` loads the manifest-declared reference cache and only the LD files
required for the requested shards. It must not parse bundled reference `.bim`
files on the default user load path. LD files and reference `.bim` files for
other shards may be present in the directory and are ignored.
`load_ld` must resolve the reference cache filename from the manifest rather
than hard-coded filenames or directory globbing. Python uses `reference_cache`
for `.npz` distributions, R uses `r_reference_cache` for `.npz` distributions,
and MATLAB/Octave uses `reference_cache` for `.mat` distributions.

## Building and conversion

LD construction and runtime conversion follow a hub-and-spoke pattern:
Python builds the `.npz` handoff distribution, and MATLAB/Octave converts from
that Python distribution into its runtime-native `.mat` format. R loads the
Python `.npz` distribution directly after the R reference-cache sidecar has
been prepared. Direct R-to-MATLAB, MATLAB-to-R, or other runtime-to-runtime LD
exchange is not a goal.

### Python LD builder

The builder constructs its reference from the same PLINK 1 bfile used to
compute LD, through `load_reference` and its
[reference validation](reference.md#reference-duplicate-key-validation).
PLINK output variants must match exactly one reference row by ordered
`(chr, bp, a1, a2)` strings, without swap fallback. The builder must fail cleanly
on any unmatched variant or ambiguity when matching LD correlations to the
reference panel.

`statgen_build_ld.py` builds the Python `.npz` LD distribution and handoff
files from a PLINK bfile. It accepts both input layouts:

- **Sharded bfile input** (`@` in path, e.g. `chr@`): one bfile per
  chromosome; the requested shard is resolved by canonical substitution.
- **Non-sharded bfile input** (single path, no `@`): the requested shard is
  selected from the single BIM.

`statgen_build_ld.py` requires a single `--shard` label (`1`-`22` or `X`) and
builds only that shard. This applies to both sharded and non-sharded bfile
inputs, so chromosomes can be built as independent parallel jobs.

For the requested shard, the builder:

1. Computes `a1` allele frequencies using PLINK2 `--freq`, aligned to BIM row
   order.
2. Fails before LD computation if any aligned frequency is exactly `0` or `1`,
   unless `--allow-monomorphic-snps` is supplied.
3. Computes pairwise signed LD using PLINK2 `--r-unphased` and `--keep-allele-order`
   with default window 10,000 kb and default `r²` storage threshold `0.05`.
4. Writes one validated `.npz` shard file for autosomes, or one validated
   `.npz` file per chrX sex label.
5. Copies the relevant rows of the source BIM into a bundled
   `reference_chr<N>.bim` file in the output directory. For chrX, writes
   `reference_chrX.bim` once, shared by all sex-label shards.

Each shard metadata records `num_monomorphic_snps`, the number of aligned SNPs
with `a1freq` exactly `0` or `1` for that shard. This field is present even
when the count is zero. When monomorphic SNPs are present and output is forced,
the builder also writes a sidecar QC table next to the shard file, named like
`ld_chr1.monomorphic.tsv` or `ld_chrX_male.monomorphic.tsv`. The sidecar is
created only when monomorphic SNPs are present and has columns:

```text
chr  snp  bp  a1  a2
```

PLINK2 tabular `.vcor` output includes only variant pairs that pass its LD
report filters. With a nonnegative `--ld-window-r2` threshold, undefined
correlations such as `nan` do not pass the filter; with a negative threshold,
PLINK2 can emit `nan` rows. `statgen_build_ld.py` therefore requires a finite,
nonnegative `ld_r2_threshold`, rejects non-finite LD values if they are
encountered in PLINK output, and treats omitted off-diagonal pairs as absent
sparse entries. Omitted pairs include low-`r²` pairs and undefined pairs; both
read as zero from the resulting sparse LD shard.

After shard jobs finish, `statgen_create_ld_manifest.py --ld <root>` creates
`ld_manifest.json` from the per-shard metadata in existing `.npz` files,
requires each shard's metadata to name its bundled `reference_bim`, validates
the resulting panel, builds a `ReferencePanel` from the bundled reference BIM
files, saves that reference with `ReferencePanel.save_cache(...)` under the LD
root, and records the cache filename and MD5 in the manifest. The command must
fail without modifying files if `ld_manifest.json` already exists, because
runtime preparation steps may add fields to that manifest.

Before publishing an `.npz` shard, the builder must validate at least:

- all required array names are present;
- `shape == [num_snp, num_snp]`;
- `len(indptr) == num_snp + 1`;
- `len(data) == len(indices) == indptr[-1]`;
- `len(a1freq) == num_snp`;
- `num_monomorphic_snps` is a non-negative integer no greater than `num_snp`;
- sparse index arrays use int32 and `nnz < 2^31`;
- sparse index bounds are valid;
- metadata schema version is supported;
- metadata and array dimensions agree.

The r² threshold, window size, PLINK2 command, PLINK version if available,
sample count, and reference checksum are build metadata and must be recorded.
The sharding of the resulting LD panel must match the sharding of the reference
used with it.

`statgen_build_ld.py` is not a general PLINK passthrough. Supported PLINK-like
controls are limited to resource controls corresponding to PLINK `--threads`
and `--memory`. For chrX sex-specific builds, statgen generates internal keep
files for female and male shards; these generated files are not user-facing
PLINK passthrough controls. Sample/family filtering, sample-level missingness
filtering, variant-level QC/filtering, and allele/reference mutation flags must
not be passed through; users who need filtering must create a filtered
bfile/reference upstream.
`num_sample` records the sample count submitted to PLINK for each shard after
statgen-owned chrX sex splitting.

By default, chrX builds produce `female` and `male` shards using FAM column 5
(1 = male, 2 = female). The script requires non-missing sex for chrX
sex-specific builds. `combined` chrX output is not the default; it is an
explicit user decision because it depends on modeling and encoding
assumptions. If users suppress sex splitting and write only a single combined
chrX shard, the manifest records `sex: "combined"`.

### MATLAB/Octave LD conversion

The MATLAB/Octave converter reads `.npz` shard files and writes `.mat` shard
files for one requested reference shard at a time. The shard argument is
required so conversion can run as independent parallel jobs without manifest
write races. For chrX, one requested shard label `X` converts all chrX
sex-label shard files present in the Python manifest (`female`, `male`, and/or
`combined`). `convert_ld_npz_to_mat(npz_root, mat_root, shard)` writes
production v7.3 MAT-files by default. Supported `format` values are `v7.3`
(the production default) and `v5`. A caller may explicitly request local v5
MAT-file output with `format = "v5"` in language-specific syntax (for example
MATLAB name-value arguments). v5 output is accepted for fixture-scale tests and
local validation, but is not a production distribution artifact. The v7
MAT-file format is intentionally not supported for LD conversion because it has
the same production-size limitation as v5: realistic LD sparse matrices can
exceed legacy MAT-file per-variable size limits.

After all requested shard conversions finish,
`create_ld_mat_manifest(npz_root, mat_root, shards)` creates the
MATLAB/Octave `ld_manifest.json`. The finalizer reads the Python
`ld_manifest.json` only to determine which shard files are expected for the
requested reference shard labels. For each requested label, absence from the
Python manifest is an error. For every expected MATLAB/Octave shard file, the
finalizer must check that the `.mat` file exists, load its `metadata` variable,
verify that metadata agrees with the expected `(chr, sex)` and corresponding
`.mat` filename, compute `file_md5` from the `.mat` file, and write a manifest
entry from the `.mat` metadata plus manifest-only `file` and `file_md5`.
Missing expected `.mat` files or bundled `reference_bim` files are errors before
`ld_manifest.json` is written. The finalizer must also build a `ReferencePanel`
from the bundled reference BIM files for the requested labels, save it with
`ReferencePanel.save_cache(...)` under `mat_root`, and record the cache filename
and MD5 in the MATLAB/Octave manifest.

MATLAB is required for production conversion because production artifacts must
be v7.3. Octave may write v5 sparse `.mat` files only when the caller
explicitly requests `format = "v5"`. MATLAB remains the normative runtime where
MATLAB and Octave differ.

The converter may assume `.npz` shards were produced and validated by
`statgen_build_ld.py`; it must validate metadata consistency before writing but
is not required to repeat expensive O(nnz) structure checks such as full
symmetry or diagonal scans. The converter copies bundled reference `.bim` files
named by each shard's `reference_bim` metadata from the `.npz` panel root into
the `.mat` output directory unchanged and carries `reference_bim` values
through to the MATLAB/Octave shard metadata.

### R LD preparation

R `.npz` support is an LD shard reader, not a general NumPy archive API. It
must support the documented LD shard payload dtypes (`float32`, `int32`,
`int64`, `uint8`, and metadata UTF-8 JSON bytes), and reject object arrays,
pickled payloads, unsupported dtypes, unsupported byte-order encodings,
malformed shapes, and archives with missing required members. The R `.npz`
reader must not require Python or external binaries. Signed 64-bit integer
payloads must preserve integer precision and must not be coerced through
double.

`statgen::prepare_ld_npz_for_r(npz_root, extract_npz = FALSE)` prepares a
Python `.npz` LD distribution for direct R loading. It reads the existing
`npz_root/ld_manifest.json`, validates that the manifest has
`runtime_format: "python_npz_csc32"`, builds an R-native `ReferencePanel` cache
from the bundled reference BIM files required by all shards described in that
manifest, saves that cache under `npz_root` as an RDS file, and updates the
existing Python manifest in place with `r_reference_cache` and
`r_reference_cache_md5`. It must not accept a shard filter, rewrite LD shard
files, copy LD shard files to another directory, create an R LD shard manifest,
or write a precomputed `ld_r2` payload.

When `extract_npz = TRUE`, `prepare_ld_npz_for_r` also creates or refreshes an
extracted cache directory next to each manifest-declared `.npz` shard. The
directory name is the shard filename plus `.d`, such as `ld_chr1.npz.d`, and it
contains the required `.npy` members plus a `.statgen-extracted-npz.json`
marker. The marker records `format: "statgen_extracted_npz"`, `source_file`,
and `source_file_md5`. R `load_ld` may use the extracted directory only when
the marker's `source_file_md5` matches the manifest `file_md5` for that shard
and all required `.npy` members are present. If the extracted directory is
absent, stale, or incomplete, `load_ld` must read the `.npz` shard normally,
using the scratch-directory behavior defined in `spec/performance-contract.md`.

`prepare_ld_npz_for_r` is idempotent. If `r_reference_cache` is already present
in the manifest, the function must validate that it is a plain relative
filename with no path separators, no `..` components, and no absolute path
prefix, then reuse that filename for the regenerated cache. If
`r_reference_cache` is absent, the function uses the default
`reference_cache.rds` filename. Each run rebuilds the R `ReferencePanel` cache
from bundled BIM files, overwrites the R cache file atomically, recomputes
`r_reference_cache_md5`, and updates the manifest in place. A stale or corrupt
existing R cache is not an error by itself; an invalid cache filename or an
unwritable destination is an error.

## In-memory objects

An `LDShard` holds one sparse signed-`r` matrix and its aligned `a1freq`. Each
`LDShard` is unambiguously identified by its `(chr, sex)` pair.

Runtime implementations may also construct and retain an internal sparse
`ld_r2` matrix at shard creation time. `ld_r2` is the element-wise square of
`ld_r` (`ld_r .^ 2` in MATLAB notation), not the matrix product
`ld_r * ld_r`. It is derived from the loaded distribution artifact, is not a
separate on-disk distribution field, and exists to make repeated
`multiply_r2` and pruning calls avoid rebuilding the same sparse squared LD
matrix.

- MATLAB/Octave and R `load_ld` support a
  `retain_ld_r = false` option for memory savings. Loading first constructs the
  raw signed `ld_r` matrix from the runtime distribution artifact, then
  constructs internal `ld_r2`, then drops `ld_r` when `retain_ld_r` is false.
  The option must default to retaining `ld_r`, preserving the ordinary LD
  object contract. When a caller disables raw-LD retention,
  `LDPanel.multiply_r2` and `fast_prune` must continue to work from `ld_r2`,
  but direct shard-level `ld_r` inspection is unavailable for that loaded
  object.

An `LDPanel` is an ordered collection of shard lists matching the paired
reference panel. For autosomal chromosome shards, each chromosome has exactly
one `LDShard` (`sex = null`). For chrX, the panel may hold any non-empty subset
of `"female"`, `"male"`, and `"combined"`. Loading always loads the full set
of shards selected from disk for each requested reference shard.

`LDPanel` carries a `default_chrX_sex` field (`"female"`, `"male"`, or
`"combined"`) that selects which chrX shard downstream operations use when they
need a single LD matrix for chrX. If chrX is present in the panel,
`default_chrX_sex` must name a shard actually present in that chrX shard group.
If chrX is absent, `default_chrX_sex` is retained as metadata and ignored by
chrX-unrelated operations. For default sex-specific panels built by
`statgen_build_ld.py`, the default is `"female"`.

`LDShard` fields:

- `chr`: contig/shard label (`1`-`22` or `X`);
- `sex`: `None`/`null` for autosomal/sex-agnostic shards; `"female"`,
  `"male"`, or `"combined"` for chrX shards;
- `num_snp`: number of SNPs in the shard;
- `ld_r`: language-native sparse signed-correlation matrix;
- `ld_r2`: optional internal language-native sparse matrix containing the
  element-wise square of `ld_r`;
- `a1freq`: vector of allele frequencies aligned to the paired
  `ReferenceShard`;
- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`.

Python stores `ld_r` as a SciPy CSC matrix from the `.npz` CSC components.
MATLAB/Octave stores `ld_r` as native sparse double loaded from `.mat`. R
stores `ld_r` in memory as `Matrix::dgCMatrix` constructed from the `.npz` CSC
components. Python retains `a1freq` as `float32` from the `.npz` payload;
MATLAB/Octave and R return `a1freq` as double.

## API

```text
load_ld(path, optional shards, optional default_chrX_sex) -> LDPanel
load_ld_reference(path, optional shards) -> ReferencePanel
validate_ld_distribution(path, optional check_payload_structure) -> report

LDPanel.num_snp -> int
LDPanel.reference -> ReferencePanel
LDPanel.a1freq(optional chrX_sex) -> num_snp float vector
LDPanel.default_chrX_sex -> "female" | "male" | "combined"
LDPanel.select_shards(shards) -> LDPanel
LDPanel.multiply_r2(M, optional chrX_sex) -> vector or matrix with same shape as M
fast_prune(logpvec, ld_panel, optional r2_threshold, optional chrX_sex) -> logpvec
```

`path` identifies a panel root directory containing `ld_manifest.json`.
`load_ld` must validate that every loaded LD shard is compatible with its
paired `ReferenceShard`: shard label, `num_snp`, and `reference_checksum` must
match. A mismatch is an error.

`load_ld` loads the manifest-declared reference cache, applies `shards`
subsetting to that reference and to the LD shard entries, and retains the loaded
reference as `LDPanel.reference`. Requesting a shard absent from the reference
cache or from the LD panel is an error. Missing required LD files or the
manifest-declared reference cache are errors. Bundled reference BIM files are
not read by `load_ld`.

`load_ld_reference` loads only the manifest-declared reference cache from an LD
distribution root, applies the same optional `shards` subsetting, and does not
load LD shard matrices. It is the convenience API for retrieving the reference
paired with an LD distribution when callers do not know the manifest-declared
cache filename.

`validate_ld_distribution` is an explicit distribution-QA path, not part of
default loading. It is self-contained: no external reference is required. It
validates only the calling runtime's LD distribution format: Python validates
`.npz` distributions with `runtime_format: "python_npz_csc32"`; MATLAB/Octave
validates `.mat` distributions with `runtime_format:
"matlab_mat_sparse_double"`; R validates `.npz` distributions with
`runtime_format: "python_npz_csc32"` plus the required
`r_reference_cache`/`r_reference_cache_md5` sidecar fields. In addition to
manifest file MD5 checksums and manifest/per-file metadata agreement, the
validator must confirm that the bundled reference representations agree with
each other and with all LD shard metadata: the manifest-declared reference
cache for the calling runtime must be present, every `reference_bim` file named
in the manifest must be present, and each LD shard's `reference_checksum` and
`num_snp` must match both the loaded reference cache and the corresponding
bundled `.bim` file. This double-check is intentional: the cache is the default
loading source, while the BIM files are portable build provenance. This
explicit QA path may parse bundled BIM files; default `load_ld` must not. When
`check_payload_structure` is true, it additionally performs expensive checks:
sparse index bounds, explicit diagonal, and symmetry validation. For
MATLAB/Octave `.mat`
distributions, v5 MAT-files are accepted for validation but must produce a
warning stating that they are fixture/local-test artifacts and not production
distribution artifacts because of MAT-file size limits.

LD-specific public cache APIs are not part of the LD contract. There is no
`save_ld_cache` or `load_ld_cache`; the runtime distribution artifacts are the
load-efficient representation.

Expected behavior:

- Panel accessors are read-only and concatenate shard data in reference panel
  order.
- `LDPanel.reference` is always populated after a successful `load_ld` call
  from the manifest-declared reference cache.
- Loading warns and proceeds when a shard metadata field
  `num_monomorphic_snps` is greater than zero. Runtime operations do not
  special-case those SNPs.
- `default_chrX_sex` defaults to `"female"` when omitted.
- `r2_threshold` defaults to `0.2` when omitted.
- `LDPanel.a1freq(optional chrX_sex)` uses `default_chrX_sex` for chrX by
  default. `chrX_sex` applies only to chrX; chr1-22 always use their
  sex-agnostic shard. If chrX is present, an override must name a chrX shard
  present in the panel. If chrX is absent, `chrX_sex` is ignored.
- `LDPanel.multiply_r2(M, optional chrX_sex)` computes
  `elementwise(LD_r .^ 2) * M` independently per reference shard, never
  materializes a dense genome-wide LD matrix, accepts both 1-D and 2-D inputs,
  and returns the same shape as `M`.
- Python should return `float32` for `float32` input and `float64` for
  `float64` input when multiplying by LD. MATLAB/Octave and R may return double
  because their native sparse LD matrices are double. Cross-runtime tests
  compare numerical values with dtype-aware tolerances, not dtype identity.
- `fast_prune` applies greedy significance-based pruning independently per
  reference shard using the LD selected by `chrX_sex`.
- `fast_prune` returns a floating vector with the same shape as `logpvec`.
  Python preserves `float32` and `float64` inputs and promotes non-floating
  inputs to `float64` so pruned values can be represented as `NaN`;
  MATLAB/Octave and R may return double.
- Shard-level matrix multiplication is an internal implementation detail, not
  part of the public API.

## r² Matrix Multiply

The core LD operation multiplies the per-shard element-wise squared LD matrix
by a genome-wide matrix or vector. This is the only numerical primitive exposed
by `LDPanel`; all downstream computations (LD scores, annotation LD weighting,
per-SNP variance) are expressed in terms of it.

```text
LDPanel.multiply_r2(M, optional chrX_sex) -> vector or matrix with same shape as M
```

`M` is a `num_snp x k` matrix or a `num_snp` vector aligned to the reference
panel row order. The result has the same shape as `M`.

`LDPanel.multiply_r2(M, optional chrX_sex)` splits `M` by reference shard,
computes `elementwise(LD_r .^ 2) * M_shard` for each shard, and concatenates
the results in reference panel order. This is not the matrix square
`LD_r * LD_r`. Omitting `chrX_sex` uses `default_chrX_sex` for the chrX shard.
Pass `"female"`, `"male"`, or `"combined"` to override in language-specific
syntax. The `chrX_sex` selector affects only chrX; autosomes are always
sex-agnostic. If chrX is absent from the panel, `chrX_sex` is ignored.

The operation must not materialize a dense genome-wide LD matrix. Callers
pre-exclude SNPs by setting the corresponding rows of `M` to zero.

## Pruning

```text
fast_prune(logpvec, ld_panel, optional r2_threshold, optional chrX_sex) -> logpvec
```

Greedy significance-based LD pruning. Input and output are genome-wide vectors
aligned to the reference panel.

Algorithm (applied independently per reference shard):

1. Stable-sort SNPs by `|logpvec|` descending; ties keep reference order. Skip
   `NaN` entries (treated as pre-excluded).
2. For each SNP in sorted order: if not already pruned, retain the current SNP
   and mark all other SNPs with `r² >= r2_threshold` as pruned (output set to
   `NaN`). The current SNP is never pruned by its own diagonal.
3. Retained positions keep their original `logpvec` values.

Omitting `chrX_sex` uses `default_chrX_sex` for the chrX shard; if chrX is
absent, the selector is ignored. Omitting `r2_threshold` uses `0.2`.
