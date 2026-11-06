# GenotypeShard and GenotypePanel

## Represents

Genotype inputs are PLINK 1 bfile-backed hardcall genotype sources. In memory,
`GenotypeShard` and `GenotypePanel` are reference-aligned access objects:
they hold sample metadata, source-row mapping, ploidy metadata, and `.bed`
paths, while hardcall genotype payloads remain on disk and are decoded on
demand.

Returned counts are oriented to reference `a1` under the allele contract in
[SPEC.md](SPEC.md).
Other genotype formats (PGEN, BGEN, VCF, dosage) are out of scope and should
be converted by project-specific preprocessing before entering `statgen`.

The names `GenotypeShard` and `GenotypePanel` are retained for consistency
with the rest of the object family, but the object is an accessor/metadata
handle rather than a resident genotype matrix.

## Disk representation

Canonical disk format is a PLINK bfile triplet identified by a
`bfile_prefix`:

```text
<bfile_prefix>.bed
<bfile_prefix>.bim
<bfile_prefix>.fam
```

The `.bim` file follows [reference.md](reference.md): six PLINK BIM columns
separated by tabs or by runs of ASCII whitespace, with no embedded whitespace in
fields.

A non-sharded single PLINK triplet combined across chromosomes is an
acceptable disk representation. For example, `bfile_prefix = "genotypes/all"`
identifies:

```text
genotypes/all.bed
genotypes/all.bim
genotypes/all.fam
```

A sharded panel uses `@` in `bfile_prefix` as the shard-label placeholder. For
example, `bfile_prefix = "genotypes/chr@"` identifies one triplet per loaded
reference shard:

```text
genotypes/chr1.bed
genotypes/chr1.bim
genotypes/chr1.fam
...
genotypes/chrX.bed
genotypes/chrX.bim
genotypes/chrX.fam
```

Each bfile may have an optional sidecar:

```text
<bfile_prefix>.ploidy
```

The `.ploidy` sidecar follows the `genomatch` BFILE `.ploidy` contract:
it is a tab-delimited file with no header, one row per BIM row in the same
order, with exactly two integer columns:

```text
male_ploidy  female_ploidy
```

Allowed ploidy values are `0`, `1`, and `2`. Missing `.ploidy` means autosomal
source rows are diploid in both sexes, equivalent to `(2, 2)`. Missing chrX
`.ploidy` defaults to `(1, 2)` and must warn, because this is an assumption;
users should provide `.ploidy` for PAR/non-PAR mixtures or other nonstandard
chrX encodings. For `@`-sharded bfiles, each resolved shard's `.ploidy` sidecar
is optional independently; it is normal for only the chrX shard to carry a
`.ploidy` file.

Converted dense genotype matrices are not a `statgen` storage format and are
not saved by genotype cache APIs.

## Loading and alignment

`load_genotype(bfile_prefix, reference, optional allow_indel_swaps=false)`
parses BIM, FAM, and optional PLOIDY
metadata, then projects the source variant rows into the supplied
`ReferencePanel` order. The result has the same shard structure and SNP axis as
`reference`.

Source-to-reference matching is defined in
[reference.md](reference.md#source-to-reference-matching). Source `chr` must
match the reference shard label exactly. Matching uses fixed-width allele
hashes and does not normalize chromosome labels or repair strand issues.
Source validation follows [contigs-and-shards.md](contigs-and-shards.md).

`load_genotype` is a reference-driven projection from a possibly larger source.
For non-sharded source input, the source BIM is split by `chr` into
per-reference-shard groups after contig-label validation. Validated source rows
whose `chr` label is outside the supplied reference shard set are ignored for
alignment, which allows loading a reference subset from a genome-wide source
file. Recognized non-supported contigs `Y` and `MT` are also ignored following
[contigs-and-shards.md](contigs-and-shards.md). Ambiguous or non-canonical
labels such as `chr1` remain errors. If a supplied reference shard has no
matching source BIM rows in the non-sharded source, the loader warns and keeps
that reference shard as entirely absent (`is_present == false`, ploidy `NaN`).

The aligned variant-presence vector is:

- `is_present`: logical vector, length equal to the paired reference shard
  or panel. `1` means a matching source BIM row is present after alignment;
  `0` means no source genotype row exists for that reference SNP.

The aligned `is_allele_swapped` logical vector has the same full reference
axis as `is_present`. It is true only for matched rows accepted by swap
fallback, and false for exact matches and absent SNPs. It records allele
orientation independently of whether any subject's hardcall is missing.

Within this genotype spec, `num_snp` in public accessor shapes means the
reference-axis SNP count of the paired `ReferenceShard` or `ReferencePanel`,
not the number of SNP rows in the source BIM. Source BIM row count is named
`source_num_snp`.

Ploidy vectors are aligned the same way:

- `ploidy_male`: numeric vector, length equal to `is_present`;
- `ploidy_female`: numeric vector, length equal to `is_present`.

For matched SNPs, ploidy values come from `.ploidy` when present. When
`.ploidy` is absent, autosomal rows use `(2, 2)` and chrX rows use `(1, 2)`.
For unmatched SNPs, `ploidy_male` and `ploidy_female` are `NaN` by definition.

Each `GenotypeShard` retains the paired `ReferenceShard` checksum. Genotype
caches save reference checksums and loaders trust those stored checksums in
the same way as sumstats and annotation caches.

## Sample metadata

FAM rows define the panel-level sample axis. Required FAM columns are exposed
as read-only panel-level accessors:

```text
fid
iid
father_id
mother_id
sex
is_male
is_female
```

The phenotype column from FAM is parsed for row-shape validation but is not a
public accessor in the genotype object contract. Sample order is the FAM file
order and defines the row order returned by genotype accessors.

FAM `sex` values are preserved as loaded and must be one of the PLINK values
`1`, `2`, or `0` for male, female, or unknown. Other values are errors.
`is_male` is `sex == 1`; `is_female` is `sex == 2`; unknown sex is represented
by `~is_male & ~is_female`.

Every FAM file must have unique `(fid, iid)` pairs. Duplicate pairs are errors.

## Sample axis reconciliation

For non-sharded genotype input, the single FAM file defines the panel-level
sample axis.

For sharded genotype input, loaded autosomal shards (`1`-`22`) define the
ordinary panel-level sample axis:

- all loaded autosomal shards must have identical FAM rows in the same order
  across `fid`, `iid`, `father_id`, `mother_id`, and `sex`;
- if one or more autosomal shards are loaded, their shared FAM row order
  defines `GenotypePanel.fid`, `iid`, `father_id`, `mother_id`, `sex`,
  `is_male`, and `is_female`.

ChrX is the only exception to identical-FAM sharded loading. The chrX FAM file
may be a subset of the panel-level sample axis:

- chrX subjects are matched to the panel sample axis by unique `(fid, iid)`;
- every chrX `(fid, iid)` pair must exist in the panel sample axis;
- for every matched chrX subject, `fid`, `iid`, `father_id`, `mother_id`, and
  `sex` must match the panel-level FAM values exactly;
- `GenotypePanel.is_subject_present("X")` returns a panel-axis logical vector
  that is true for subjects present in the chrX FAM and false otherwise.

If chrX is loaded without any autosomal shard, the chrX FAM file defines the
panel-level sample axis and `GenotypePanel.is_subject_present("X")` is all
true. This is the only case where chrX does not reconcile against an autosomal
sample axis.

For autosomal shards, `GenotypePanel.is_subject_present(shard)` always returns
an all-true vector aligned to the panel-level sample axis.

Genotype fetch accessors always return rows in the panel-level sample order.
For chrX SNPs, genotype rows read from the chrX `.bed` file are expanded onto
the panel-level sample axis by `(fid, iid)`; subjects absent from the chrX FAM
receive missing calls. This is the only sample-axis filling performed by
genotype fetch accessors.

## In-memory objects

A `GenotypeShard` is aligned to one `ReferenceShard` and contains:

- `chr`: paired reference shard label;
- `num_snp`: reference-axis SNP count for the paired reference shard, equal to
  `length(is_present)`;
- `bed_path`: source `.bed` path for that shard;
- `bed_file_size`: `.bed` file size in bytes recorded at load or cache-save
  time; at `GenotypeShard` construction (source load or cache load),
  implementations validate
  `3 + ceil(source_num_sample / 4) * source_num_snp == bed_file_size`
  as a metadata integrity check;
- `source_num_snp`: physical-source metadata: the total number of SNP rows in
  the `.bim`/`.bed` file identified by this shard's `bed_path`. This is the
  row-count bound for `source_row0` and the value used for BED file-size
  validation. It is not the number of source rows that matched this reference
  shard. For non-sharded bfile input, every `GenotypeShard` pointing to the
  shared `.bed` records the same full-file source SNP count;
- `source_num_sample`: number of subject rows in the source FAM shard;
- `source_row0`: integer vector, length equal to the paired reference shard,
  containing the zero-based source BIM row index for each matched SNP and `-1`
  for unmatched SNPs;
- `subject_present`: logical vector, length `num_sample`, indicating which
  panel-level subjects are present in the source FAM for this shard;
- `source_subject_row0`: integer vector, length `num_sample`, containing the
  zero-based source FAM row index for each panel-level subject present in this
  shard and `-1` for subjects absent from the shard.

The shard field is named `subject_present` (without `is_` prefix) because it
is an internal per-shard attribute, not a top-level panel accessor. The
panel-level method is `is_subject_present(shard)`, which takes a shard label
and returns the corresponding shard's vector. The asymmetry with the SNP axis
(`is_present` is a flat panel-wide vector) reflects the difference: SNP
presence is uniform in meaning across shards, while sample presence is
intrinsically per-shard because only chrX may carry a FAM subset.
- `is_present`, `is_allele_swapped`, `ploidy_male`, and `ploidy_female`
  aligned to the paired reference shard;
- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`. Genotype shards expose this field by name for reference
  compatibility checks.

`GenotypeShard` does not store FAM columns (`fid`, `iid`, `father_id`,
`mother_id`, `sex`). The source `.fam` is parsed temporarily during loading
to compute `subject_present` and `source_subject_row0`, then discarded.
Panel-level FAM accessors are held by `GenotypePanel` alone.

A `GenotypeShard` is the unit of BED access. It is paired 1:1 with one
`ReferenceShard` and owns the source BED mapping for that reference shard:
`bed_path`, `source_row0`, `source_num_snp`, `source_num_sample`,
`subject_present`, and `source_subject_row0`. For non-sharded bfile input,
multiple `GenotypeShard`s may point to the same source `bed_path`; their
`source_row0` vectors index different rows of the same source BED/BIM. For
`@`-sharded bfile input, each `GenotypeShard` normally points to the resolved
BED file for its paired reference shard. `GenotypePanel` operations split
panel-global SNP requests by `GenotypeShard`, call the corresponding
shard-level BED reader, and reassemble outputs in requested panel-global order.

A `GenotypePanel` is an ordered collection of `GenotypeShard` objects aligned
to a `ReferencePanel`. It retains paths and metadata only; dense genotype
matrix access is deferred to explicit accessor calls. For non-sharded bfile
input, `GenotypePanel` still produces one `GenotypeShard` per reference shard,
all pointing to the same source `.bed`; the in-memory and cache representation
is always shard-structured regardless of whether the source bfile was sharded.
`GenotypePanel.source_layout` records the physical source layout used when the
metadata was built: `"non_sharded"` for a single bfile source or `"sharded"` for
an `@` source. Shard subsetting and cache subsetting preserve this value.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: each shard's vector length and row order correspond to the
paired `ReferenceShard`, and any cache is valid only for that reference.

Source BIM and FAM parsing is performance-sensitive tabular input:

- Python implementations use `pandas` readers with explicit column names and
  dtypes where needed.
- MATLAB/Octave implementations use `textscan` with explicit schemas, matching
  [matlab.md](matlab.md).
- R implementations use language-native tabular readers with explicit column
  schemas.

Genotype hardcall decoding must be vectorized across subjects for each SNP.

## Cache layout

Genotype caches save metadata and aligned vectors only. They do not save
genotype hardcalls. Python genotype caches use NumPy `.npz` format.

MATLAB/Octave genotype caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
is_present
is_allele_swapped
ploidy_male
ploidy_female
source_row0
subject_present
source_subject_row0
fid
iid
father_id
mother_id
sex
is_male
is_female
```

`metadata` is a struct with:

```text
schema = "genotype_cache/0.2"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
source_layout
bed_paths
bed_file_sizes
source_num_snp
source_num_sample
num_sample
```

`is_present` and `is_allele_swapped` are panel-wide logical vectors of equal
length. A true swap flag requires `is_present == true`; cache loaders validate
this invariant. Writers emit `/0.2` and persist swap state. Loaders
also accept `genotype_cache/0.1`, initializing all swap flags to false without
reparsing sources or recomputing matches. Shard selection preserves the stored
mapping and swap flags on the selected reference axis. Caches store the
resulting orientation, not the `allow_indel_swaps` loading option.

`ploidy_male` and `ploidy_female` are panel-wide numeric vectors with `NaN`
for `is_present == false`.
`source_row0` is a panel-wide integer vector with `-1` for
`is_present == false`. FAM columns and sex masks are panel-level vectors with
length `num_sample`. `subject_present` and `source_subject_row0` are stored per
shard as `num_sample × n_shards` panel-axis matrices, or an equivalent
language-native representation that preserves the same values. `bed_paths`,
`bed_file_sizes`, `source_num_snp`, and `source_num_sample` in `metadata` are
vectors of length `n_shards` because they are stored on reconstructed
`GenotypeShard`s, but their meaning is physical-source metadata for the
corresponding `bed_path`. For non-sharded source bfiles, the repeated entries
record the shared `.bed` path and the same full-file source SNP count for each
shard, not shard-local match counts. `source_layout` is panel-level metadata and
must be preserved on cache load, including optional shard subsetting. Shard
offsets are zero-based half-open intervals into the panel-wide SNP-axis vectors
and are sufficient to reconstruct `GenotypeShard` objects. Cache metadata
validation should be cheap, depending on shard count, vector dimensions, sample
count, stored physical-source `source_num_snp`, `source_num_sample`, recorded
`.bed` file sizes, and `source_layout` rather than reparsing BIM/FAM or PLOIDY
source files.

R genotype caches are RDS files containing one named list with the same logical
top-level fields and metadata fields as the MATLAB/Octave cache layout above.
R `metadata` is a named list. `is_present`, `is_allele_swapped`,
`subject_present`, `is_male`, and `is_female` are logical vectors or matrices.
`ploidy_male` and `ploidy_female` are numeric vectors. `source_row0` and `source_subject_row0` are integer vectors
or matrices with `-1` for absent SNPs or subjects. FAM columns are character
vectors; `sex` is an integer or numeric vector preserving FAM sex codes.

`load_genotype_cache(path, optional shards)` restores the metadata/accessor
object from cache. It must not parse source `.bim`, `.fam`, or `.ploidy`
files, and it should not require the original `.bed` files to be present at
cache-load time because callers may later supply an alternate `bed_path` to
`fetch_genotypes_int8` or `fetch_genotypes`. BED existence and file-size checks
are performed lazily by genotype fetch accessors for the addressed shards. The
cache does not verify BED payload checksums because PLINK `.bed` has no
manifest and full checksum verification is not part of the genotype cache
contract.

## Genotype access

```text
GenotypePanel.fetch_genotypes_int8(snp_indices, optional bed_path) -> matrix
GenotypePanel.fetch_genotypes(snp_indices, optional bed_path, optional haploid_mode) -> matrix
```

`snp_indices` are panel-global SNP indices in the host language's ordinary
index base. The returned matrix has shape:

```text
num_sample × length(snp_indices)
```

Columns are returned in exactly the requested `snp_indices` order.
`fetch_genotypes_int8` is the core hardcall accessor. It returns an `int8`
matrix with reference-`a1` counts on the raw diploid scale, encoded as `0`,
`1`, and `2`, and missing hardcalls encoded as `-1`. After decoding source-`a1`
counts, replace each nonmissing call `x` with `2 - x` where
`is_allele_swapped` is true. Missing calls and absent subjects remain missing.
This orientation step occurs before any ploidy scaling on each fetch; it does
not mutate BED data or panel metadata.

`fetch_genotypes` wraps `fetch_genotypes_int8` and returns a double-precision
matrix. Its `haploid_mode` argument defaults to `"raw"`. In `"raw"` mode,
`0`, `1`, and `2` are converted to double and missing hardcalls are represented
as `NaN`. In `"ploidy_scaled"` mode, non-missing calls are transformed as
`raw * ploidy / 2`, using `ploidy_male` for FAM sex `1` and `ploidy_female`
for FAM sex `2`. It maps PLINK diploid-style hardcall encodings onto the
declared biological ploidy. So male chrX with ploidy `1` turns `0/2` into
`0/1`, while diploid calls stay `0/1/2`.

`haploid_mode` is a per-call output option, not object state. In
`"ploidy_scaled"` mode, the call fails if any addressed subject has unknown sex
and a non-missing addressed call for a SNP where `ploidy_male` differs from
`ploidy_female`. Unknown-sex subjects absent from a shard do not trigger this
guard because their fetched calls are missing before and after scaling. When a
`bed_path` override is supplied, scaling still uses the loaded or cached ploidy
metadata, because PLINK `.bed` files do not carry ploidy metadata.

Every requested SNP must have `is_present == true`; otherwise the call fails
clearly before reading genotype payloads. Repeated SNP indices are allowed and
return repeated genotype columns.

`fetch_genotypes_int8` reads the underlying `.bed` file on each call. It
validates the PLINK BED magic bytes and requires SNP-major mode. For each SNP,
it reads only the corresponding packed row from disk and decodes it vectorized
across all subjects. Under the `statgen` allele contract,
the BED decoder first converts PLINK two-bit values to source-`a1` counts as:

```text
00 -> 2
01 -> -1
10 -> 1
11 -> 0
```

The optional `bed_path` argument lets callers read genotype payloads from a
different location while reusing the cached/aligned metadata. It may identify a
single `.bed` file or a sharded path containing `@`, resolved by the same
shard-label substitution rules as source loaders. A single `.bed` override is
used for all addressed shards only when those addressed shards share the same
recorded physical source layout, defined as equal `source_num_snp` and equal
`source_num_sample` across all addressed shards. An `@` override resolves one
`.bed` path per addressed shard and preserves sharded physical layout.

Overrides must not convert between source layouts. Enforcement uses
`GenotypePanel.source_layout`, not equality of stored `bed_path` values, because
single-shard subsets of sharded metadata still have sharded physical layout:

- `source_layout == "non_sharded"`: a flat override is accepted; an
  `@`-containing override is an error.
- `source_layout == "sharded"`: an `@`-sharded override is accepted; a flat
  override is accepted only when all addressed shards pass the layout equality
  check (`source_num_snp` and `source_num_sample` equal across all addressed
  shards), otherwise it is an error.

The alternate `.bed`
file(s) must have the same FAM/BIM structure and row order as the metadata
object, including source allele orientation and any chrX sample subset.
Stored swap flags still apply to a BED override. Since `.bed` files contain no
manifest, the accessor can only validate basic structure: file existence, PLINK
magic bytes, SNP-major mode, and file size for each addressed shard. The
accessor validates that the actual disk size of each addressed `.bed` file
equals the `bed_file_size` stored in the shard metadata; size mismatch is an
error.

## API

```text
load_genotype(bfile_prefix, reference, optional allow_indel_swaps=false)
    -> GenotypePanel
save_genotype_cache(panel, path) -> void
load_genotype_cache(path, optional shards) -> GenotypePanel

GenotypePanel.num_snp -> int
GenotypePanel.is_present -> num_snp logical vector
GenotypePanel.is_allele_swapped -> num_snp logical vector
GenotypePanel.ploidy_male -> num_snp float vector
GenotypePanel.ploidy_female -> num_snp float vector
GenotypePanel.num_sample -> int
GenotypePanel.fid -> num_sample string vector
GenotypePanel.iid -> num_sample string vector
GenotypePanel.father_id -> num_sample string vector
GenotypePanel.mother_id -> num_sample string vector
GenotypePanel.sex -> num_sample vector preserving FAM values 1, 2, or 0
GenotypePanel.is_male -> num_sample logical vector
GenotypePanel.is_female -> num_sample logical vector
GenotypePanel.is_subject_present(shard) -> num_sample logical vector
GenotypePanel.fetch_genotypes_int8(snp_indices, optional bed_path)
    -> num_sample × len(snp_indices) int8 matrix
GenotypePanel.fetch_genotypes(snp_indices, optional bed_path, optional haploid_mode)
    -> num_sample × len(snp_indices) double matrix
GenotypePanel.source_layout -> "non_sharded" | "sharded"
GenotypePanel.select_shards(shards) -> GenotypePanel
GenotypePanel.save_cache(path) -> void
```

Expected behavior:

- `reference` is required; rows are projected into reference order on load.
- `bfile_prefix` is a PLINK bfile prefix without the `.bed`, `.bim`, or `.fam`
  suffix. Without `@`, it identifies one non-sharded PLINK bfile triplet
  spanning the requested reference shards. With `@`, it is a sharded bfile
  prefix template; `@` is replaced by each reference shard label. The resolved
  prefixes identify `.bed`, `.bim`, `.fam`, and optional `.ploidy` sidecars.
  Every requested reference shard must have a corresponding resolved source
  bfile shard. Missing `@` source shards are errors; users should call
  `ReferencePanel.select_shards` before `load_genotype` when loading a shard
  subset.
- Shard discovery, contig validation, and row order follow
  [contigs-and-shards.md](contigs-and-shards.md). Source loading does not
  accept a separate `shards` selector; the supplied `reference` defines the
  requested shard set.
- Source loaders validate `.bed`, `.bim`, and `.fam` presence for each source
  shard present on disk; optional `.ploidy` presence; BED magic bytes;
  SNP-major mode; exact expected BED file size; BIM/FAM/PLOIDY row counts; and
  FAM sample-axis reconciliation rules.
- Matching follows the shared source-to-reference matching contract in
  [reference.md](reference.md#source-to-reference-matching). The loader does
  not normalize chromosome labels or perform strand handling.
- `GenotypePanel` does not expose `bim` as a public accessor. BIM is parsed
  only to align source rows to the reference and to construct internal
  `source_row0`.
- `GenotypePanel.num_snp` returns the total reference-axis SNP count across all
  loaded shards, equal to the length of `is_present`.
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native vectors or matrices.
- `GenotypePanel.is_subject_present(shard)` accepts one loaded shard label and
  returns a panel-axis logical vector. It fails on unknown or unloaded shard
  labels.
- `is_present`, `ploidy_male`, and `ploidy_female` are full reference-axis
  vectors; missing source variants are represented by `is_present == false` and
  ploidy `NaN`, not by row drops.
- Cache save/load preserves only metadata needed to avoid reparsing BIM, FAM,
  and PLOIDY files and to fetch hardcalls from BED on demand.
- `GenotypePanel.save_cache(...)` is a thin convenience method equivalent to
  `save_genotype_cache(panel, ...)`.
- `load_genotype_cache` performs cache-internal validation only and supports
  optional `shards` subsetting; per-shard reference checksums and BED file-size
  metadata are trusted from cache metadata until a genotype fetch addresses the
  corresponding shard.
- `GenotypePanel.fetch_genotypes_int8` reads hardcalls on demand, returns raw
  diploid-scale `int8` reference-`a1` counts with `-1` for missing calls, and
  fails on any requested SNP where `is_present == false`.
- `GenotypePanel.fetch_genotypes` is a double-precision wrapper over
  `fetch_genotypes_int8`, converting `-1` to `NaN` by default. Its
  `haploid_mode` argument accepts `"raw"` and `"ploidy_scaled"`; the latter
  applies `raw * ploidy / 2` using per-SNP ploidy metadata and FAM sex. Unknown
  sex is an error only for non-missing fetched calls where male and female
  ploidy differ.
- Both genotype fetch accessors always return rows in panel-level sample order.
  For chrX shards with a FAM subset, subjects absent from the chrX source FAM
  receive missing calls (`-1` for `fetch_genotypes_int8`, `NaN` for
  `fetch_genotypes`).
- `GenotypePanel.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md) and preserves the sample axis
  unchanged.
- Cache payloads store per-shard reference checksums so compatibility with a
  `ReferencePanel` can be checked after load via
  `ReferencePanel.is_object_compatible`.
