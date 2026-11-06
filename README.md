# statgen

`statgen` provides Python, MATLAB/Octave, and R tools for working with
reference-aligned statistical genetics data: reference variants, genotypes, LD,
annotations, and GWAS summary statistics.

All panel objects use reference coordinates. Genome build, contig naming, and
strand should be resolved upstream by
[genomatch](https://github.com/precimed/genomatch) or an equivalent pipeline.

## Where to Start

- Preparing inputs and caches from PLINK, summary statistics, and BED files:
  [Tutorial 1](docs/TUTORIAL_1_PREPARE_DATA.md) builds an LD distribution and
  saves Python and MATLAB caches from source files.
- Running Python analyses:
  [Tutorial 2](docs/TUTORIAL_2_PYTHON.md) runs LD operations, prunes variants,
  and fetches genotypes.
- Running MATLAB analyses:
  [Tutorial 3](docs/TUTORIAL_3_MATLAB.md) runs the same analysis workflow in
  MATLAB.
- Running R analyses:
  [Tutorial 4](docs/TUTORIAL_4_R.md) runs the same workflow with the R package.
- Looking up MATLAB API help interactively:
  `help statgen`.
- Finding installation commands, core concepts, or the API map:
  [Input Layouts](#input-layouts), [Concepts That Matter](#concepts-that-matter),
  [Setup](#setup), and [Key APIs](#key-apis).
- Following changes or getting help:
  [Changelog](CHANGELOG.md) and [issue tracker](https://github.com/precimed/statgen/issues).
- Developer docs for implementing or changing `statgen` behavior:
  [spec/SPEC.md](spec/SPEC.md).

## Object Overview

| Object | Meaning | Typical input |
| --- | --- | --- |
| `ReferencePanel` | Ordered reference variant table | External PLINK `.bim`, either one file or an `@`-sharded path |
| `LDPanel` | Sparse LD distribution tied to a matching reference | Built LD distribution directory produced by `statgen` from PLINK bfiles |
| `AnnotationPanel` | Numeric annotation matrix painted onto the reference | External BED or BED-like annotation files |
| `GenotypePanel` | Reference-aligned PLINK genotype metadata with on-demand hardcall access | External PLINK 1 bfile prefix, either one bfile or an `@`-sharded prefix |
| `Sumstats` | One aligned summary-statistics trait or source | External `.tsv.gz` with `chr`, `bp`, `a1`, `a2`, `p`, and optional columns such as `beta`, `se`, `z`, and `n` |

## Input Layouts

- Reference input may be one `.bim` file or an `@`-sharded `.bim` template.
- Genotype input may be one PLINK bfile prefix or an `@`-sharded bfile prefix.
- LD input is a distribution directory with `ld_manifest.json`, not a single
  shard file or an `@` template.
- Sumstats input is one TSV/TSV.GZ file per trait or source, not a sharded
  input.
- Annotation input is one BED file per annotation column; multiple BED files
  are annotations, not chromosome shards.
- Cache loaders read one cache file and may optionally return a logical shard
  subset.

## Concepts That Matter

- All analysis objects are aligned to a `ReferencePanel`.
- Supported contig labels are `1`-`22` and `X`; labels are not normalized, so
  `chr1`/`chrX` inputs should be fixed upstream.
- LD is a directory distribution with an `ld_manifest.json` file, not a normal
  object cache.
- LD should be built from an unrelated, single-ancestry genotype reference panel
  appropriate for the downstream analyses.
- Python builds the canonical `ld_npz/` distribution. R prepares that same
  `.npz` distribution in place with an R reference-cache sidecar and optional
  extracted shard caches. MATLAB/Octave uses a separate `ld_mat/` distribution
  converted from the Python-built LD shards because it needs native sparse
  `.mat` files.
- LD panels preserve the full reference SNP axis. If a forced LD build contains
  monomorphic SNPs, LD involving those SNPs is undefined and represented by
  omitted sparse entries.
- Genotype caches store metadata and source BED paths, not dense genotypes.
- chrX LD sex-specific shards use FAM sex labels.
- Genotype `.ploidy` sidecars belong to genotype access, not LD reference
  inputs.
- Missing genotype `.ploidy` sidecars imply `(2, 2)` ploidy for autosomes and
  `(1, 2)` for chrX, with a warning for chrX. Provide `.ploidy` for PAR/non-PAR
  mixtures or nonstandard chrX encodings.
- Shard objects are internal; users work with panel objects and `select_shards`.

## Allele Matching

Results describe the reference panel's `a1` (ALT) allele. Loading automatically
adjusts reversed SNP alleles. Swapping insertion/deletion alleles requires the
`allow_indel_swaps` option.

See the [allele matching guide](docs/ALLELE_MATCHING.md) for examples, allele
roles, and genotype cache migration.

## Setup

### Python

Install from a cloned repository:

```sh
pip install -e python/
```

Then import APIs from the relevant submodule, for example:

```python
from statgen.reference import load_reference
from statgen.sumstats import load_sumstats
```

### MATLAB

From a cloned repository, add the MATLAB package folder to the path:

```matlab
addpath('/path/to/statgen/matlab')
```

Alternatively, download the MATLAB bundle from a GitHub release and add that
bundle's `matlab` folder to the path.

MATLAB help text is available for public functions, public panel classes, and
their methods:

```matlab
help statgen                           % package overview
help statgen.load_reference            % function help
help statgen.ReferencePanel            % class help
help statgen.ReferencePanel.select_shards
help statgen.LDPanel.multiply_r2       % method help
```

Octave is supported for small-scale testing and CI compatibility checks. Use
MATLAB for production-scale MATLAB workflows, especially large LD `.mat`
distributions. In Octave, the package overview is available with
`help statgen.Contents`.

## Key APIs

### ReferencePanel

- `load_reference(path)` loads a `.bim` file or an `@`-sharded `.bim` pattern.
- `panel.chr`, `panel.snp`, `panel.bp`, `panel.a1`, and `panel.a2` expose
  reference-coordinate variant fields.
- `panel.select_shards(shards)` subsets by shard label.
- `panel.is_object_compatible(other)` checks reference compatibility and returns
  `false` on mismatch rather than raising.
- `panel.save_cache(path)` and `load_reference_cache(...)` save and reload
  reference caches.

### LDPanel

- `load_ld(path, shards=None, default_chrX_sex=None)` loads a sparse LD
  distribution and its manifest-declared reference cache.
- `load_ld_reference(path, shards=None)` loads only the manifest-declared
  reference cache from an LD distribution root.
- `validate_ld_distribution(path)` checks manifest, checksums, bundled
  references, and optionally sparse payload structure.
- `panel.a1freq(chrX_sex=None)` returns reference-aligned allele frequencies.
- `panel.multiply_r2(M, chrX_sex=None)` multiplies by LD `r²`.
- `panel.select_shards(shards)` subsets by shard label.
- `fast_prune(logpvec, ld_panel, r2_threshold=0.2, chrX_sex=None)` performs LD
  pruning using aligned scores.

### Building LD Distributions

The LD workflow intentionally differs slightly by runtime:

```text
Python builds ld_npz/
R prepares ld_npz/ in place
MATLAB/Octave converts ld_npz/ to ld_mat/
```

`ld_npz/` is the canonical build and interchange format. R reads its `.npz`
shards directly after adding `reference_cache.rds`; with `extract_npz = TRUE`,
R can also create sibling `*.npz.d` extracted caches for faster repeated loads.
MATLAB/Octave uses a separate `ld_mat/` directory because its runtime sparse
matrix representation is stored in `.mat` files.

- `python script/statgen_build_ld.py ... --shard SHARD` builds one LD shard.
- `python script/statgen_create_ld_manifest.py --ld PATH` finalizes a Python LD
  distribution.
- R loads the Python `.npz` LD distribution directly after one R-specific
  preparation step: `statgen::prepare_ld_npz_for_r(npz_root)`. Use
  `extract_npz = TRUE` to create sibling `*.npz.d` extracted caches for faster
  repeated R loads.
- `statgen.convert_ld_npz_to_mat(input_root, output_root, shard)` converts a
  Python LD shard for MATLAB.
- `statgen.create_ld_mat_manifest(input_root, output_root, shards)` finalizes a
  MATLAB LD distribution after conversion.

### AnnotationPanel

- `load_annotations(bed_paths, reference)` paints BED intervals onto a
  `ReferencePanel`; BED basenames become unique annotation names and columns are
  marked binary.
- `load_annotation(path, reference, ...)` loads one BED-like annotation source;
  3-column inputs produce one binary column and selected numeric value columns
  produce continuous annotations.
- `create_annotations(reference, annomat, annonames, is_binary=None,
  annotation_metadata=None)` creates an annotation panel from an already aligned
  matrix, not raw BED input.
- `create_annotation(reference, annovec, annotation_name, is_binary=None,
  annotation_metadata=None)` creates a one-column annotation panel from an
  already aligned vector.
- `panel.annomat` exposes the reference-aligned sparse numeric annotation
  matrix.
- `panel.is_binary` marks binary columns, and `panel.annotation_metadata`
  exposes per-column metadata strings stored with the annotations.
- `panel.select_shards(shards)` and `panel.select_annotations(names)` subset an
  annotation panel.
- `panel.union_annotations(other, mode="by_name")` combines annotation panels.
- `panel.save_cache(path)` and `load_annotations_cache(...)` save and reload
  painted annotation caches.

### GenotypePanel

- `load_genotype(bfile_prefix, reference, allow_indel_swaps=False)` loads PLINK 1
  genotype metadata and aligns it to a `ReferencePanel`. `bfile_prefix` may be a single prefix or an
  `@`-sharded prefix.
- `panel.fetch_genotypes(snp_indices)` returns selected genotype calls without
  loading the full BED file. Counts are oriented to the reference panel's `a1`
  (ALT) allele on the raw diploid scale by default; use `ploidy_scaled` to map
  male chrX ploidy `1` calls from `0/2` to `0/1`. In MATLAB, `snp_indices` are one-based panel SNP indices; requested
  SNPs must be present in the genotype source.
- `panel.fetch_genotypes_int8(snp_indices)` returns selected hardcalls as compact
  integer calls.
- `panel.is_present` is a boolean mask in reference coordinates for variants
  present after alignment; `panel.is_allele_swapped` marks accepted swaps.
- `panel.fid`, `panel.iid`, `panel.sex`, `panel.is_male`, and
  `panel.is_female` expose sample metadata.
- `panel.save_cache(path)` and `load_genotype_cache(...)` save and reload
  genotype metadata caches for faster repeated loading.

### Sumstats

- `load_sumstats(path, reference, allow_indel_swaps=False)` loads a `.tsv.gz`
  summary-statistics file using the [allele matching policy](#allele-matching).
- Input files require `p`; `z` and `n` are optional and warn when absent or
  incomplete. Common aliases such as `POS`, `EffectAllele`, and `OtherAllele`
  are accepted.
- `create_sumstats(reference, pvec, zvec=None, nvec=None, ...)` creates a
  `Sumstats` object from already aligned vectors, not raw TSV input.
- `sumstats.logpvec`, `sumstats.zvec`, and `sumstats.nvec` expose aligned
  vectors. Optional vectors are `None` in Python or empty in MATLAB when
  absent.
- `sumstats.beta_vec`, `sumstats.zvec`, and `sumstats.eaf_vec` are oriented
  to the reference panel's `a1` (ALT) allele; EAF is its frequency.
- `sumstats.is_present` is a boolean mask in reference coordinates for variants
  present after alignment.
- `sumstats.select_shards(shards)` subsets by shard label.
- `sumstats.save_cache(path)` and `load_sumstats_cache(...)` save and reload
  aligned summary-statistics caches.

## License

`statgen` is licensed under the MIT License. Contributions are accepted under
the same license. See [LICENSE](LICENSE) and [CONTRIBUTING.md](CONTRIBUTING.md).
