# Tutorial 1: Prepare LD and Object Caches

This tutorial prepares the files used by the Python, MATLAB/Octave, and R
analysis tutorials:

- build a sparse LD distribution from PLINK bfiles;
- convert the LD distribution for MATLAB/Octave;
- load reference, summary statistics, annotations, and genotype metadata from
  disk representations;
- save Python and MATLAB/Octave caches for repeated analyses.

The repository includes small synthetic source fixtures under
`docs/tutorial_1_fixtures/source/`. Run the tutorial commands from
`docs/tutorial_1_fixtures/` to use those fixtures unchanged. Tutorial 1 reads
from `source/`, writes reusable artifacts to `derived/`, and uses `scratch/`
for disposable LD build files.

GitHub releases also attach a `statgen-tutorial-1-fixtures-<version>.zip`
archive with the same `source/` directory. Extract that archive into an empty
working directory and run the tutorial commands from that directory.

For example, without cloning the repository, set `VERSION` to the release
version you want to use:

```bash
VERSION="<release-version>"
mkdir statgen-tutorial-1-fixtures
cd statgen-tutorial-1-fixtures
curl -LO "https://github.com/precimed/statgen/releases/download/v${VERSION}/statgen-tutorial-1-fixtures-${VERSION}.zip"
unzip "statgen-tutorial-1-fixtures-${VERSION}.zip"
```

## Contents

- [Prerequisites](#prerequisites)
- [Inputs and Paths](#inputs-and-paths)
- [Build the LD Distribution](#build-the-ld-distribution)
- [Convert LD for MATLAB/Octave](#convert-ld-for-matlaboctave)
- [Build Python Caches](#build-python-caches)
- [Build MATLAB/Octave Caches](#build-matlaboctave-caches)
- [Outputs](#outputs)

## Prerequisites

Start from the repository root and move into the checked-in fixture directory:

```bash
cd docs/tutorial_1_fixtures
REPO_ROOT="../.."
```

The LD build step assumes `plink2` is installed and available on `PATH`:

```bash
plink2 --version
```

Use a Python environment with `statgen` installed from the repository:

```bash
conda create -n statgen python=3.11 numpy scipy pandas pytest -y
conda activate statgen
pip install -e "$REPO_ROOT/python/"
```

For MATLAB/Octave cache preparation, use MATLAB or Octave with the repository's
MATLAB package folder on the path:

```matlab
addpath('../../matlab')
```

## Inputs and Paths

```text
PLINK LD reference bfiles:
source/ld_reference/chr21.{bed,bim,fam}
source/ld_reference/chr22.{bed,bim,fam}
source/ld_reference/chrX.{bed,bim,fam}

PLINK genotype bfiles for on-demand genotype access:
source/genotypes/chr21.{bed,bim,fam}
source/genotypes/chr22.{bed,bim,fam}
source/genotypes/chrX.{bed,bim,fam}

Optional genotype chrX ploidy sidecar:
source/genotypes/chrX.ploidy

Summary statistics:
source/sumstats/trait_a.tsv.gz
source/sumstats/trait_b.tsv.gz

Annotations:
source/annotations/coding_exon.bed
source/annotations/conservation.annot
source/annotations/conservation.meta
source/annotations/exon.bed
source/annotations/functional_groups.annot
source/annotations/intron.bed
source/annotations/utr3.bed
source/annotations/utr5.bed
source/annotations/whole_gene.bed
```

The optional chrX `.ploidy` sidecar belongs to genotype bfiles, not to LD
reference bfiles. It supplies male/female ploidy metadata for chrX genotype
access.

The bfile prefixes are sharded by chromosome, so commands can use an `@`
placeholder:

```bash
BFILE="source/ld_reference/chr@"
GENOTYPE_BFILE="source/genotypes/chr@"
LD_NPZ_DIR="derived/ld_npz"
LD_MAT_DIR="derived/ld_mat"
PY_CACHE="derived/python_cache"
MAT_CACHE="derived/matlab_cache"
SCRATCH="scratch/ld_build"
```

Prepare inputs so genome build, chromosome names, and DNA strand agree across
files before using `statgen`.

Summary-statistics files must contain `chr`, `bp`, `a1`, `a2`, and `p` columns.
The `p` column is required, but empty fields, `NA`, and `NaN` are accepted as
missing values. Rows with missing `p` are validated, then dropped before
matching; their reference slots remain missing unless another valid row matches.
At `info` verbosity the loader reports the number dropped. Non-missing `p`
values must be finite and in `[0, 1]`; `p == 0` becomes `Inf` in `logpvec`.
`z` and `n` are optional; when present, missing values among variants with valid
`p` values will produce warnings.

Annotation examples include binary 3-column BED files, a grouped BED-like file,
and a headered BED-like continuous annotation file. `functional_groups.annot`
uses one group column to produce multiple binary annotation columns.
`conservation.meta` is a column-aligned metadata sidecar for
`conservation.annot`; headerless continuous files are also supported by
`load_annotation`.

## Build the LD Distribution

Build one chromosome per job. The `--shard` argument makes the same command
suitable for scheduler array jobs.

```bash
mkdir -p "$LD_NPZ_DIR" "$SCRATCH"

for shard in 21 22 X; do
  python "$REPO_ROOT/script/statgen_build_ld.py" \
    --bfile "$BFILE" \
    --out "$LD_NPZ_DIR" \
    --shard "$shard" \
    --scratch "$SCRATCH/shard_${shard}"
done
```

Autosomes write files such as `ld_chr21.npz` and `ld_chr22.npz`. chrX writes
sex-specific files by default, such as `ld_chrX_female.npz` and
`ld_chrX_male.npz`. The LD build uses FAM sex labels for chrX sex-specific
shards; a `.ploidy` sidecar is not part of LD reference preparation.

After all shard jobs finish, create the Python LD manifest:

```bash
python "$REPO_ROOT/script/statgen_create_ld_manifest.py" --ld "$LD_NPZ_DIR"
```

The manifest step writes `ld_manifest.json` and validates the completed LD
distribution.

## Convert LD for MATLAB/Octave

MATLAB/Octave loads native sparse `.mat` LD distributions. Convert the Python
`.npz` handoff files first.

```matlab
addpath('../../matlab')

ld_npz_dir = 'derived/ld_npz';
ld_mat_dir = 'derived/ld_mat';

shards = {'21', '22', 'X'};
for i = 1:numel(shards)
    statgen.convert_ld_npz_to_mat(ld_npz_dir, ld_mat_dir, shards{i});
end

statgen.create_ld_mat_manifest(ld_npz_dir, ld_mat_dir, shards);
statgen.validate_ld_distribution(ld_mat_dir, false);
```

The default conversion writes v7.3 MAT-files under MATLAB. Octave users may use
`'format', 'v5'` for fixture-scale validation.

## Build Python Caches

This pass loads source BIM, TSV, BED/BED-like annotation, and genotype metadata
files and writes reference-aligned caches.

```python
from pathlib import Path

from statgen.annotations import load_annotation, load_annotations
from statgen.genotype import load_genotype
from statgen.reference import load_reference
from statgen.sumstats import load_sumstats

cache_root = Path("derived/python_cache")
cache_root.mkdir(parents=True, exist_ok=True)

reference = load_reference(
    "source/ld_reference/chr@.bim",
    shards=["21", "22", "X"],
)

trait_a = load_sumstats("source/sumstats/trait_a.tsv.gz", reference)
trait_b = load_sumstats("source/sumstats/trait_b.tsv.gz", reference)

annot_root = Path("source/annotations")
binary_annotations = load_annotations(
    [
        annot_root / "coding_exon.bed",
        annot_root / "exon.bed",
        annot_root / "intron.bed",
        annot_root / "utr3.bed",
        annot_root / "utr5.bed",
        annot_root / "whole_gene.bed",
    ],
    reference,
)
grouped_annotations = load_annotation(
    annot_root / "functional_groups.annot",
    reference,
    has_header=True,
    group_column="group",
)
continuous_annotations = load_annotation(
    annot_root / "conservation.annot",
    reference,
    has_header=True,
    value_columns=["conservation", "promoter_activity"],
    annotation_metadata_path=annot_root / "conservation.meta",
)
annotations = binary_annotations.union_annotations(grouped_annotations).union_annotations(continuous_annotations)

genotype = load_genotype("source/genotypes/chr@", reference)

assert reference.is_object_compatible(trait_a)
assert reference.is_object_compatible(trait_b)
assert reference.is_object_compatible(annotations)
assert reference.is_object_compatible(genotype)

reference.save_cache(cache_root / "reference.npz")
trait_a.save_cache(cache_root / "trait_a.sumstats.npz")
trait_b.save_cache(cache_root / "trait_b.sumstats.npz")
annotations.save_cache(cache_root / "annotations.npz")
genotype.save_cache(cache_root / "genotype.npz")
```

The loaders above use the default `allow_indel_swaps=False`. If indel swaps
are appropriate for your inputs, set the option on each source load before
saving caches, for example:

```python
trait_a = load_sumstats("source/sumstats/trait_a.tsv.gz", reference, allow_indel_swaps=True)
genotype = load_genotype("source/genotypes/chr@", reference, allow_indel_swaps=True)
trait_a.save_cache(cache_root / "trait_a.sumstats.npz")
genotype.save_cache(cache_root / "genotype.npz")
```

See the [allele matching guide](ALLELE_MATCHING.md) for when to enable this
option and [cache migration](ALLELE_MATCHING.md#genotype-caches-and-bed-overrides) for
existing genotype caches.

The original BED files must still be available when later fetching genotype
calls, unless a replacement `bed_path` is passed to `fetch_genotypes`.

## Build MATLAB/Octave Caches

This pass writes MATLAB/Octave-native caches for the same reference-aligned
objects. The LD panel itself is already stored as the converted `.mat`
distribution from the previous step.

```matlab
addpath('../../matlab')

mat_cache = 'derived/matlab_cache';
if exist(mat_cache, 'dir') ~= 7
    mkdir(mat_cache);
end

reference = statgen.load_reference( ...
    'source/ld_reference/chr@.bim', ...
    {'21', '22', 'X'});

trait_a = statgen.load_sumstats( ...
    'source/sumstats/trait_a.tsv.gz', ...
    reference);
trait_b = statgen.load_sumstats( ...
    'source/sumstats/trait_b.tsv.gz', ...
    reference);

annot_root = 'source/annotations';
bed_paths = {
    fullfile(annot_root, 'coding_exon.bed')
    fullfile(annot_root, 'exon.bed')
    fullfile(annot_root, 'intron.bed')
    fullfile(annot_root, 'utr3.bed')
    fullfile(annot_root, 'utr5.bed')
    fullfile(annot_root, 'whole_gene.bed')
};
binary_annotations = statgen.load_annotations(bed_paths, reference);
grouped_annotations = statgen.load_annotation( ...
    fullfile(annot_root, 'functional_groups.annot'), ...
    reference, ...
    'has_header', true, ...
    'group_column', 'group');
continuous_annotations = statgen.load_annotation( ...
    fullfile(annot_root, 'conservation.annot'), ...
    reference, ...
    'has_header', true, ...
    'value_columns', {'conservation', 'promoter_activity'}, ...
    'annotation_metadata_path', fullfile(annot_root, 'conservation.meta'));
annotations = binary_annotations.union_annotations(grouped_annotations).union_annotations(continuous_annotations);

genotype = statgen.load_genotype( ...
    'source/genotypes/chr@', ...
    reference);

assert(reference.is_object_compatible(trait_a))
assert(reference.is_object_compatible(trait_b))
assert(reference.is_object_compatible(annotations))
assert(reference.is_object_compatible(genotype))

reference.save_cache(fullfile(mat_cache, 'reference.mat'));
trait_a.save_cache(fullfile(mat_cache, 'trait_a.sumstats.mat'));
trait_b.save_cache(fullfile(mat_cache, 'trait_b.sumstats.mat'));
annotations.save_cache(fullfile(mat_cache, 'annotations.mat'));
genotype.save_cache(fullfile(mat_cache, 'genotype.mat'));
```

MATLAB/Octave uses the same default matching policy. To enable indel swaps,
pass `'allow_indel_swaps', true` to each applicable source loader before saving
its cache:

```matlab
trait_a = statgen.load_sumstats('source/sumstats/trait_a.tsv.gz', reference, ...
    'allow_indel_swaps', true);
genotype = statgen.load_genotype('source/genotypes/chr@', reference, ...
    'allow_indel_swaps', true);
trait_a.save_cache(fullfile(mat_cache, 'trait_a.sumstats.mat'));
genotype.save_cache(fullfile(mat_cache, 'genotype.mat'));
```

See the [allele matching guide](ALLELE_MATCHING.md) for details.

## Outputs

This tutorial prepares these outputs for the next tutorials:

```text
derived/ld_npz/
derived/ld_mat/
derived/python_cache/
derived/matlab_cache/
```

The `scratch/ld_build/` directory contains temporary LD build files and can be
deleted after successful builds.

Continue with:

- [Tutorial 2: Python Analysis from Caches](TUTORIAL_2_PYTHON.md)
- [Tutorial 3: MATLAB Analysis from Caches](TUTORIAL_3_MATLAB.md)
- [Tutorial 4: R Analysis](TUTORIAL_4_R.md)
