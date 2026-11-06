# Tutorial 3: MATLAB Analysis

This tutorial uses the MATLAB/Octave artifacts prepared by
[Tutorial 1](TUTORIAL_1_PREPARE_DATA.md) and runs common operations:

- check object compatibility;
- compute LD-weighted annotations with `LDPanel.multiply_r2(...)`;
- prune summary-statistics scores with `statgen.fast_prune(...)`;
- fetch genotype calls for pruned SNPs with `GenotypePanel.fetch_genotypes(...)`.

The example uses chromosomes 21, 22, and X. Paths are relative to the same
tutorial working directory used in Tutorial 1.

## Contents

- [Prerequisites](#prerequisites)
- [Prepared Inputs](#prepared-inputs)
- [Load Caches](#load-caches)
- [Check Compatibility](#check-compatibility)
- [Run LD Operations](#run-ld-operations)
- [Fetch Genotypes for Pruned SNPs](#fetch-genotypes-for-pruned-snps)
- [Subset to One Chromosome](#subset-to-one-chromosome)

## Prerequisites

This tutorial assumes Tutorial 1 has already prepared the reusable artifacts
under `derived/`:

```text
derived/ld_mat/
derived/matlab_cache/reference.mat
derived/matlab_cache/trait_a.sumstats.mat
derived/matlab_cache/trait_b.sumstats.mat
derived/matlab_cache/annotations.mat
derived/matlab_cache/genotype.mat
```

Use MATLAB or Octave with the repository's MATLAB package folder on the path:

```matlab
addpath('/path/to/statgen/matlab')
```

If you have not built the caches yet, start with
[Tutorial 1](TUTORIAL_1_PREPARE_DATA.md).

## Prepared Inputs

```matlab
mat_cache = 'derived/matlab_cache';
ld_mat_dir = 'derived/ld_mat';
```

## Load Caches

```matlab
reference = statgen.load_ld_reference(ld_mat_dir);

ld = statgen.load_ld(ld_mat_dir, 'female');

trait_a = statgen.load_sumstats_cache( ...
    fullfile(mat_cache, 'trait_a.sumstats.mat'));
trait_b = statgen.load_sumstats_cache( ...
    fullfile(mat_cache, 'trait_b.sumstats.mat'));
annotations = statgen.load_annotations_cache( ...
    fullfile(mat_cache, 'annotations.mat'));
genotype = statgen.load_genotype_cache( ...
    fullfile(mat_cache, 'genotype.mat'));
```

Use `'male'` as the second `load_ld` argument when an analysis should use the
male chrX LD shard by default.

## Check Compatibility

All loaded objects should use the same reference coordinates.

```matlab
assert(reference.is_object_compatible(ld))
assert(reference.is_object_compatible(trait_a))
assert(reference.is_object_compatible(trait_b))
assert(reference.is_object_compatible(annotations))
assert(reference.is_object_compatible(genotype))
```

## Run LD Operations

`annotations.annomat` is a sparse numeric matrix aligned to the reference, so
it can be multiplied by LD `r²` directly. Use `annotations.is_binary` to
distinguish binary BED-derived columns from continuous annotation columns.

```matlab
a1freq = ld.a1freq();
annomat = annotations.annomat;
continuous_names = annotations.annonames(~annotations.is_binary);
continuous_meta = annotations.annotation_metadata(~annotations.is_binary);

ld_weighted_annotations = ld.multiply_r2(annomat);
trait_a_pruned = statgen.fast_prune(trait_a.logpvec, ld, 0.2);
```

`trait_a_pruned` is still in reference coordinates. Pruned variants are set to
`NaN`.

## Fetch Genotypes for Pruned SNPs

Fetch genotype calls for retained pruned SNPs that are also present in the
genotype source.

```matlab
pruned_snp_indices = find(isfinite(trait_a_pruned) & genotype.is_present);
pruned_genotypes = genotype.fetch_genotypes(pruned_snp_indices);

disp(size(pruned_snp_indices))
disp(size(pruned_genotypes))
disp(genotype.is_allele_swapped(pruned_snp_indices))
```

`pruned_snp_indices` are one-based reference-coordinate indices in
MATLAB/Octave. `pruned_genotypes` has samples as rows and requested SNPs as
columns.
Calls use the raw diploid scale by default; pass `'ploidy_scaled'` to adjust
for biological ploidy.

Counts, Z, beta, and EAF describe the reference panel's `a1` (ALT) allele;
see the [allele matching guide](ALLELE_MATCHING.md).

When using a genotype metadata cache, the original BED files must still be
available. If they moved, pass a replacement BED path as the second argument to
`fetch_genotypes`.

## Subset to One Chromosome

For analyses restricted to a chromosome subset, subset the reference first and
then load or select compatible objects against that subset.

```matlab
ref_chr21 = reference.select_shards({'21'});
ld_chr21 = statgen.load_ld(ld_mat_dir, {'21'});
trait_a_chr21 = trait_a.select_shards({'21'});
ann_chr21 = annotations.select_shards({'21'});
genotype_chr21 = genotype.select_shards({'21'});
```
