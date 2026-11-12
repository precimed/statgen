# Tutorial 2: Python Analysis

This tutorial uses the Python artifacts prepared by
[Tutorial 1](TUTORIAL_1_PREPARE_DATA.md) and runs common operations:

- check object compatibility;
- compute LD-weighted annotations with `LDPanel.multiply_r2(...)`;
- prune summary-statistics scores with `fast_prune(...)`;
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
derived/ld_npz/
derived/python_cache/reference.npz
derived/python_cache/trait_a.sumstats.npz
derived/python_cache/trait_b.sumstats.npz
derived/python_cache/annotations.npz
derived/python_cache/genotype.npz
```

Use a Python environment with `statgen` installed:

```bash
conda activate statgen
```

If you have not created the environment yet, see the setup commands in
[Tutorial 1](TUTORIAL_1_PREPARE_DATA.md#prerequisites).

## Prepared Inputs

```python
from pathlib import Path

import numpy as np

from statgen.annotations import load_annotations_cache
from statgen.genotype import load_genotype_cache
from statgen.ld import fast_prune, load_ld, load_ld_reference
from statgen.sumstats import load_sumstats_cache

cache_root = Path("derived/python_cache")
ld_npz_dir = Path("derived/ld_npz")
```

## Load Caches

```python
reference = load_ld_reference(ld_npz_dir)

ld = load_ld(
    ld_npz_dir,
    default_chrX_sex="female",
)

trait_a = load_sumstats_cache(cache_root / "trait_a.sumstats.npz")
trait_b = load_sumstats_cache(cache_root / "trait_b.sumstats.npz")
annotations = load_annotations_cache(cache_root / "annotations.npz")
genotype = load_genotype_cache(cache_root / "genotype.npz")
```

Use `default_chrX_sex="male"` when an analysis should use the male chrX LD
shard by default.

## Check Compatibility

All loaded objects should use the same reference coordinates.

```python
assert reference.is_object_compatible(ld)
assert reference.is_object_compatible(trait_a)
assert reference.is_object_compatible(trait_b)
assert reference.is_object_compatible(annotations)
assert reference.is_object_compatible(genotype)
```

## Run LD Operations

`annotations.annomat` is a sparse numeric matrix aligned to the reference, so
it can be multiplied by LD `r²` directly. Use `annotations.is_binary` to
distinguish binary BED-derived columns from continuous annotation columns.

```python
a1freq = ld.a1freq()
annomat = annotations.annomat
continuous_names = annotations.annonames[~annotations.is_binary]
continuous_meta = annotations.annotation_metadata[~annotations.is_binary]

ld_weighted_annotations = ld.multiply_r2(annomat)
trait_a_pruned = fast_prune(trait_a.logpvec, ld, r2_threshold=0.2)
```

`trait_a_pruned` is still in reference coordinates. Pruned variants are set to
`NaN`.

## Fetch Genotypes for Pruned SNPs

Fetch genotype calls for retained pruned SNPs that are also present in the
genotype source.

```python
pruned_snp_indices = np.flatnonzero(
    np.isfinite(trait_a_pruned) & genotype.is_present
)
pruned_genotypes = genotype.fetch_genotypes(pruned_snp_indices)

print(pruned_snp_indices.shape)
print(pruned_genotypes.shape)
print(genotype.is_allele_swapped[pruned_snp_indices])
```

`pruned_snp_indices` are zero-based reference-coordinate indices in Python.
`pruned_genotypes` has samples as rows and requested SNPs as columns.
Calls use the raw diploid scale by default; pass
`haploid_mode="ploidy_scaled"` to adjust for biological ploidy.

Counts, Z, beta, and EAF describe the reference panel's `a1` (ALT) allele;
see the [allele matching guide](ALLELE_MATCHING.md).

When using a genotype metadata cache, BED files must be available to fetch calls.
If they moved, pass a replacement `bed_path` to `load_genotype_cache`.

## Subset to One Chromosome

For analyses restricted to a chromosome subset, subset the reference first and
then load or select compatible objects against that subset.

```python
ref_chr21 = reference.select_shards(["21"])
ld_chr21 = load_ld(ld_npz_dir, shards=["21"])
trait_a_chr21 = trait_a.select_shards(["21"])
ann_chr21 = annotations.select_shards(["21"])
genotype_chr21 = genotype.select_shards(["21"])
```
