# Tutorial 4: R Analysis

This tutorial uses the same source fixtures and Python-built LD distribution as
[Tutorial 1](TUTORIAL_1_PREPARE_DATA.md), but runs the analysis workflow with
the R package.

R differs from Python and MATLAB/Octave in one important way:

```text
Python builds derived/ld_npz/
R prepares derived/ld_npz/ in place
MATLAB/Octave converts derived/ld_npz/ to derived/ld_mat/
```

R reads Python `.npz` LD shards directly after `prepare_ld_npz_for_r()` adds an
R reference-cache sidecar. With `extract_npz = TRUE`, R also creates sibling
`*.npz.d` extracted caches for faster repeated loads.

## Contents

- [Prerequisites](#prerequisites)
- [Prepare The LD Directory For R](#prepare-the-ld-directory-for-r)
- [Load Objects](#load-objects)
- [Check Compatibility](#check-compatibility)
- [Run LD Operations](#run-ld-operations)
- [Fetch Genotypes For Pruned SNPs](#fetch-genotypes-for-pruned-snps)
- [Save R Caches](#save-r-caches)
- [Subset To One Chromosome](#subset-to-one-chromosome)

## Prerequisites

Start from the Tutorial 1 working directory after Tutorial 1 has built at least
`derived/ld_npz/`.

If you cloned the repository and used the checked-in fixtures, that directory
is:

```bash
cd docs/tutorial_1_fixtures
```

If you are using GitHub release assets instead, extract
`statgen-tutorial-1-fixtures-<version>.zip`, run Tutorial 1 from the extracted
directory, then run this tutorial from that same directory.

Install the R package. After `statgen` is available on CRAN, use:

```r
install.packages("statgen")
```

When working from a cloned repository, install the local source tree instead:

```r
install.packages("../../R-package", repos = NULL, type = "source")
```

When working from GitHub release assets before the CRAN package is available,
install the attached R source package tarball:

```r
install.packages("/path/to/statgen_<version>.RENAME_TO_CRAN_FILENAME.tar.gz", repos = NULL, type = "source")
```

The `RENAME_TO_CRAN_FILENAME` marker is only for CRAN submission safety; local
installation reads the package metadata from inside the tarball.

Then load it:

```r
library(statgen)
set_verbosity("quiet")
```

## Prepare The LD Directory For R

Prepare the Python LD distribution in place:

```r
ld_npz_dir <- "derived/ld_npz"

prepare_ld_npz_for_r(ld_npz_dir, extract_npz = TRUE)
validate_ld_distribution(ld_npz_dir)
```

The manifest still names `.npz` LD shards. R adds `reference_cache.rds` and, when
requested, extracted `*.npz.d` cache directories. If extracted caches are absent
or stale, `load_ld()` falls back to reading the `.npz` shards normally.

## Load Objects

Use the LD distribution reference as the analysis reference.

```r
reference <- load_ld_reference(ld_npz_dir)

ld <- load_ld(
  ld_npz_dir,
  default_chrX_sex = "female"
)
```

Load summary statistics, annotations, and genotype metadata against the same
reference.

```r
trait_a <- load_sumstats("source/sumstats/trait_a.tsv.gz", reference)
trait_b <- load_sumstats("source/sumstats/trait_b.tsv.gz", reference)

annotation_paths <- file.path(
  "source/annotations",
  c(
    "coding_exon.bed",
    "exon.bed",
    "intron.bed",
    "utr3.bed",
    "utr5.bed",
    "whole_gene.bed"
  )
)
binary_annotations <- load_annotations(annotation_paths, reference)
grouped_annotations <- load_annotation(
  "source/annotations/functional_groups.annot",
  reference,
  has_header = TRUE,
  group_column = "group"
)
continuous_annotations <- load_annotation(
  "source/annotations/conservation.annot",
  reference,
  has_header = TRUE,
  value_columns = c("conservation", "promoter_activity"),
  annotation_metadata_path = "source/annotations/conservation.meta"
)
annotations <- union_annotations(
  union_annotations(binary_annotations, grouped_annotations),
  continuous_annotations
)

genotype <- load_genotype("source/genotypes/chr@", reference)
```

The loaders above use `match_by = "genomic_position"` by default. To enable
indel swaps when loading your data:

```r
trait_a <- load_sumstats("source/sumstats/trait_a.tsv.gz", reference,
                         allow_indel_swaps = TRUE)
genotype <- load_genotype("source/genotypes/chr@", reference,
                          allow_indel_swaps = TRUE)
```

For marker-ID and allele matching, use `match_by = "id"` with indel swaps
disabled:

```r
trait_a <- load_sumstats("source/sumstats/trait_a.tsv.gz", reference, match_by = "id")
genotype <- load_genotype("source/genotypes/chr@", reference, match_by = "id")
```

ID mode accepts non-palindromic SNVs and matches IDs case-sensitively. Sumstats
require `SNP`, the allele columns, and `P`, ignore coordinates, and search all
source rows. Genotypes search only the corresponding chromosome.

See the [allele matching guide](ALLELE_MATCHING.md), or `?allele_matching` in R,
for the loading options and how results are adjusted.

Use `default_chrX_sex = "male"` when the analysis should use male chrX LD by
default.

## Check Compatibility

All objects should share the same reference coordinates.

```r
stopifnot(is_object_compatible(reference, ld))
stopifnot(is_object_compatible(reference, trait_a))
stopifnot(is_object_compatible(reference, trait_b))
stopifnot(is_object_compatible(reference, annotations))
stopifnot(is_object_compatible(reference, genotype))
```

## Run LD Operations

`annomat(annotations)` is a sparse numeric matrix aligned to the reference, so
it can be multiplied by LD `r^2` directly. Use `is_binary(annotations)` to
distinguish binary membership columns from continuous numeric annotations.

```r
a1 <- a1freq(ld)
continuous_names <- annonames(annotations)[!is_binary(annotations)]
continuous_meta <- annotation_metadata(annotations)[!is_binary(annotations)]
ld_weighted_annotations <- multiply_r2(ld, annomat(annotations))
trait_a_pruned <- fast_prune(logpvec(trait_a), ld, r2_threshold = 0.2)
```

`trait_a_pruned` remains in reference coordinates. Pruned variants are set to
`NaN`.

## Fetch Genotypes For Pruned SNPs

Fetch genotype calls for retained pruned SNPs that are also present in the
genotype source.

```r
pruned_snp_indices <- which(is.finite(trait_a_pruned) & is_present(genotype))
pruned_genotypes <- fetch_genotypes(genotype, pruned_snp_indices)

dim(pruned_genotypes)
is_allele_swapped(genotype)[pruned_snp_indices]
```

R SNP indices are one-based. Genotype matrices have samples as rows and
requested SNPs as columns. Calls use the raw diploid scale by default; pass
`haploid_mode = "ploidy_scaled"` to adjust for biological ploidy.

When using a genotype metadata cache, the original BED files must still be
available. If they moved, pass a replacement `bed_path` to the fetch function.

## Save R Caches

R caches are runtime-native RDS files. They are useful when the same source
files are loaded repeatedly in R.

```r
r_cache <- "derived/r_cache"
dir.create(r_cache, recursive = TRUE, showWarnings = FALSE)

save_reference_cache(reference, file.path(r_cache, "reference.rds"))
save_sumstats_cache(trait_a, file.path(r_cache, "trait_a.sumstats.rds"))
save_sumstats_cache(trait_b, file.path(r_cache, "trait_b.sumstats.rds"))
save_annotations_cache(annotations, file.path(r_cache, "annotations.rds"))
save_genotype_cache(genotype, file.path(r_cache, "genotype.rds"))
```

Load the caches later with the corresponding `load_*_cache()` functions.

## Subset To One Chromosome

For analyses restricted to a chromosome subset, load or select compatible
objects against that subset.

```r
ref_chr21 <- select_shards(reference, "21")
ld_chr21 <- load_ld(ld_npz_dir, shards = "21")
trait_a_chr21 <- select_shards(trait_a, "21")
ann_chr21 <- select_shards(annotations, "21")
genotype_chr21 <- select_shards(genotype, "21")
```
