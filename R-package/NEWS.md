# statgen 0.4.0

* Summary-statistics and genotype loaders now match reversed SNVs, including A/T and C/G under the same-strand assumption. Set `allow_indel_swaps = TRUE` to also match reversed indels with exactly one single-base allele.
* Matched allele swaps orient Z scores, BETA, EAF, and genotype calls to the reference `a1` (ALT) allele. Genotype panels expose `is_allele_swapped`; genotype caches now use schema `genotype_cache/0.2`.
* `load_sumstats()` now validates and drops rows with missing `p` before matching. `create_sumstats()` accepts missing `p`, setting corresponding log-p values and supplied optional statistics to missing.
* `load_reference()` now rejects duplicate variants and variants listed in both allele orders.
* Fixed variant matching for large allele hashes that could lose exact matches or silently accept duplicate keys.

# statgen 0.3.4

* `load_annotation()` gains a `group_column` argument to create one binary annotation per distinct group value in a source column.

# statgen 0.3.3

* New `load_annotation()` function loads a single BED-like annotation file, supporting multiple value columns, optional headers, annotation name overrides, and sidecar metadata files.
* `AnnotationPanel` now exposes `is_binary` and `annotation_metadata` accessors.
* Annotation matrices are now sparse numeric matrices (previously binary-only).

# statgen 0.3.2

* Initial CRAN submission.
