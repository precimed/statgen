test_that("allele swaps agree with the portable independent expectations", {
  root <- system.file("extdata", "allele_swaps", package = "statgen", mustWork = TRUE)
  oracle <- jsonlite::fromJSON(file.path(root, "expected.json"), simplifyMatrix = TRUE)
  ref <- load_reference(file.path(root, "reference.bim"))
  previous <- get_verbosity()
  on.exit(set_verbosity(previous))
  set_verbosity("quiet")
  # JSON nulls represent missing numeric values; normalize them to NaN.
  numeric_values <- function(x) {
    if (!is.list(x)) {
      x[is.na(x)] <- NaN
      return(x)
    }
    vapply(x, function(v) if (is.null(v)) NaN else as.numeric(v), numeric(1))
  }
  numeric_matrix <- function(x) {
    if (is.matrix(x)) {
      x[is.na(x)] <- NaN
      return(x)
    }
    do.call(rbind, lapply(x, numeric_values))
  }
  for (enabled in c(FALSE, TRUE)) {
    expected <- oracle$options[[tolower(as.character(enabled))]]
    g <- load_genotype(file.path(root, "genotype", "@"), ref, allow_indel_swaps = enabled)
    s <- suppressWarnings(load_sumstats(file.path(root, "sumstats.tsv.gz"), ref, allow_indel_swaps = enabled))
    expect_identical(source_row0(g), expected$source_row0)
    expect_identical(is_allele_swapped(g), expected$is_allele_swapped)
    expect_identical(is_present(g), expected$source_row0 >= 0L)
    indices <- c(rev(which(is_present(g))), 2L, 14L, 2L)
    raw <- numeric_matrix(expected$raw_counts)[, indices, drop = FALSE]
    raw[is.na(raw)] <- -1L
    storage.mode(raw) <- "integer"
    expect_identical(fetch_genotypes_int8(g, indices), raw)
    expect_equal(fetch_genotypes(g, indices, haploid_mode = "ploidy_scaled"),
                 numeric_matrix(expected$scaled_counts)[, indices, drop = FALSE])
    for (field in c("z", "n", "beta", "se", "eaf", "info")) {
      accessor <- switch(field, z = zvec, n = nvec, beta = beta_vec,
                         se = se_vec, eaf = eaf_vec, info = info_vec)
      expect_equal(accessor(s), numeric_values(expected$aligned_sumstats[[field]]))
    }
    cache <- tempfile(fileext = ".rds")
    on.exit(unlink(cache), add = TRUE)
    save_genotype_cache(g, cache)
    expect_identical(readRDS(cache)$metadata$schema, "genotype_cache/0.2")
    h <- load_genotype_cache(cache)
    expect_identical(is_allele_swapped(h), is_allele_swapped(g))
    expect_identical(fetch_genotypes_int8(h, indices), raw)
    hx <- load_genotype_cache(cache, "X")
    expect_identical(is_allele_swapped(hx), is_allele_swapped(g)[14:16])
  }
})
