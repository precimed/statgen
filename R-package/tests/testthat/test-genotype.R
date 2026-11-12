.prepare_sharded_genotype_fixture <- function() {
  root <- tempfile("statgen-genotype-")
  dir.create(root)
  file.copy(system.file("extdata", "reference_chr1.bim", package = "statgen", mustWork = TRUE), file.path(root, "1.bim"))
  file.copy(system.file("extdata", "reference_chrX.bim", package = "statgen", mustWork = TRUE), file.path(root, "X.bim"))
  for (label in c("1", "X")) {
    for (suffix in c(".bed", ".bim", ".fam")) {
      file.copy(
        system.file("extdata", paste0("genotype_", label, suffix), package = "statgen", mustWork = TRUE),
        file.path(root, paste0("g", label, suffix))
      )
    }
  }
  file.copy(
    system.file("extdata", "genotype_X.ploidy", package = "statgen", mustWork = TRUE),
    file.path(root, "gX.ploidy")
  )
  root
}

test_that("load_genotype reads bundled PLINK fixture and fetches slices", {
  ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen", mustWork = TRUE)
  bed_path <- system.file("extdata", "genotype_1.bed", package = "statgen", mustWork = TRUE)
  prefix <- sub("\\.bed$", "", bed_path)

  ref <- load_reference(ref_path)
  g <- load_genotype(prefix, ref)

  expect_s3_class(g, "GenotypePanel")
  expect_equal(num_snp(g), 5L)
  expect_equal(num_sample(g), 4L)
  expect_equal(source_layout(g), "non_sharded")
  expect_equal(fid(g), c("FAM1", "FAM1", "FAM2", "FAM2"))
  expect_equal(sex(g), c(1L, 2L, 1L, 2L))
  expect_equal(source_row0(g), 0:4)

  gi <- fetch_genotypes_int8(g, c(1, 5))
  expect_equal(dim(gi), c(4L, 2L))
  expect_equal(as.integer(gi), rep(2L, 8L))
  gi_empty <- fetch_genotypes_int8(g, integer(0))
  expect_equal(dim(gi_empty), c(4L, 0L))
  expect_equal(typeof(gi_empty), "integer")

  gf <- fetch_genotypes(g, c(1, 5))
  expect_equal(dim(gf), c(4L, 2L))
  expect_equal(as.numeric(gf), rep(2, 8L))
})

test_that("sharded genotype exposes chrX ploidy, subject presence, and shard subset", {
  root <- .prepare_sharded_genotype_fixture()
  ref <- load_reference(file.path(root, "@.bim"))
  g <- load_genotype(file.path(root, "g@"), ref)

  expect_equal(source_layout(g), "sharded")
  expect_equal(ploidy_male(g), c(rep(2, 5), rep(1, 3)))
  expect_equal(ploidy_female(g), rep(2, 8))
  expect_equal(is_subject_present(g, "X"), rep(TRUE, 4))

  gx <- select_shards(g, "X")
  expect_equal(num_snp(gx), 3L)
  expect_equal(num_sample(gx), 4L)
  expect_equal(vapply(shards(gx), function(s) s$chr, character(1)), "X")
})

test_that("genotype metadata cache round-trips", {
  root <- .prepare_sharded_genotype_fixture()
  cache <- tempfile(fileext = ".rds")
  override <- file.path(root, "override@.bed")

  ref <- load_reference(file.path(root, "@.bim"))
  g <- load_genotype(file.path(root, "g@"), ref)
  save_genotype_cache(g, cache)
  file.copy(file.path(root, "g1.bed"), file.path(root, "override1.bed"))
  file.copy(file.path(root, "gX.bed"), file.path(root, "overrideX.bed"))
  unlink(file.path(root, "g1.bed"))
  unlink(file.path(root, "gX.bed"))
  expect_warning(loaded <- load_genotype_cache(cache), "Missing or incorrectly sized")
  loaded_x <- load_genotype_cache(cache, shards = "X", bed_path = override)
  relocated <- load_genotype_cache(cache, bed_path = override)

  expect_equal(source_layout(loaded), source_layout(g))
  expect_equal(source_row0(loaded), source_row0(g))
  expect_equal(fid(loaded), fid(g))
  expect_error(fetch_genotypes_int8(loaded, 1), "File not found")
  expect_equal(fetch_genotypes_int8(relocated, c(1, 6)), matrix(2L, 4, 2))
  expect_equal(num_snp(loaded_x), 3L)
  expect_equal(source_layout(loaded_x), "sharded")
  expect_equal(fetch_genotypes_int8(loaded_x, 1), fetch_genotypes_int8(relocated, 6))
})

test_that("genotype validation rejects bad sex and malformed cache", {
  root <- .prepare_sharded_genotype_fixture()
  writeLines("FAM1\tIND1\t0\t0\t3\t-9", file.path(root, "bad.fam"))
  file.copy(file.path(root, "g1.bim"), file.path(root, "bad.bim"))
  file.copy(file.path(root, "g1.bed"), file.path(root, "bad.bed"))
  ref1 <- load_reference(file.path(root, "1.bim"))

  expect_error(load_genotype(file.path(root, "bad"), ref1), "FAM sex")

  bad_cache <- tempfile(fileext = ".rds")
  saveRDS(list(metadata = list(schema = "genotype_cache/0.1", n_shards = 0L)), bad_cache)
  expect_error(load_genotype_cache(bad_cache), "n_shards")
})
