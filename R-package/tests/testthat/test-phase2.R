.phase2_reference_template <- function() {
  ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen", mustWork = TRUE)
  file.path(dirname(ref_path), "reference_chr@.bim")
}

.phase2_capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

.phase2_write_lf <- function(text, path) {
  writeBin(charToRaw(text), path)
}

test_that("missing p keeps optional-column presence when all source rows drop", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chr\tbp\ta1\ta2\tp\tz\tn\tbeta\tse\teaf\tinfo",
    "1\t100\tA\tG\tNA\t1\t100\t0.1\t0.2\t0.3\t0.9"
  ), path)
  captured <- .phase2_capture_warnings(load_sumstats(path, ref))
  s <- captured$value
  expect_true(any(grepl("no matching source rows", captured$warnings)))
  expect_false(any(grepl("is absent", captured$warnings)))
  expect_false(any(is_present(s)))
  for (accessor in list(logpvec, zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)) {
    expect_length(accessor(s), num_snp(ref))
    expect_true(all(is.nan(accessor(s))))
  }
  cache <- tempfile(fileext = ".rds")
  save_sumstats_cache(s, cache)
  loaded <- load_sumstats_cache(cache)
  for (accessor in list(logpvec, zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)) {
    expect_identical(accessor(loaded), accessor(s))
  }

  writeLines(c("chr\tbp\ta1\ta2\tp", "1\t100\tA\tG\tNaN"), path)
  s <- suppressWarnings(load_sumstats(path, ref))
  for (accessor in list(zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)) {
    expect_null(accessor(s))
  }
})

test_that("factory NA and NaN p mask optional values without changing inputs", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  p <- c(NA_real_, NaN, 0, 1, 0.1)
  values <- as.numeric(seq_len(num_snp(ref)))
  s <- create_sumstats(ref, p, values, values, values, values, values, values)
  expect_identical(is_present(s), c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_true(all(is.nan(logpvec(s)[1:2])))
  for (accessor in list(zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)) {
    expect_identical(accessor(s), c(NaN, NaN, 3, 4, 5))
  }
  expect_identical(values, as.numeric(1:5))
  expect_identical(p, c(NA_real_, NaN, 0, 1, 0.1))
  cache <- tempfile(fileext = ".rds")
  save_sumstats_cache(s, cache)
  loaded <- load_sumstats_cache(cache)
  for (accessor in list(logpvec, zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)) {
    expect_identical(accessor(loaded), accessor(s))
  }
  expect_null(zvec(create_sumstats(ref, p)))
  expect_error(create_sumstats(ref, c("bad", rep("0.1", 4))), "must be numeric")
  for (value in c(Inf, -Inf)) {
    expect_error(create_sumstats(ref, c(value, rep(0.1, 4))), "finite numeric")
  }
})

test_that("numeric errors report the original row and value without fread warnings", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  columns <- c("chr", "bp", "a1", "a2", "p", "z", "n", "beta", "se", "eaf", "info")
  for (field in c("bp", "p", "z", "n", "beta", "se", "eaf", "info")) {
    path <- tempfile(fileext = ".tsv")
    row <- c("1", "100", "A", "G", "NA", "1", "2", "3", "4", "0.5", "6")
    row[match(field, columns)] <- "abc"
    writeLines(c(
      paste(columns, collapse = "\t"),
      "Y\t100\tA\tG\tNA\t1\t2\t3\t4\t0.5\t6",
      paste(row, collapse = "\t"),
      paste(replace(row, row == "abc", "def"), collapse = "\t")
    ), path)
    expect_no_warning(msg <- tryCatch(load_sumstats(path, ref), error = conditionMessage))
    expect_match(msg, "row 3:", fixed = TRUE)
    expect_match(msg, sQuote("abc"), fixed = TRUE)
    expect_false(grepl("def", msg, fixed = TRUE))
    if (field == "bp") {
      expect_match(msg, "bp is not an integer:", fixed = TRUE)
      expect_false(grepl("or missing", msg, fixed = TRUE))
    }
  }
})

test_that("load_sumstats rejects files where fread stops at an extra field", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chr\tbp\ta1\ta2\tp",
    "1\t100\tA\tG\t0.1",
    "1\t200\tC\tT\t0.2\textra",
    "1\t300\tA\tC\t0.3"
  ), path)
  for (attempt in seq_len(2L)) {
    expect_no_warning(expect_error(load_sumstats(path, ref), "malformed TSV"))
  }
})

test_that("load_sumstats aligns TSV rows and cache round-trips", {
  sum_path <- system.file("extdata", "traits.tsv.gz", package = "statgen", mustWork = TRUE)
  ref <- load_reference(.phase2_reference_template())

  expect_warning(s <- load_sumstats(sum_path, ref), "zvec has missing values")
  expect_s3_class(s, "Sumstats")
  expect_equal(num_snp(s), 8L)
  expect_equal(is_present(s), c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(zvec(s), c(2.5, NaN, 1.8, -1.2, NaN, 3.0, 0.5, NaN))
  expect_equal(nvec(s), c(1000, 1000, 1000, 1000, NaN, 500, 500, NaN))
  expect_true(is.infinite(logpvec(s)[[3]]))
  expect_true(is.null(beta_vec(s)))

  cache <- tempfile(fileext = ".rds")
  save_sumstats_cache(s, cache)
  expect_warning(loaded <- load_sumstats_cache(cache), "zvec has missing values")
  expect_equal(logpvec(loaded), logpvec(s))
  expect_equal(nvec(loaded), nvec(s))
  expect_equal(is_present(loaded), is_present(s))

  x_only <- select_shards(loaded, "X")
  expect_equal(num_snp(x_only), 3L)
  expect_equal(is_present(x_only), c(TRUE, TRUE, FALSE))
})

test_that("load_sumstats accepts documented column aliases", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tN",
    "1\t100\trs1001\tA\tG\t0\t2\t100"
  ), path)

  captured <- .phase2_capture_warnings(load_sumstats(path, ref))
  s <- captured$value
  expect_true(is.infinite(logpvec(s)[[1]]))
  expect_equal(zvec(s)[[1]], 2)
  expect_equal(nvec(s)[[1]], 100)
  expect_true(is_present(s)[[1]])
})

test_that("create_sumstats preserves optional fields and save_cache convenience", {
  ref <- load_reference(.phase2_reference_template())
  n <- num_snp(ref)
  s <- create_sumstats(
    ref,
    p = rep(0.5, n),
    z = seq_len(n),
    n = rep(100, n),
    beta = rep(0.1, n),
    se = rep(0.2, n),
    eaf = rep(0.3, n),
    info = rep(0.9, n)
  )

  expect_equal(nvec(s), rep(100, n))
  expect_equal(beta_vec(s), rep(0.1, n))
  expect_equal(se_vec(s), rep(0.2, n))
  expect_equal(eaf_vec(s), rep(0.3, n))
  expect_equal(info_vec(s), rep(0.9, n))

  cache <- tempfile(fileext = ".rds")
  save_cache(s, cache)
  loaded <- load_sumstats_cache(cache, shards = "X")
  expect_equal(num_snp(loaded), 3L)
  expect_equal(nvec(loaded), rep(100, 3L))
})

test_that("load_sumstats warns for absent zvec and nvec", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chr\tbp\ta1\ta2\tp",
    "1\t100\tA\tG\t0.1"
  ), path)

  captured <- .phase2_capture_warnings(load_sumstats(path, ref))
  expect_true(any(grepl("zvec is absent", captured$warnings)))
  expect_true(any(grepl("nvec is absent", captured$warnings)))
  expect_null(zvec(captured$value))
  expect_null(nvec(captured$value))

  cache <- tempfile(fileext = ".rds")
  save_sumstats_cache(captured$value, cache)
  loaded <- .phase2_capture_warnings(load_sumstats_cache(cache))
  expect_true(any(grepl("zvec is absent", loaded$warnings)))
  expect_true(any(grepl("nvec is absent", loaded$warnings)))
})

test_that("load_annotations paints BED intervals and cache round-trips", {
  bed1 <- system.file("extdata", "anno1.bed", package = "statgen", mustWork = TRUE)
  bed2 <- system.file("extdata", "anno2.bed", package = "statgen", mustWork = TRUE)
  ref <- load_reference(.phase2_reference_template())
  ann <- load_annotations(c(bed1, bed2), ref)

  expect_s3_class(ann, "AnnotationPanel")
  expect_equal(num_snp(ann), 8L)
  expect_equal(num_annot(ann), 2L)
  expect_equal(annonames(ann), c("anno1", "anno2"))
  expect_s4_class(annomat(ann), "dgCMatrix")
  expect_equal(is_binary(ann), c(TRUE, TRUE))
  expect_true(all(grepl("source_file", annotation_metadata(ann), fixed = TRUE)))
  ann_dense <- as.matrix(annomat(ann))
  expect_equal(colnames(annomat(ann)), c("anno1", "anno2"))
  expect_equal(as.vector(ann_dense[, 1L]), c(1, 1, 1, 1, 1, 0, 0, 0))
  expect_equal(as.vector(ann_dense[, 2L]), c(0, 0, 0, 0, 0, 1, 1, 1))

  cache <- tempfile(fileext = ".rds")
  save_cache(ann, cache)
  loaded <- load_annotations_cache(cache, shards = "X")
  expect_equal(annonames(loaded), c("anno1", "anno2"))
  expect_equal(colnames(annomat(loaded)), c("anno1", "anno2"))
  expect_equal(as.matrix(annomat(loaded)), as.matrix(annomat(ann))[6:8, ])
})

test_that("load_annotation paints continuous values and preserves metadata", {
  ref <- load_reference(.phase2_reference_template())
  annot <- tempfile(fileext = ".annot")
  writeLines(c(
    "chrom\tstart0\tend0\tscore\tweight",
    "1\t99\t200\t0.5\t10",
    "1\t299\t400\t1.5\t20",
    "X\t99\t301\t2.5\t30"
  ), annot)
  sidecar <- tempfile(fileext = ".meta")
  writeLines(c("chrom meta", "start meta", "end meta", "score meta", "weight meta"), sidecar)

  ann <- load_annotation(
    annot,
    ref,
    has_header = TRUE,
    value_columns = c("weight", "score"),
    annotation_metadata_path = sidecar
  )
  expect_equal(annonames(ann), c("weight", "score"))
  expect_equal(is_binary(ann), c(FALSE, FALSE))
  expect_equal(annotation_metadata(ann), c("weight meta", "score meta"))
  expect_s4_class(annomat(ann), "dgCMatrix")
  expect_equal(
    as.vector(as.matrix(annomat(ann))[, 1L]),
    c(10, 10, 20, 20, 0, 30, 30, 30)
  )
  expect_equal(
    as.vector(as.matrix(annomat(ann))[, 2L]),
    c(0.5, 0.5, 1.5, 1.5, 0, 2.5, 2.5, 2.5)
  )

  cache <- tempfile(fileext = ".rds")
  save_annotations_cache(ann, cache)
  loaded <- load_annotations_cache(cache)
  expect_equal(is_binary(loaded), c(FALSE, FALSE))
  expect_equal(annotation_metadata(loaded), c("weight meta", "score meta"))
  expect_equal(as.matrix(annomat(loaded)), as.matrix(annomat(ann)))
})

test_that("load_annotation supports grouped binary annotations", {
  ref <- load_reference(.phase2_reference_template())
  grouped <- system.file("extdata", "grouped.annot", package = "statgen", mustWork = TRUE)
  ann <- load_annotation(grouped, ref, has_header = TRUE, group_column = "group")
  indexed <- load_annotation(grouped, ref, has_header = TRUE, group_column = 4L)

  expect_equal(annonames(ann), c("coding", "regulatory"))
  expect_equal(is_binary(ann), c(TRUE, TRUE))
  expect_equal(
    unname(as.matrix(annomat(ann))),
    matrix(
      c(1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1),
      nrow = 8L,
      ncol = 2L
    )
  )
  meta <- lapply(annotation_metadata(ann), jsonlite::fromJSON)
  expect_equal(vapply(meta, function(x) x$group_value, character(1)), c("coding", "regulatory"))
  expect_equal(vapply(meta, function(x) x$num_source_intervals, integer(1)), c(2L, 2L))
  indexed_meta <- lapply(annotation_metadata(indexed), jsonlite::fromJSON)
  expect_equal(annonames(indexed), c("coding", "regulatory"))
  expect_equal(vapply(indexed_meta, function(x) as.integer(x$group_column), integer(1)), c(4L, 4L))

  overlap <- tempfile(fileext = ".annot")
  writeLines(c(
    "chrom\tstart0\tend0\tgroup",
    "1\t99\t250\tcoding",
    "1\t150\t401\tcoding",
    "X\t99\t301\tregulatory"
  ), overlap)
  overlapped <- load_annotation(overlap, ref, has_header = TRUE, group_column = "group")
  expect_equal(
    unname(as.matrix(annomat(overlapped))),
    matrix(
      c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1),
      nrow = 8L,
      ncol = 2L
    )
  )
  overlap_meta <- lapply(annotation_metadata(overlapped), jsonlite::fromJSON)
  expect_equal(vapply(overlap_meta, function(x) x$num_source_intervals, integer(1)), c(2L, 1L))

  conflicts <- tempfile(fileext = ".annot")
  writeLines(c(
    "chrom\tstart0\tend0\tgroup\tscore",
    "1\t99\t200\tcoding\t1"
  ), conflicts)
  sidecar <- tempfile(fileext = ".meta")
  .phase2_write_lf("meta\n", sidecar)
  expect_error(load_annotation(conflicts, ref, has_header = TRUE, group_column = "group", value_columns = "score"), "mutually exclusive")
  expect_error(load_annotation(conflicts, ref, has_header = TRUE, group_column = "group", annotation_names = "coding"), "annotation_names is invalid")
  expect_error(load_annotation(conflicts, ref, has_header = TRUE, group_column = "group", annotation_metadata = "meta"), "annotation_metadata is invalid")
  expect_error(load_annotation(conflicts, ref, has_header = TRUE, group_column = "group", annotation_metadata_path = sidecar), "annotation_metadata_path is invalid")
})

test_that("load_annotation supports binary sidecar, name override, and headerless 4-column default", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  bed <- tempfile(fileext = ".bed")
  writeLines("1\t99\t200", bed)
  sidecar <- tempfile(fileext = ".meta")
  .phase2_write_lf("binary\nmetadata\n", sidecar)

  binary <- load_annotation(
    bed,
    ref,
    annotation_names = "renamed_binary",
    annotation_metadata_path = sidecar
  )
  expect_equal(annonames(binary), "renamed_binary")
  expect_equal(is_binary(binary), TRUE)
  expect_equal(annotation_metadata(binary), "binary\nmetadata\n")
  expect_equal(as.vector(as.matrix(annomat(binary))), c(1, 1, 0, 0, 0))

  annot <- tempfile(fileext = ".annot")
  writeLines(c("1\t99\t200\t0.5", "1\t299\t400\t1.5"), annot)
  continuous <- load_annotation(annot, ref)
  expect_equal(annonames(continuous), tools::file_path_sans_ext(basename(annot)))
  expect_equal(is_binary(continuous), FALSE)
  expect_equal(as.vector(as.matrix(annomat(continuous))), c(0.5, 0.5, 1.5, 1.5, 0))
})

test_that("load_annotations supports batch metadata sidecars", {
  bed1 <- system.file("extdata", "anno1.bed", package = "statgen", mustWork = TRUE)
  bed2 <- system.file("extdata", "anno2.bed", package = "statgen", mustWork = TRUE)
  ref <- load_reference(.phase2_reference_template())
  sidecar <- tempfile(fileext = ".meta")
  .phase2_write_lf("sidecar\nmetadata\n", sidecar)

  ann <- load_annotations(
    c(bed1, bed2),
    ref,
    annotation_metadata_paths = c(sidecar, NA_character_)
  )
  expect_equal(annotation_metadata(ann)[[1]], "sidecar\nmetadata\n")
  expect_true(grepl("source_file", annotation_metadata(ann)[[2]], fixed = TRUE))
  expect_true(grepl("anno2.bed", annotation_metadata(ann)[[2]], fixed = TRUE))
})

test_that("load_annotations accepts scalar path and leading comments", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  bed <- tempfile(fileext = ".bed")
  writeLines(c("# leading comment", "", "1\t99\t200"), bed)

  ann <- load_annotations(bed, ref)
  expect_equal(num_annot(ann), 1L)
  expect_equal(as.vector(as.matrix(annomat(ann))), c(1, 1, 0, 0, 0))
})

test_that("annotation factories, selection, and union work", {
  ref <- load_reference(.phase2_reference_template())
  n <- num_snp(ref)
  a <- create_annotation(ref, rep(c(1, 0), length.out = n), "a", annotation_metadata = "meta a")
  b <- create_annotations(
    ref,
    annotation_matrix = matrix(rep(c(0, 1), length.out = n), ncol = 1L),
    annotation_names = "b",
    annotation_metadata = "meta b"
  )
  u <- union_annotations(a, b)

  expect_equal(annonames(u), c("a", "b"))
  expect_equal(annotation_metadata(u), c("meta a", "meta b"))
  expect_equal(is_binary(u), c(TRUE, TRUE))
  expect_equal(colnames(annomat(u)), c("a", "b"))
  selected <- select_annotations(u, names = c("b", "a"))
  expect_equal(annonames(selected), c("b", "a"))
  expect_equal(annotation_metadata(selected), c("meta b", "meta a"))
  expect_equal(colnames(annomat(selected)), c("b", "a"))
  expect_error(select_annotations(u, "missing"), "unknown annotation")
  expect_error(select_annotations(u, c("a", "a")), "names must be unique")

  cont <- create_annotations(ref, 2 * matrix(1, nrow = n, ncol = 1L), "cont", is_binary = FALSE)
  expect_equal(is_binary(cont), FALSE)
  expect_error(create_annotations(ref, 2 * matrix(1, nrow = n, ncol = 1L), "bad", is_binary = TRUE), "non-binary")
})

test_that("annotations cache validation rejects bad schema and zero shards", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  ann <- create_annotation(ref, rep(TRUE, num_snp(ref)), "a")
  cache <- tempfile(fileext = ".rds")
  save_annotations_cache(ann, cache)
  payload <- readRDS(cache)

  bad_schema <- tempfile(fileext = ".rds")
  payload$metadata$schema <- "annotations_cache/bad"
  saveRDS(payload, bad_schema)
  expect_error(load_annotations_cache(bad_schema), "Unsupported annotations cache schema")

  zero_shards <- tempfile(fileext = ".rds")
  saveRDS(
    list(
      metadata = list(
        schema = "annotations_cache/0.1",
        n_shards = 0L,
        shard_labels = character(),
        shard_checksums = character(),
        shard_start0 = numeric(),
        shard_stop0 = numeric()
      ),
      annomat = Matrix::Matrix(matrix(logical(), nrow = 0L, ncol = 1L), sparse = TRUE),
      annonames = "a"
    ),
    zero_shards
  )
  expect_error(load_annotations_cache(zero_shards), "n_shards must be at least 1")
})
