.parse_bim <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  df <- tryCatch(
    utils::read.table(
      path,
      header = FALSE,
      sep = "",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      colClasses = rep("character", 6L),
      fill = FALSE
    ),
    error = function(e) {
      stop(sprintf("%s: expected 6 whitespace-delimited columns", path), call. = FALSE)
    }
  )
  if (ncol(df) != 6L) {
    stop(sprintf("%s: expected 6 whitespace-delimited columns, got %d", path, ncol(df)), call. = FALSE)
  }
  names(df) <- c("chr", "snp", "cm", "bp", "a1", "a2")
  bad_chr <- is.na(df$chr) | df$chr == ""
  if (any(bad_chr)) {
    stop(sprintf("%s:%d: chr must be non-empty", path, which(bad_chr)[[1]]), call. = FALSE)
  }
  chr_style <- startsWith(tolower(df$chr), "chr")
  if (any(chr_style)) {
    stop(sprintf("%s:%d: chr-style labels (e.g., chr1/chrX) are not allowed", path, which(chr_style)[[1]]), call. = FALSE)
  }
  known_chr <- df$chr %in% c(.canonical_chr_order, .ignored_chr)
  if (!all(known_chr)) {
    idx <- which(!known_chr)[[1]]
    stop(sprintf("%s:%d: unsupported chr label %s; expected 1-22, X (Y/MT are ignored)", path, idx, sQuote(df$chr[[idx]])), call. = FALSE)
  }
  .warn_allele_cap(df$a1, df$a2, path, seq_len(nrow(df)))
  bad_allele <- is.na(df$a1) | is.na(df$a2) | df$a1 == "" | df$a2 == ""
  if (any(bad_allele)) {
    stop(sprintf("%s:%d: a1 and a2 must be non-empty", path, which(bad_allele)[[1]]), call. = FALSE)
  }
  bad_a1 <- !grepl("^[ACGT]+$", df$a1)
  if (any(bad_a1)) {
    idx <- which(bad_a1)[[1]]
    stop(sprintf("%s:%d: a1 must be uppercase DNA bases (A/C/G/T): %s", path, idx, sQuote(df$a1[[idx]])), call. = FALSE)
  }
  bad_a2 <- !grepl("^[ACGT]+$", df$a2)
  if (any(bad_a2)) {
    idx <- which(bad_a2)[[1]]
    stop(sprintf("%s:%d: a2 must be uppercase DNA bases (A/C/G/T): %s", path, idx, sQuote(df$a2[[idx]])), call. = FALSE)
  }
  same_allele <- df$a1 == df$a2
  if (any(same_allele)) {
    stop(sprintf("%s:%d: a1 and a2 must differ", path, which(same_allele)[[1]]), call. = FALSE)
  }
  cm <- suppressWarnings(as.numeric(df$cm))
  if (any(is.na(cm) | !is.finite(cm))) {
    idx <- which(is.na(cm) | !is.finite(cm))[[1]]
    stop(sprintf("%s:%d: cm is not a number: %s", path, idx, sQuote(df$cm[[idx]])), call. = FALSE)
  }
  bp_num <- suppressWarnings(as.numeric(df$bp))
  bad_bp <- is.na(bp_num) | !is.finite(bp_num) | floor(bp_num) != bp_num
  if (any(bad_bp)) {
    idx <- which(bad_bp)[[1]]
    stop(sprintf("%s:%d: bp is not an integer: %s", path, idx, sQuote(df$bp[[idx]])), call. = FALSE)
  }
  df$bp <- as.integer(bp_num)
  df$a1_hash64 <- .allele_hash64(df$a1)
  df$a2_hash64 <- .allele_hash64(df$a2)
  .validate_bim_order_and_duplicates(df[df$chr %in% .canonical_chr_order, , drop = FALSE], path)
  df
}

.validate_bim_order_and_duplicates <- function(df, path) {
  if (!nrow(df)) {
    return(invisible(NULL))
  }
  rank <- unname(.chr_rank[df$chr])
  bad_order <- c(FALSE, rank[-1L] < rank[-length(rank)] | (rank[-1L] == rank[-length(rank)] & df$bp[-1L] < df$bp[-nrow(df)]))
  if (any(bad_order)) {
    stop(sprintf("%s:%d: rows must be sorted by (chr_rank, bp) in canonical contig order", path, which(bad_order)[[1]]), call. = FALSE)
  }

  for (label in unique(df$chr)) {
    keep <- df$chr == label
    .check_unique_variant_keys(
      df$bp[keep],
      df$a1_hash64[keep],
      df$a2_hash64[keep],
      sprintf("%s: shard %s", path, label),
      "variant"
    )
  }
  invisible(NULL)
}

.parse_fam <- function(path) {
  df <- tryCatch(
    utils::read.table(
      path,
      header = FALSE,
      sep = "",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      colClasses = rep("character", 6L),
      fill = FALSE
    ),
    error = function(e) {
      stop(sprintf("%s: expected 6 whitespace-delimited columns", path), call. = FALSE)
    }
  )
  if (ncol(df) != 6L) {
    stop(sprintf("%s: expected 6 whitespace-delimited columns, got %d", path, ncol(df)), call. = FALSE)
  }
  names(df) <- c("fid", "iid", "father_id", "mother_id", "sex", "pheno")
  bad_required <- df == ""
  if (any(bad_required)) {
    stop(sprintf("%s:%d: FAM fields must be non-empty; 6 whitespace-delimited columns are required", path, which(rowSums(bad_required) > 0L)[[1]]), call. = FALSE)
  }
  sex_num <- suppressWarnings(as.numeric(df$sex))
  bad_sex <- is.na(sex_num) | !is.finite(sex_num) | floor(sex_num) != sex_num | !(sex_num %in% c(0, 1, 2))
  if (any(bad_sex)) {
    stop(sprintf("%s:%d: FAM sex must be one of PLINK values 1, 2, or 0", path, which(bad_sex)[[1]]), call. = FALSE)
  }
  out <- df[c("fid", "iid", "father_id", "mother_id")]
  out$sex <- as.integer(sex_num)
  dup <- duplicated(paste(out$fid, out$iid, sep = "\r"))
  if (any(dup)) {
    stop(sprintf("%s:%d: duplicate FAM subject pair (fid, iid)", path, which(dup)[[1]]), call. = FALSE)
  }
  out
}

.parse_ploidy <- function(path, source_num_snp, source_chr = NULL) {
  if (is.null(path)) {
    male <- rep.int(2, source_num_snp)
    female <- rep.int(2, source_num_snp)
    if (!is.null(source_chr)) {
      if (length(source_chr) != source_num_snp) {
        stop("source_chr length must match source_num_snp", call. = FALSE)
      }
      male[as.character(source_chr) == "X"] <- 1
    }
    return(list(male = as.numeric(male), female = as.numeric(female)))
  }
  df <- tryCatch(
    utils::read.table(
      path,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      colClasses = rep("character", 2L),
      fill = FALSE
    ),
    error = function(e) {
      stop(sprintf("%s: expected 2 tab-delimited columns", path), call. = FALSE)
    }
  )
  if (ncol(df) != 2L) {
    stop(sprintf("%s: expected 2 tab-delimited columns, got %d", path, ncol(df)), call. = FALSE)
  }
  if (nrow(df) != source_num_snp) {
    stop(sprintf("%s: PLOIDY row count mismatch: expected %d, got %d", path, source_num_snp, nrow(df)), call. = FALSE)
  }
  values <- matrix(
    c(
      suppressWarnings(as.numeric(df[[1L]])),
      suppressWarnings(as.numeric(df[[2L]]))
    ),
    ncol = 2L
  )
  bad <- is.na(values) | !is.finite(values) | floor(values) != values | !(values %in% c(0, 1, 2))
  if (any(bad)) {
    stop(sprintf("%s:%d: ploidy values must be integers 0, 1, or 2", path, which(rowSums(bad) > 0L)[[1]]), call. = FALSE)
  }
  list(male = as.numeric(values[, 1]), female = as.numeric(values[, 2]))
}

.expected_bed_size <- function(source_num_sample, source_num_snp) {
  sample_count <- as.numeric(source_num_sample)
  snp_count <- as.numeric(source_num_snp)
  3 + floor((sample_count + 3) / 4) * snp_count
}

.validate_bed_file <- function(path, source_num_sample, source_num_snp) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 3L)
  if (length(magic) != 3L || !identical(magic[1:2], as.raw(c(0x6c, 0x1b)))) {
    stop(sprintf("%s: invalid PLINK BED magic bytes", path), call. = FALSE)
  }
  if (!identical(magic[[3]], as.raw(0x01))) {
    stop(sprintf("%s: PLINK BED must be SNP-major mode", path), call. = FALSE)
  }
  size <- file.info(path)$size
  expected <- .expected_bed_size(source_num_sample, source_num_snp)
  if (!identical(as.numeric(size), as.numeric(expected))) {
    stop(sprintf("%s: BED file size mismatch: expected %.0f, got %.0f", path, expected, size), call. = FALSE)
  }
  as.numeric(size)
}

.bed_lookup_int <- local({
  values <- 0:255
  out <- matrix(0L, nrow = 256L, ncol = 4L)
  for (j in 0:3) {
    two_bit <- bitwAnd(bitwShiftR(values, 2L * j), 3L)
    decoded <- integer(256L)
    decoded[two_bit == 0L] <- 2L
    decoded[two_bit == 1L] <- -1L
    decoded[two_bit == 2L] <- 1L
    decoded[two_bit == 3L] <- 0L
    out[, j + 1L] <- decoded
  }
  out
})

.bed_max_read_bytes <- 64 * 1024^2

.read_bed_rows_int8_codes <- function(path, source_rows0, source_num_sample) {
  source_rows0 <- as.integer(source_rows0)
  sample_count <- as.numeric(source_num_sample)
  bytes_per_snp <- floor((sample_count + 3) / 4)
  if (!length(source_rows0)) {
    return(matrix(0L, nrow = source_num_sample, ncol = 0L))
  }
  unique_rows0 <- sort(unique(source_rows0))
  inverse <- match(source_rows0, unique_rows0)
  unique_decoded <- matrix(0L, nrow = source_num_sample, ncol = length(unique_rows0))
  rows_per_read <- if (bytes_per_snp == 0) {
    length(unique_rows0)
  } else {
    max(1L, as.integer(floor(.bed_max_read_bytes / bytes_per_snp)))
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)

  breaks <- c(1L, which(diff(unique_rows0) != 1L) + 1L, length(unique_rows0) + 1L)
  for (g in seq_len(length(breaks) - 1L)) {
    group_cols <- seq.int(breaks[[g]], breaks[[g + 1L]] - 1L)
    chunk_starts <- seq.int(1L, length(group_cols), by = rows_per_read)
    for (chunk_start in chunk_starts) {
      chunk_stop <- min(chunk_start + rows_per_read - 1L, length(group_cols))
      block_cols <- group_cols[chunk_start:chunk_stop]
      block_rows0 <- unique_rows0[block_cols]
      byte_offset <- 3 + as.numeric(block_rows0[[1]]) * bytes_per_snp
      expected_bytes <- bytes_per_snp * length(block_rows0)
      if (expected_bytes > .Machine$integer.max) {
        stop(sprintf("%s: BED read block is too large: %.0f bytes", path, expected_bytes), call. = FALSE)
      }
      seek(con, where = byte_offset, origin = "start", rw = "read")
      packed <- readBin(con, what = "integer", n = as.integer(expected_bytes), size = 1L, signed = FALSE, endian = "little")
      if (length(packed) != expected_bytes) {
        stop(
          sprintf("%s: short BED read for source SNP rows %d-%.0f: expected %.0f bytes, got %d", path, block_rows0[[1]], block_rows0[[length(block_rows0)]], expected_bytes, length(packed)),
          call. = FALSE
        )
      }
      decoded_flat <- as.vector(t(.bed_lookup_int[packed + 1L, , drop = FALSE]))
      decoded_block <- matrix(decoded_flat, nrow = bytes_per_snp * 4, ncol = length(block_rows0))
      unique_decoded[, block_cols] <- decoded_block[seq_len(sample_count), , drop = FALSE]
    }
  }
  unique_decoded[, inverse, drop = FALSE]
}
