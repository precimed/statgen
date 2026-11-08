.sumstats_cache_schema <- "sumstats_cache/0.1"
.sumstats_numeric_na_values <- c("", "NA", "NaN")
.sumstats_required_cols <- c("chr", "bp", "a1", "a2", "p")
.sumstats_optional_cols <- c("z", "n", "beta", "se", "eaf", "info")
.sumstats_allowed_cols <- c(.sumstats_required_cols, .sumstats_optional_cols, "snp")
.sumstats_col_map <- c(pos = "bp", effectallele = "a1", otherallele = "a2")

logpvec <- function(x, ...) UseMethod("logpvec")
zvec <- function(x, ...) UseMethod("zvec")
nvec <- function(x, ...) UseMethod("nvec")
beta_vec <- function(x, ...) UseMethod("beta_vec")
se_vec <- function(x, ...) UseMethod("se_vec")
eaf_vec <- function(x, ...) UseMethod("eaf_vec")
info_vec <- function(x, ...) UseMethod("info_vec")

load_sumstats <- function(path, reference, allow_indel_swaps = FALSE, match_by = "genomic_position") {
  allow_indel_swaps <- .validate_matching_options(match_by, allow_indel_swaps)
  path <- .validate_path_scalar(path, "path")
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  df <- .parse_sumstats_tsv(path, match_by)
  .build_sumstats_panel(df, reference, path, allow_indel_swaps, match_by)
}

create_sumstats <- function(reference, p, z = NULL, n = NULL,
                            beta = NULL, se = NULL,
                            eaf = NULL, info = NULL) {
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  num_snp_value <- num_snp(reference)
  p <- .coerce_aligned_numeric_vector(p, num_snp_value, "p", allow_inf = FALSE)
  .validate_p_values(p, "p")
  aligned <- list(
    z = .coerce_optional_aligned_numeric_vector(z, num_snp_value, "z"),
    n = .coerce_optional_aligned_numeric_vector(n, num_snp_value, "n"),
    beta = .coerce_optional_aligned_numeric_vector(beta, num_snp_value, "beta"),
    se = .coerce_optional_aligned_numeric_vector(se, num_snp_value, "se"),
    eaf = .coerce_optional_aligned_numeric_vector(eaf, num_snp_value, "eaf"),
    info = .coerce_optional_aligned_numeric_vector(info, num_snp_value, "info")
  )
  for (field in .sumstats_optional_cols) {
    if (!is.null(aligned[[field]])) aligned[[field]][is.na(p)] <- NaN
  }
  .build_sumstats_from_aligned(reference, .derive_logp(p), aligned)
}

save_sumstats_cache <- function(sumstats, path) {
  if (!inherits(sumstats, "Sumstats")) {
    stop("sumstats must be a Sumstats object", call. = FALSE)
  }
  path <- .validate_path_scalar(path, "path")
  offsets <- shard_offsets(sumstats)
  metadata <- list(
    schema = .sumstats_cache_schema,
    n_shards = length(sumstats$shards),
    shard_labels = vapply(sumstats$shards, function(s) s$label, character(1)),
    shard_checksums = vapply(sumstats$shards, function(s) s$reference_checksum, character(1)),
    shard_start0 = offsets$start0,
    shard_stop0 = offsets$stop0,
    has_z = !is.null(zvec(sumstats)),
    has_n = !is.null(nvec(sumstats)),
    has_beta = !is.null(beta_vec(sumstats)),
    has_se = !is.null(se_vec(sumstats)),
    has_eaf = !is.null(eaf_vec(sumstats)),
    has_info = !is.null(info_vec(sumstats))
  )
  payload <- list(
    metadata = metadata,
    logpvec = logpvec(sumstats),
    zvec = zvec(sumstats),
    nvec = nvec(sumstats),
    beta_vec = beta_vec(sumstats),
    se_vec = se_vec(sumstats),
    eaf_vec = eaf_vec(sumstats),
    info_vec = info_vec(sumstats)
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, path)
  invisible(NULL)
}

load_sumstats_cache <- function(path, shards = NULL) {
  path <- .validate_path_scalar(path, "path")
  payload <- readRDS(path)
  .validate_sumstats_cache_payload(payload)
  meta <- payload$metadata
  selected <- .validate_requested_shards(shards, meta$shard_labels, "load_sumstats_cache")
  label_to_index <- stats::setNames(seq_along(meta$shard_labels), meta$shard_labels)

  out <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    label <- selected[[i]]
    idx <- label_to_index[[label]]
    start <- meta$shard_start0[[idx]] + 1L
    stop <- meta$shard_stop0[[idx]]
    out[[i]] <- .new_sumstats_shard(
      label = label,
      reference_checksum = meta$shard_checksums[[idx]],
      logpvec = payload$logpvec[start:stop],
      zvec = if (isTRUE(meta$has_z)) payload$zvec[start:stop] else NULL,
      nvec = if (isTRUE(meta$has_n)) payload$nvec[start:stop] else NULL,
      beta_vec = if (isTRUE(meta$has_beta)) payload$beta_vec[start:stop] else NULL,
      se_vec = if (isTRUE(meta$has_se)) payload$se_vec[start:stop] else NULL,
      eaf_vec = if (isTRUE(meta$has_eaf)) payload$eaf_vec[start:stop] else NULL,
      info_vec = if (isTRUE(meta$has_info)) payload$info_vec[start:stop] else NULL
    )
  }
  panel <- .new_sumstats_panel(out)
  .warn_optional_zn_completeness(zvec(panel), nvec(panel), logpvec(panel), "load_sumstats_cache")
  panel
}

num_snp.SumstatsShard <- function(x, ...) x$num_snp
num_snp.Sumstats <- function(x, ...) x$num_snp
shards.Sumstats <- function(x, ...) x$shards
shard_offsets.Sumstats <- function(x, ...) x$shard_offsets

logpvec.SumstatsShard <- function(x, ...) x$logpvec
logpvec.Sumstats <- function(x, ...) unlist(lapply(x$shards, logpvec), use.names = FALSE)

zvec.SumstatsShard <- function(x, ...) x$zvec
zvec.Sumstats <- function(x, ...) .optional_sumstats_concat(x, "zvec")

nvec.SumstatsShard <- function(x, ...) x$nvec
nvec.Sumstats <- function(x, ...) .optional_sumstats_concat(x, "nvec")

is_present.SumstatsShard <- function(x, ...) !is.nan(x$logpvec)
is_present.Sumstats <- function(x, ...) !is.nan(logpvec(x))

beta_vec.SumstatsShard <- function(x, ...) x$beta_vec
beta_vec.Sumstats <- function(x, ...) .optional_sumstats_concat(x, "beta_vec")

se_vec.SumstatsShard <- function(x, ...) x$se_vec
se_vec.Sumstats <- function(x, ...) .optional_sumstats_concat(x, "se_vec")

eaf_vec.SumstatsShard <- function(x, ...) x$eaf_vec
eaf_vec.Sumstats <- function(x, ...) .optional_sumstats_concat(x, "eaf_vec")

info_vec.SumstatsShard <- function(x, ...) x$info_vec
info_vec.Sumstats <- function(x, ...) .optional_sumstats_concat(x, "info_vec")

select_shards.Sumstats <- function(x, shards, ...) {
  available <- vapply(x$shards, function(s) s$label, character(1))
  selected <- .validate_requested_shards(shards, available, "Sumstats.select_shards")
  by_label <- stats::setNames(x$shards, available)
  .new_sumstats_panel(unname(by_label[selected]))
}

save_cache.Sumstats <- function(x, path, ...) save_sumstats_cache(x, path)

print.Sumstats <- function(x, ...) {
  labels <- vapply(x$shards, function(s) s$label, character(1))
  cat(sprintf("<Sumstats: %d SNPs across %d shard(s): %s>\n", x$num_snp, length(labels), paste(labels, collapse = ", ")))
  invisible(x)
}

.optional_sumstats_concat <- function(x, field) {
  present <- vapply(x$shards, function(s) !is.null(s[[field]]), logical(1))
  if (!any(present)) {
    return(NULL)
  }
  if (!all(present)) {
    stop(sprintf("Invalid Sumstats object: inconsistent %s presence across shards", field), call. = FALSE)
  }
  unlist(lapply(x$shards, function(s) s[[field]]), use.names = FALSE)
}

.new_sumstats_shard <- function(label, reference_checksum, logpvec,
                                zvec = NULL, nvec = NULL, beta_vec = NULL,
                                se_vec = NULL, eaf_vec = NULL, info_vec = NULL) {
  logpvec <- as.numeric(logpvec)
  n <- length(logpvec)
  structure(
    list(
      label = as.character(label),
      reference_checksum = as.character(reference_checksum),
      num_snp = as.integer(n),
      logpvec = logpvec,
      zvec = .coerce_optional_shard_numeric_vector(zvec, n, "zvec"),
      nvec = .coerce_optional_shard_numeric_vector(nvec, n, "nvec"),
      beta_vec = .coerce_optional_shard_numeric_vector(beta_vec, n, "beta_vec"),
      se_vec = .coerce_optional_shard_numeric_vector(se_vec, n, "se_vec"),
      eaf_vec = .coerce_optional_shard_numeric_vector(eaf_vec, n, "eaf_vec"),
      info_vec = .coerce_optional_shard_numeric_vector(info_vec, n, "info_vec")
    ),
    class = "SumstatsShard"
  )
}

.new_sumstats_panel <- function(shard_list) {
  if (!length(shard_list)) {
    stop("Sumstats requires at least one shard", call. = FALSE)
  }
  for (field in c("zvec", "nvec", "beta_vec", "se_vec", "eaf_vec", "info_vec")) {
    present <- vapply(shard_list, function(s) !is.null(s[[field]]), logical(1))
    if (any(present) && !all(present)) {
      stop(sprintf("Sumstats shards must have consistent %s presence", field), call. = FALSE)
    }
  }
  pos <- 0L
  offsets <- data.frame(
    shard_label = character(length(shard_list)),
    start0 = integer(length(shard_list)),
    stop0 = integer(length(shard_list)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(shard_list)) {
    offsets$shard_label[[i]] <- shard_list[[i]]$label
    offsets$start0[[i]] <- pos
    pos <- pos + shard_list[[i]]$num_snp
    offsets$stop0[[i]] <- pos
  }
  structure(list(shards = shard_list, num_snp = pos, shard_offsets = offsets), class = "Sumstats")
}

.parse_sumstats_tsv <- function(path, match_by = "genomic_position") {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  raw_names <- .read_sumstats_header(path)
  canonical_names <- .canonicalize_sumstats_columns(raw_names, path)
  id_mode <- match_by == "id"
  required <- if (id_mode) c("snp", "a1", "a2", "p") else .sumstats_required_cols
  missing <- setdiff(required, canonical_names)
  if (length(missing)) {
    stop(sprintf("%s: missing required columns: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }
  selected <- canonical_names %in% c(required, .sumstats_optional_cols)
  selected_raw <- raw_names[selected]
  selected_names <- canonical_names[selected]
  char_cols <- selected_raw[selected_names %in% c("chr", "snp", "a1", "a2")]
  numeric_cols <- selected_raw[!(selected_names %in% c("chr", "snp", "a1", "a2"))]

  fread_args <- c(
    .sumstats_fread_input(path),
    list(
      select = selected_raw,
      colClasses = list(character = char_cols, numeric = numeric_cols),
      data.table = FALSE,
      showProgress = FALSE,
      na.strings = .sumstats_numeric_na_values
    )
  )
  stopped_early <- FALSE
  df <- tryCatch(
    withCallingHandlers(
      do.call(data.table::fread, fread_args),
      warning = function(w) {
        message <- conditionMessage(w)
        if (grepl("Stopped early", message, fixed = TRUE)) {
          stopped_early <<- TRUE
          invokeRestart("muffleWarning")
        }
        if (startsWith(message, "Attempt to override column")) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {
      stop(sprintf("%s: malformed TSV", path), call. = FALSE)
    }
  )
  if (stopped_early) {
    stop(sprintf("%s: malformed TSV", path), call. = FALSE)
  }
  names(df) <- selected_names

  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(sprintf("%s: missing required columns: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }

  if (!id_mode) {
    bad_chr <- is.na(df$chr) | df$chr == ""
    if (any(bad_chr)) {
      stop(sprintf("%s: row %d: chr must be non-empty", path, which(bad_chr)[[1]] + 1L), call. = FALSE)
    }
    .validate_source_chr_labels(df$chr, path)
    df <- df[df$chr %in% .canonical_chr_order, , drop = FALSE]
  }

  for (field in c("a1", "a2")) {
    bad <- is.na(df[[field]]) | df[[field]] == ""
    if (any(bad)) {
      row <- as.integer(rownames(df)[which(bad)[[1]]]) + 1L
      stop(sprintf("%s: row %d: %s must be non-empty", path, row, field), call. = FALSE)
    }
  }
  .warn_allele_cap(df$a1, df$a2, path, as.integer(rownames(df)) + 1L)
  bad_a1 <- !grepl("^[ACGT]+$", df$a1)
  if (any(bad_a1)) {
    idx <- which(bad_a1)[[1]]
    row <- as.integer(rownames(df)[idx]) + 1L
    stop(sprintf("%s: row %d: a1 must be uppercase DNA bases (A/C/G/T): %s", path, row, sQuote(df$a1[[idx]])), call. = FALSE)
  }
  bad_a2 <- !grepl("^[ACGT]+$", df$a2)
  if (any(bad_a2)) {
    idx <- which(bad_a2)[[1]]
    row <- as.integer(rownames(df)[idx]) + 1L
    stop(sprintf("%s: row %d: a2 must be uppercase DNA bases (A/C/G/T): %s", path, row, sQuote(df$a2[[idx]])), call. = FALSE)
  }
  same_allele <- df$a1 == df$a2
  if (any(same_allele)) {
    row <- as.integer(rownames(df)[which(same_allele)[[1]]]) + 1L
    stop(sprintf("%s: row %d: a1 and a2 must differ", path, row), call. = FALSE)
  }

  if (!id_mode) {
    bp <- .parse_sumstats_numeric(df$bp, path, "bp", rownames(df))
    bad_bp <- is.na(bp) | !is.finite(bp) | floor(bp) != bp
    if (any(bad_bp)) {
      idx <- which(bad_bp)[[1]]
      row <- as.integer(rownames(df)[idx]) + 1L
      stop(sprintf("%s: row %d: bp is not an integer: %s", path, row, sQuote(df$bp[[idx]])), call. = FALSE)
    }
    df$bp <- as.integer(bp)
  }

  p <- .parse_sumstats_numeric(df$p, path, "p", rownames(df))
  bad_p <- !is.na(p) & (!is.finite(p) | p < 0 | p > 1)
  if (any(bad_p)) {
    idx <- which(bad_p)[[1]]
    row <- as.integer(rownames(df)[idx]) + 1L
    stop(sprintf("%s: row %d: p must be finite numeric in [0, 1]: %s", path, row, sQuote(df$p[[idx]])), call. = FALSE)
  }
  df$p <- p

  for (field in .sumstats_optional_cols) {
    if (field %in% names(df)) {
      value <- .parse_sumstats_numeric(df[[field]], path, field, rownames(df))
      value[!is.finite(value)] <- NaN
      df[[field]] <- value
    }
  }
  missing_p <- is.na(df$p)
  if (any(missing_p)) {
    .info(sprintf("%s: dropped %d row(s) with missing p", path, sum(missing_p)))
  }
  df[!missing_p, , drop = FALSE]
}

.parse_sumstats_numeric <- function(values, path, field, source_rows) {
  out <- suppressWarnings(as.numeric(values))
  if (!is.numeric(values)) {
    missing <- is.na(values) | values %in% .sumstats_numeric_na_values
    bad <- is.na(out) & !missing
    if (any(bad)) {
      idx <- which(bad)[[1]]
      row <- as.integer(source_rows[idx]) + 1L
      reason <- if (field == "bp") "is not an integer" else "must be numeric or missing"
      stop(sprintf("%s: row %d: %s %s: %s", path, row, field, reason, sQuote(values[[idx]])), call. = FALSE)
    }
  }
  out[is.na(out)] <- NaN
  out
}

.read_sumstats_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- tryCatch(readLines(con, n = 1L, warn = FALSE), error = function(e) character())
  if (!length(header)) {
    stop(sprintf("%s: malformed TSV", path), call. = FALSE)
  }
  strsplit(header[[1]], "\t", fixed = TRUE)[[1]]
}

.sumstats_fread_input <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE) || requireNamespace("R.utils", quietly = TRUE)) {
    return(list(input = path))
  }
  gzip <- Sys.which("gzip")
  if (!nzchar(gzip)) {
    stop(sprintf("%s: gzip input requires R.utils or a gzip executable", path), call. = FALSE)
  }
  cmd <- paste(shQuote(gzip), "-dc", shQuote(path))
  list(cmd = cmd)
}

.canonicalize_sumstats_columns <- function(columns, path) {
  out <- tolower(as.character(columns))
  mapped <- .sumstats_col_map[out]
  out[!is.na(mapped)] <- unname(mapped[!is.na(mapped)])
  duplicates <- unique(out[duplicated(out)])
  if (length(duplicates)) {
    stop(sprintf("%s: duplicate columns after column normalization: %s", path, paste(sort(duplicates), collapse = ", ")), call. = FALSE)
  }
  out
}

.validate_source_chr_labels <- function(chr_values, path) {
  chr_style <- startsWith(tolower(chr_values), "chr")
  if (any(chr_style)) {
    stop(sprintf("%s: row %d: chr-style labels (e.g., chr1/chrX) are not allowed", path, which(chr_style)[[1]] + 1L), call. = FALSE)
  }
  known <- chr_values %in% c(.canonical_chr_order, .ignored_chr)
  if (!all(known)) {
    idx <- which(!known)[[1]]
    stop(sprintf("%s: row %d: unsupported chr label %s; expected 1-22, X (Y/MT are ignored)", path, idx + 1L, sQuote(chr_values[[idx]])), call. = FALSE)
  }
}

.build_sumstats_panel <- function(df, reference, source_path, allow_indel_swaps, match_by = "genomic_position") {
  n <- num_snp(reference)
  aligned_p <- rep.int(NaN, n)
  aligned_optional <- lapply(.sumstats_optional_cols, function(field) {
    if (field %in% names(df)) rep.int(NaN, n) else NULL
  })
  names(aligned_optional) <- .sumstats_optional_cols

  id_mode <- match_by == "id"
  src_a1_hash64 <- .allele_hash64(df$a1)
  src_a2_hash64 <- .allele_hash64(df$a2)
  if (id_mode) {
    id_match <- .match_id_panel(reference, list(list(
      context = source_path, ids = df$snp, a1 = src_a1_hash64, a2 = src_a2_hash64
    )))
  } else {
    src_chr <- as.character(df$chr)
    keep_supported <- src_chr %in% .canonical_chr_order
    source_rows_by_chr <- split(which(keep_supported), src_chr[keep_supported])
  }

  ref_shards <- shards(reference)
  offsets <- shard_offsets(reference)
  for (i in seq_along(ref_shards)) {
    ref_shard <- ref_shards[[i]]
    start <- offsets$start0[[i]] + 1L
    stop <- offsets$stop0[[i]]
    if (id_mode) {
      ix <- seq.int(start, stop)
      local_match <- id_match$source_rows[ix]
      swapped <- id_match$is_allele_swapped[ix]
      has_match <- !is.na(local_match)
      matched_src <- local_match[has_match]
    } else {
      src_idx <- source_rows_by_chr[[ref_shard$label]]
      if (is.null(src_idx)) src_idx <- integer()
      match <- .match_variant_keys(
        bp(ref_shard), a1_hash64(ref_shard), a2_hash64(ref_shard),
        df$bp[src_idx], src_a1_hash64[src_idx], src_a2_hash64[src_idx],
        ref_shard$label, "sumstats", allow_indel_swaps
      )
      local_match <- match$source_rows
      swapped <- match$is_allele_swapped
      .report_variant_swaps(source_path, ref_shard$label, sum(swapped),
                            match$num_strand_ambiguous, match$num_disabled_indels)
      has_match <- !is.na(local_match)
      matched_src <- src_idx[local_match[has_match]]
    }
    aligned_idx <- seq.int(start, stop)[has_match]
    aligned_p[aligned_idx] <- df$p[matched_src]
    for (field in .sumstats_optional_cols) {
      if (!is.null(aligned_optional[[field]])) {
        aligned_optional[[field]][aligned_idx] <- df[[field]][matched_src]
      }
    }

    swapped_idx <- seq.int(start, stop)[swapped]
    for (field in c("z", "beta", "eaf")) {
      if (!is.null(aligned_optional[[field]])) {
        values <- aligned_optional[[field]][swapped_idx]
        aligned_optional[[field]][swapped_idx] <- if (field == "eaf") 1 - values else -values
      }
    }
    if (!any(has_match)) {
      warning(sprintf("%s: no matching source rows for reference shard '%s'", source_path, ref_shard$label), call. = FALSE)
    }
  }

  logp <- .derive_logp(aligned_p)
  .warn_optional_zn_completeness(aligned_optional$z, aligned_optional$n, logp, "load_sumstats")
  .build_sumstats_from_aligned(reference, logp, aligned_optional)
}

.build_sumstats_from_aligned <- function(reference, aligned_logp, aligned_optional) {
  ref_shards <- shards(reference)
  offsets <- shard_offsets(reference)
  out <- vector("list", length(ref_shards))
  for (i in seq_along(ref_shards)) {
    ref_shard <- ref_shards[[i]]
    start <- offsets$start0[[i]] + 1L
    stop <- offsets$stop0[[i]]
    out[[i]] <- .new_sumstats_shard(
      label = ref_shard$label,
      reference_checksum = ref_shard$checksum,
      logpvec = aligned_logp[start:stop],
      zvec = if (is.null(aligned_optional$z)) NULL else aligned_optional$z[start:stop],
      nvec = if (is.null(aligned_optional$n)) NULL else aligned_optional$n[start:stop],
      beta_vec = if (is.null(aligned_optional$beta)) NULL else aligned_optional$beta[start:stop],
      se_vec = if (is.null(aligned_optional$se)) NULL else aligned_optional$se[start:stop],
      eaf_vec = if (is.null(aligned_optional$eaf)) NULL else aligned_optional$eaf[start:stop],
      info_vec = if (is.null(aligned_optional$info)) NULL else aligned_optional$info[start:stop]
    )
  }
  .new_sumstats_panel(out)
}

.derive_logp <- function(p) {
  out <- rep.int(NaN, length(p))
  zero <- !is.na(p) & p == 0
  positive <- !is.na(p) & p > 0
  out[zero] <- Inf
  out[positive] <- -log10(p[positive])
  out
}

.validate_p_values <- function(p, name) {
  bad <- !is.na(p) & (!is.finite(p) | p < 0 | p > 1)
  if (any(bad)) {
    stop(sprintf("%s[%d] must be finite numeric in [0, 1]", name, which(bad)[[1]]), call. = FALSE)
  }
}

.coerce_aligned_numeric_vector <- function(vec, n, name, allow_inf = FALSE) {
  arr <- suppressWarnings(as.numeric(vec))
  if (any(is.na(arr) & !is.na(vec))) {
    stop(sprintf("%s must be numeric", name), call. = FALSE)
  }
  if (length(arr) != n) {
    stop(sprintf("%s length mismatch: expected %d, got %d", name, n, length(arr)), call. = FALSE)
  }
  if (!allow_inf) {
    bad <- is.infinite(arr)
    if (any(bad)) {
      stop(sprintf("%s[%d] must be finite numeric or NaN", name, which(bad)[[1]]), call. = FALSE)
    }
  }
  arr
}

.coerce_optional_aligned_numeric_vector <- function(vec, n, name) {
  if (is.null(vec)) {
    return(NULL)
  }
  .coerce_aligned_numeric_vector(vec, n, name, allow_inf = FALSE)
}

.coerce_optional_shard_numeric_vector <- function(vec, n, name) {
  if (is.null(vec)) {
    return(NULL)
  }
  arr <- as.numeric(vec)
  if (length(arr) != n) {
    stop(sprintf("%s length mismatch: expected %d, got %d", name, n, length(arr)), call. = FALSE)
  }
  arr
}

.warn_optional_zn_completeness <- function(zvec, nvec, logpvec, context) {
  present <- !is.nan(logpvec)
  for (entry in list(list(name = "zvec", value = zvec), list(name = "nvec", value = nvec))) {
    if (is.null(entry$value)) {
      warning(sprintf("%s: %s is absent", context, entry$name), call. = FALSE)
    } else if (any(is.nan(entry$value[present]))) {
      warning(sprintf("%s: %s has missing values among present sumstats variants", context, entry$name), call. = FALSE)
    }
  }
}

.validate_sumstats_cache_payload <- function(payload) {
  cache <- .validate_cache_payload_metadata(payload, .sumstats_cache_schema, "sumstats")
  meta <- cache$meta
  total <- cache$total
  if (is.null(payload$logpvec) || length(payload$logpvec) != total) {
    stop("Invalid sumstats cache: logpvec length mismatch", call. = FALSE)
  }
  for (field in c("z", "n", "beta", "se", "eaf", "info")) {
    flag <- paste0("has_", field)
    vector_name <- if (field %in% c("z", "n")) paste0(field, "vec") else paste0(field, "_vec")
    if (isTRUE(meta[[flag]]) && (is.null(payload[[vector_name]]) || length(payload[[vector_name]]) != total)) {
      stop(sprintf("Invalid sumstats cache: %s length mismatch", vector_name), call. = FALSE)
    }
  }
}
