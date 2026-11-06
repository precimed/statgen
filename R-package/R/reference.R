.reference_cache_schema <- "reference_cache/0.1"

num_snp <- function(x, ...) UseMethod("num_snp")
shards <- function(x, ...) UseMethod("shards")
chr <- function(x, ...) UseMethod("chr")
snp <- function(x, ...) UseMethod("snp")
bp <- function(x, ...) UseMethod("bp")
a1 <- function(x, ...) UseMethod("a1")
a2 <- function(x, ...) UseMethod("a2")
a1_hash64 <- function(x, ...) UseMethod("a1_hash64")
a2_hash64 <- function(x, ...) UseMethod("a2_hash64")
is_single_nucleotide_variant <- function(x, ...) UseMethod("is_single_nucleotide_variant")
is_strand_ambiguous <- function(x, ...) UseMethod("is_strand_ambiguous")
shard_offsets <- function(x, ...) UseMethod("shard_offsets")
select_shards <- function(x, shards, ...) UseMethod("select_shards")
is_object_compatible <- function(reference, object, ...) UseMethod("is_object_compatible")
validate_checksums <- function(x, ...) UseMethod("validate_checksums")
save_cache <- function(x, path, ...) UseMethod("save_cache")

load_reference <- function(path, shards = NULL) {
  path <- .validate_path_scalar(path, "path")
  if (grepl("@", path, fixed = TRUE)) {
    available_labels <- character()
    available_paths <- character()
    for (label in .canonical_chr_order) {
      candidate <- gsub("@", label, path, fixed = TRUE)
      if (file.exists(candidate)) {
        available_labels <- c(available_labels, label)
        available_paths <- c(available_paths, candidate)
      }
    }
    if (!length(available_labels)) {
      stop(sprintf("No BIM shards found matching template: %s", path), call. = FALSE)
    }
    selected <- .validate_requested_shards(shards, available_labels, "load_reference")
    out <- vector("list", length(selected))
    for (i in seq_along(selected)) {
      label <- selected[[i]]
      bim <- .parse_bim(available_paths[[match(label, available_labels)]])
      out[[i]] <- .new_reference_shard(
        label, bim$chr, bim$snp, bim$bp, bim$a1, bim$a2,
        bim$a1_hash64, bim$a2_hash64
      )
    }
    for (shard in out) .validate_reference_keys(shard)
    return(.new_reference_panel(out))
  }

  bim <- .parse_bim(path)
  available_labels <- .canonical_chr_order[vapply(
    .canonical_chr_order,
    function(label) any(bim$chr == label),
    logical(1)
  )]
  selected <- .validate_requested_shards(shards, available_labels, "load_reference")
  out <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    label <- selected[[i]]
    keep <- bim$chr == label
    out[[i]] <- .new_reference_shard(
      label, bim$chr[keep], bim$snp[keep], bim$bp[keep], bim$a1[keep], bim$a2[keep],
      bim$a1_hash64[keep], bim$a2_hash64[keep]
    )
  }
  for (shard in out) .validate_reference_keys(shard)
  .new_reference_panel(out)
}

save_reference_cache <- function(panel, path) {
  if (!inherits(panel, "ReferencePanel")) {
    stop("panel must be a ReferencePanel", call. = FALSE)
  }
  path <- .validate_path_scalar(path, "path")
  offsets <- shard_offsets(panel)
  metadata <- list(
    schema = .reference_cache_schema,
    n_shards = length(panel$shards),
    shard_labels = vapply(panel$shards, function(s) s$label, character(1)),
    shard_checksums = vapply(panel$shards, function(s) s$checksum, character(1)),
    shard_start0 = offsets$start0,
    shard_stop0 = offsets$stop0
  )
  payload <- list(
    metadata = metadata,
    bp = bp(panel),
    snp_text_by_shard = vapply(panel$shards, function(s) .encode_text_payload(snp(s)), character(1)),
    a1_text_by_shard = vapply(panel$shards, function(s) .encode_text_payload(a1(s)), character(1)),
    a2_text_by_shard = vapply(panel$shards, function(s) .encode_text_payload(a2(s)), character(1)),
    a1_hash64 = a1_hash64(panel),
    a2_hash64 = a2_hash64(panel)
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, path)
  invisible(NULL)
}

load_reference_cache <- function(path, shards = NULL) {
  path <- .validate_path_scalar(path, "path")
  payload <- readRDS(path)
  .validate_reference_cache_payload(payload)
  meta <- payload$metadata
  selected <- .validate_requested_shards(shards, meta$shard_labels, "load_reference_cache")
  label_to_index <- stats::setNames(seq_along(meta$shard_labels), meta$shard_labels)

  out <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    label <- selected[[i]]
    idx <- label_to_index[[label]]
    start <- meta$shard_start0[[idx]] + 1L
    stop <- meta$shard_stop0[[idx]]
    row_count <- stop - meta$shard_start0[[idx]]
    out[[i]] <- .new_reference_shard_from_cache(
      label = label,
      num_snp = row_count,
      bp = payload$bp[start:stop],
      snp_text = payload$snp_text_by_shard[[idx]],
      a1_text = payload$a1_text_by_shard[[idx]],
      a2_text = payload$a2_text_by_shard[[idx]],
      checksum = meta$shard_checksums[[idx]],
      a1_hash64 = payload$a1_hash64[start:stop],
      a2_hash64 = payload$a2_hash64[start:stop]
    )
  }
  .new_reference_panel(out)
}

num_snp.ReferenceShard <- function(x, ...) x$num_snp
num_snp.ReferencePanel <- function(x, ...) x$num_snp
shards.ReferencePanel <- function(x, ...) x$shards
shards.default <- function(x, ...) {
  stop("object has no shards accessor", call. = FALSE)
}

chr.ReferenceShard <- function(x, ...) rep.int(x$label, x$num_snp)
chr.ReferencePanel <- function(x, ...) unlist(lapply(x$shards, chr), use.names = FALSE)

snp.ReferenceShard <- function(x, ...) .shard_string_vector(x, "snp")
snp.ReferencePanel <- function(x, ...) unlist(lapply(x$shards, snp), use.names = FALSE)

bp.ReferenceShard <- function(x, ...) x$bp
bp.ReferencePanel <- function(x, ...) unlist(lapply(x$shards, bp), use.names = FALSE)

a1.ReferenceShard <- function(x, ...) .shard_string_vector(x, "a1")
a1.ReferencePanel <- function(x, ...) unlist(lapply(x$shards, a1), use.names = FALSE)

a2.ReferenceShard <- function(x, ...) .shard_string_vector(x, "a2")
a2.ReferencePanel <- function(x, ...) unlist(lapply(x$shards, a2), use.names = FALSE)

a1_hash64.ReferenceShard <- function(x, ...) x$a1_hash64
a1_hash64.ReferencePanel <- function(x, ...) do.call(c, lapply(x$shards, a1_hash64))

a2_hash64.ReferenceShard <- function(x, ...) x$a2_hash64
a2_hash64.ReferencePanel <- function(x, ...) do.call(c, lapply(x$shards, a2_hash64))

is_single_nucleotide_variant.ReferenceShard <- function(x, ...) {
  .is_single_nucleotide_variant_hash(a1_hash64(x), a2_hash64(x))
}

is_single_nucleotide_variant.ReferencePanel <- function(x, ...) {
  unlist(lapply(x$shards, is_single_nucleotide_variant), use.names = FALSE)
}

is_strand_ambiguous.ReferenceShard <- function(x, ...) {
  .is_strand_ambiguous_hash(a1_hash64(x), a2_hash64(x))
}

is_strand_ambiguous.ReferencePanel <- function(x, ...) {
  unlist(lapply(x$shards, is_strand_ambiguous), use.names = FALSE)
}

shard_offsets.ReferencePanel <- function(x, ...) x$shard_offsets

select_shards.ReferencePanel <- function(x, shards, ...) {
  available <- vapply(x$shards, function(s) s$label, character(1))
  selected <- .validate_requested_shards(shards, available, "ReferencePanel.select_shards")
  by_label <- stats::setNames(x$shards, available)
  .new_reference_panel(unname(by_label[selected]))
}

validate_checksums.ReferencePanel <- function(x, ...) {
  for (shard in x$shards) {
    computed <- .checksum_from_arrays(chr(shard), bp(shard), a1(shard), a2(shard))
    if (!identical(computed, shard$checksum)) {
      stop(sprintf("Reference checksum mismatch for shard %s", shard$label), call. = FALSE)
    }
  }
  TRUE
}

is_object_compatible.ReferencePanel <- function(reference, object, ...) {
  obj_shards <- tryCatch(shards(object), error = function(e) NULL)
  if (is.null(obj_shards)) {
    warning("statgen: is_object_compatible: object has no shards accessor", call. = FALSE)
    return(FALSE)
  }
  if (length(obj_shards) != length(reference$shards)) {
    warning(
      sprintf(
        "statgen: is_object_compatible: shard count mismatch: reference has %d, object has %d",
        length(reference$shards), length(obj_shards)
      ),
      call. = FALSE
    )
    return(FALSE)
  }

  ok <- TRUE
  for (i in seq_along(reference$shards)) {
    ref_s <- reference$shards[[i]]
    obj_s <- obj_shards[[i]]
    obj_label <- obj_s$label
    if (is.null(obj_label) || !identical(obj_label, ref_s$label)) {
      warning(
        sprintf(
          "statgen: is_object_compatible: shard label mismatch: reference %s, object %s",
          ref_s$label,
          if (is.null(obj_label)) "<missing>" else as.character(obj_label)
        ),
        call. = FALSE
      )
      ok <- FALSE
      next
    }
    obj_num <- tryCatch(num_snp(obj_s), error = function(e) obj_s$num_snp)
    if (is.null(obj_num) || !identical(as.integer(obj_num), as.integer(ref_s$num_snp))) {
      warning(
        sprintf("statgen: is_object_compatible: shard %s: row count mismatch", ref_s$label),
        call. = FALSE
      )
      ok <- FALSE
      next
    }
    obj_checksum <- obj_s$reference_checksum
    if (is.null(obj_checksum)) {
      obj_checksum <- obj_s$checksum
    }
    if (!is.null(obj_checksum) && !identical(obj_checksum, ref_s$checksum)) {
      warning(
        sprintf("statgen: is_object_compatible: shard %s: reference_checksum mismatch", ref_s$label),
        call. = FALSE
      )
      ok <- FALSE
    }
  }
  ok
}

is_object_compatible.default <- function(reference, object, ...) {
  warning("statgen: is_object_compatible: reference must be a ReferencePanel", call. = FALSE)
  FALSE
}

save_cache.ReferencePanel <- function(x, path, ...) save_reference_cache(x, path)

print.ReferenceShard <- function(x, ...) {
  cat(sprintf("<ReferenceShard %s: %d SNPs>\n", x$label, x$num_snp))
  invisible(x)
}

print.ReferencePanel <- function(x, ...) {
  labels <- vapply(x$shards, function(s) s$label, character(1))
  cat(sprintf("<ReferencePanel: %d SNPs across %d shard(s): %s>\n", x$num_snp, length(labels), paste(labels, collapse = ", ")))
  invisible(x)
}

.new_reference_shard <- function(label, chr, snp, bp, a1, a2, a1_hash64 = NULL, a2_hash64 = NULL) {
  chr <- as.character(chr)
  snp <- as.character(snp)
  bp <- as.integer(bp)
  a1 <- as.character(a1)
  a2 <- as.character(a2)
  .validate_equal_lengths(list(chr = chr, snp = snp, bp = bp, a1 = a1, a2 = a2), "ReferenceShard vector lengths must match")
  if (length(chr) && any(chr != label)) {
    stop(sprintf("Reference shard %s contains multiple chr labels", label), call. = FALSE)
  }
  if (is.null(a1_hash64)) {
    a1_hash64 <- .allele_hash64(a1)
  } else {
    a1_hash64 <- bit64::as.integer64(a1_hash64)
  }
  if (is.null(a2_hash64)) {
    a2_hash64 <- .allele_hash64(a2)
  } else {
    a2_hash64 <- bit64::as.integer64(a2_hash64)
  }
  .validate_equal_lengths(
    list(bp = bp, a1_hash64 = a1_hash64, a2_hash64 = a2_hash64),
    "ReferenceShard allele hash vector lengths must match"
  )
  structure(
    list(
      label = label,
      num_snp = length(chr),
      snp = snp,
      bp = bp,
      a1 = a1,
      a2 = a2,
      a1_hash64 = a1_hash64,
      a2_hash64 = a2_hash64,
      checksum = .checksum_from_arrays(chr, bp, a1, a2)
    ),
    class = "ReferenceShard"
  )
}

.new_reference_shard_from_cache <- function(label, num_snp, bp, snp_text, a1_text, a2_text, checksum, a1_hash64, a2_hash64) {
  bp <- as.integer(bp)
  a1_hash64 <- bit64::as.integer64(a1_hash64)
  a2_hash64 <- bit64::as.integer64(a2_hash64)
  if (length(bp) != num_snp || length(a1_hash64) != num_snp || length(a2_hash64) != num_snp) {
    stop("Invalid reference cache: shard vector lengths mismatch", call. = FALSE)
  }
  structure(
    list(
      label = label,
      num_snp = as.integer(num_snp),
      snp_text = snp_text,
      bp = bp,
      a1_text = a1_text,
      a2_text = a2_text,
      a1_hash64 = a1_hash64,
      a2_hash64 = a2_hash64,
      checksum = checksum
    ),
    class = "ReferenceShard"
  )
}

.new_reference_panel <- function(shard_list) {
  if (!length(shard_list)) {
    stop("ReferencePanel requires at least one shard", call. = FALSE)
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
  structure(
    list(
      shards = shard_list,
      num_snp = pos,
      shard_offsets = offsets
    ),
    class = "ReferencePanel"
  )
}

.validate_reference_cache_payload <- function(payload) {
  cache <- .validate_cache_payload_metadata(payload, .reference_cache_schema, "reference")
  meta <- cache$meta
  n_shards <- cache$n_shards
  total <- cache$total
  payload_fields <- c("bp", "snp_text_by_shard", "a1_text_by_shard", "a2_text_by_shard", "a1_hash64", "a2_hash64")
  for (field in payload_fields) {
    if (is.null(payload[[field]])) {
      stop(sprintf("Invalid reference cache: missing %s", field), call. = FALSE)
    }
  }
  if (length(payload$snp_text_by_shard) != n_shards || length(payload$a1_text_by_shard) != n_shards || length(payload$a2_text_by_shard) != n_shards) {
    stop("Invalid reference cache: string payload shard count mismatch", call. = FALSE)
  }
  if (length(payload$bp) != total || length(payload$a1_hash64) != total || length(payload$a2_hash64) != total) {
    stop("Invalid reference cache: panel-wide vector lengths mismatch", call. = FALSE)
  }
}

.checksum_from_arrays <- function(chr, bp, a1, a2) {
  text <- paste0(chr, ":", as.integer(bp), ":", a1, ":", a2, "\n", collapse = "")
  digest::digest(enc2utf8(text), algo = "md5", serialize = FALSE)
}

.encode_text_payload <- function(values) {
  paste(as.character(values), collapse = "\n")
}

.decode_text_payload <- function(payload, expected_len, field) {
  if (!is.character(payload) || length(payload) != 1L || is.na(payload)) {
    stop(sprintf("Invalid reference cache: %s text payload must be a character scalar", field), call. = FALSE)
  }
  values <- if (identical(payload, "")) character() else strsplit(payload, "\n", fixed = TRUE)[[1]]
  if (length(values) != expected_len) {
    stop(sprintf("Invalid reference cache: decoded %s length mismatch", field), call. = FALSE)
  }
  values
}

.shard_string_vector <- function(x, field) {
  direct <- x[[field]]
  if (!is.null(direct)) {
    return(direct)
  }
  .decode_text_payload(x[[paste0(field, "_text")]], x$num_snp, field)
}

.validate_reference_keys <- function(shard) {
  # Ordered-key duplicates have already been rejected by the BIM parser.
  snv <- is_single_nucleotide_variant(shard)
  positions <- bp(shard)[snv]
  a <- a1_hash64(shard)[snv]
  b <- a2_hash64(shard)[snv]
  lo <- a
  hi <- b
  reverse <- a > b
  lo[reverse] <- b[reverse]
  hi[reverse] <- a[reverse]
  ord <- .variant_key_order(positions, lo, hi)
  left <- ord[-length(ord)]
  right <- ord[-1L]
  duplicate <- positions[left] == positions[right] & lo[left] == lo[right] & hi[left] == hi[right]
  if (any(duplicate)) {
    pos <- right[which(duplicate)[[1]]]
    stop(sprintf("reference shard %s: reciprocal SNV keys are not allowed at bp=%d",
                 shard$label, positions[pos]), call. = FALSE)
  }
  invisible(NULL)
}
