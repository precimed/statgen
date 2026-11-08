# Sort exact 32-bit halves so multi-key order() never interprets integer64 hashes as doubles.
.variant_key_order <- function(bp, a1_hash64, a2_hash64) {
  two32 <- bit64::as.integer64(4294967296)
  order(
    bp,
    as.double(a1_hash64 %/% two32), as.double(a1_hash64 %% two32),
    as.double(a2_hash64 %/% two32), as.double(a2_hash64 %% two32)
  )
}

.check_unique_variant_keys <- function(bp, a1_hash64, a2_hash64, where, object_name = "variant") {
  .validate_equal_lengths(
    list(bp = bp, a1_hash64 = a1_hash64, a2_hash64 = a2_hash64),
    "numeric variant key vector lengths mismatch"
  )
  if (length(bp) < 2L) {
    return(invisible(NULL))
  }
  ord <- .variant_key_order(bp, a1_hash64, a2_hash64)
  duplicate <- bp[ord][-1L] == bp[ord][-length(ord)] &
    a1_hash64[ord][-1L] == a1_hash64[ord][-length(ord)] &
    a2_hash64[ord][-1L] == a2_hash64[ord][-length(ord)]
  if (any(duplicate)) {
    stop(sprintf("%s: duplicate %s matching key", where, object_name), call. = FALSE)
  }
  invisible(NULL)
}

# Reference payloads are already validated by load_reference; caches are trusted.
.match_variant_keys <- function(ref_bp, ref_a1_hash64, ref_a2_hash64,
                                src_bp, src_a1_hash64, src_a2_hash64,
                                shard_label, source_name, allow_indel_swaps = FALSE) {
  rows <- .lookup_variant_candidates(
    ref_bp, ref_a1_hash64, ref_a2_hash64,
    src_bp, src_a1_hash64, src_a2_hash64, seq_along(src_bp),
    rep.int(TRUE, length(src_bp)), shard_label, source_name
  )$rows
  swapped <- rep.int(FALSE, length(ref_bp))
  unused <- rep.int(TRUE, length(src_bp))
  unused[rows[!is.na(rows)]] <- FALSE
  snv <- .is_single_nucleotide_variant_hash(src_a1_hash64, src_a2_hash64)
  indel <- xor(.is_single_base_hash(src_a1_hash64), .is_single_base_hash(src_a2_hash64))
  eligible <- snv | (indel & allow_indel_swaps)
  source_rows <- which(unused & (snv | indel))
  missing <- which(is.na(rows))
  fallback <- .lookup_variant_candidates(
    ref_bp[missing], ref_a1_hash64[missing], ref_a2_hash64[missing],
    src_bp[source_rows], src_a2_hash64[source_rows], src_a1_hash64[source_rows],
    source_rows, eligible[source_rows], shard_label, source_name
  )
  accepted <- missing[!is.na(fallback$rows)]
  rows[accepted] <- fallback$rows[!is.na(fallback$rows)]
  swapped[accepted] <- TRUE
  selected <- rows[swapped]
  list(
    source_rows = rows,
    is_allele_swapped = swapped,
    num_strand_ambiguous = sum(.is_strand_ambiguous_hash(src_a1_hash64[selected], src_a2_hash64[selected])),
    num_disabled_indels = fallback$excluded_hits
  )
}

.lookup_variant_candidates <- function(ref_bp, ref_a1, ref_a2, src_bp, src_a1, src_a2,
                                       source_rows, eligible, label, context) {
  n_ref <- length(ref_bp)
  rows <- rep.int(NA_integer_, n_ref)
  if (!n_ref || !length(src_bp)) {
    return(list(rows = rows, excluded_hits = 0L))
  }
  bp <- c(ref_bp, src_bp)
  a1 <- c(bit64::as.integer64(ref_a1), bit64::as.integer64(src_a1))
  a2 <- c(bit64::as.integer64(ref_a2), bit64::as.integer64(src_a2))
  is_ref <- c(rep.int(TRUE, n_ref), rep.int(FALSE, length(src_bp)))
  eligible_all <- c(rep.int(FALSE, n_ref), eligible)
  local_index <- c(seq_len(n_ref), source_rows)
  ord <- .variant_key_order(bp, a1, a2)
  left <- ord[-length(ord)]
  right <- ord[-1L]
  group <- cumsum(c(TRUE, !(bp[left] == bp[right] & a1[left] == a1[right] & a2[left] == a2[right])))
  n_groups <- group[length(group)]
  ref_count <- tabulate(group[is_ref[ord]], n_groups)
  src_count <- tabulate(group[!is_ref[ord]], n_groups)
  eligible_count <- tabulate(group[eligible_all[ord]], n_groups)
  ambiguous <- ref_count > 0L & eligible_count > 1L
  if (any(ambiguous)) {
    pos <- ord[match(which(ambiguous)[[1]], group)]
    stop(sprintf("Ambiguous duplicate %s/reference matching key in shard %s: bp=%d, a1_hash64=%s, a2_hash64=%s",
                 context, label, bp[pos], as.character(a1[pos]), as.character(a2[pos])), call. = FALSE)
  }
  src_by_group <- rep.int(NA_integer_, n_groups)
  src <- eligible_all[ord]
  src_by_group[group[src]] <- local_index[ord[src]]
  ref <- is_ref[ord]
  matched <- eligible_count[group[ref]] == 1L
  rows[local_index[ord[ref]][matched]] <- src_by_group[group[ref][matched]]
  list(rows = rows, excluded_hits = sum(src_count[ref_count > 0L & eligible_count == 0L]))
}

.id_matching_eligible <- function(ids, a1, a2, context) {
  snv <- .is_single_nucleotide_variant_hash(a1, a2)
  palindromic <- .is_strand_ambiguous_hash(a1, a2)
  if (any(palindromic)) .info(sprintf("%s: excluded %d palindromic SNV(s)", context, sum(palindromic)))
  if (any(!snv)) .info(sprintf("%s: excluded %d multi-base variant(s)", context, sum(!snv)))
  alleles_ok <- snv & !palindromic & a1 != a2
  missing <- is.na(ids) | ids %in% c("", ".", "NA", "NaN")
  valid <- !missing & grepl("\\A[!-~]{1,32}\\z", ids, perl = TRUE, useBytes = TRUE)
  n_missing <- sum(alleles_ok & missing)
  n_invalid <- sum(alleles_ok & !valid & !missing)
  if (n_missing) warning(sprintf("%s: excluded %d unusable missing, empty, or placeholder ID(s)", context, n_missing), call. = FALSE)
  if (n_invalid) warning(sprintf("%s: excluded %d unusable overlength or invalid-character ID(s)", context, n_invalid), call. = FALSE)
  alleles_ok & valid
}

.match_id_panel <- function(reference, sources, chromosome_scope = FALSE) {
  ref_a1 <- a1_hash64(reference)
  ref_a2 <- a2_hash64(reference)
  ref_ids <- snp(reference)  # Decode packed reference IDs once, only in ID mode.
  ri <- which(.id_matching_eligible(ref_ids, ref_a1, ref_a2, "reference panel"))
  labels <- vapply(shards(reference), function(s) s$label, character(1))
  eligible <- lapply(sources, function(source) {
    scope <- if (chromosome_scope) source$chr %in% labels else rep.int(TRUE, length(source$ids))
    ok <- rep.int(FALSE, length(source$ids))
    ok[scope] <- .id_matching_eligible(source$ids[scope], source$a1[scope], source$a2[scope], source$context)
    ok
  })
  src_ids <- unlist(lapply(sources, function(s) s$ids), use.names = FALSE)
  src_a1 <- do.call(c, lapply(sources, function(s) s$a1))
  src_a2 <- do.call(c, lapply(sources, function(s) s$a2))
  si <- which(unlist(eligible, use.names = FALSE))
  rows <- rep.int(NA_integer_, num_snp(reference))
  swapped <- rep.int(FALSE, length(rows))
  if (!length(ri)) return(list(source_rows = rows, is_allele_swapped = swapped))

  ids <- c(ref_ids[ri], src_ids[si])
  unique_ids <- unique(ids)
  codes <- match(ids, unique_ids)
  n_codes <- as.double(length(unique_ids))
  # Packed integer codes stay exact in double; allele hashes remain integer64.
  if (n_codes > (2^53 - 1) / (2 * (length(labels) + 1))) {
    stop("Too many ID codes for exact numeric matching keys", call. = FALSE)
  }
  r1 <- bit64::match(ref_a1[ri], .single_base_hash64)
  r2 <- bit64::match(ref_a2[ri], .single_base_hash64)
  s1 <- bit64::match(src_a1[si], .single_base_hash64)
  s2 <- bit64::match(src_a2[si], .single_base_hash64)
  ref_keys <- as.double(codes[seq_along(ri)]) * 2 + (r1 + r2) %% 2L
  src_keys <- as.double(codes[length(ri) + seq_along(si)]) * 2 + (s1 + s2) %% 2L
  duplicate <- duplicated(ref_keys) | duplicated(ref_keys, fromLast = TRUE)
  if (any(duplicate)) {
    warning(sprintf("reference panel: excluded %d slots in %d duplicate ID key groups",
                    sum(duplicate), length(unique(ref_keys[duplicate]))), call. = FALSE)
  }
  ri <- ri[!duplicate]
  ref_keys <- ref_keys[!duplicate]
  if (!length(ri) || !length(si)) return(list(source_rows = rows, is_allele_swapped = swapped))
  if (chromosome_scope) {
    ref_chr <- rep.int(seq_along(labels), vapply(shards(reference), num_snp, integer(1)))
    src_chr <- match(unlist(lapply(sources, function(s) s$chr), use.names = FALSE), labels)
    ref_keys <- ref_chr[ri] * (n_codes * 2) + ref_keys
    src_keys <- src_chr[si] * (n_codes * 2) + src_keys
  }
  loc <- match(src_keys, ref_keys)
  hit <- !is.na(loc)
  counts <- tabulate(loc[hit], length(ri))
  if (any(counts > 1L)) {
    slot <- which(counts > 1L)[[1]]
    contexts <- vapply(sources, function(s) s$context, character(1))
    stop(sprintf("Ambiguous duplicate source/reference ID matching key for ID %s in %s",
                 sQuote(ref_ids[ri[slot]]), paste(contexts, collapse = ", ")), call. = FALSE)
  }
  rows[ri[loc[hit]]] <- si[hit]
  slots <- which(!is.na(rows))
  source_first <- bit64::match(src_a1[rows[slots]], .single_base_hash64)
  target_first <- bit64::match(ref_a1[slots], .single_base_hash64)
  target_other <- bit64::match(ref_a2[slots], .single_base_hash64)
  swapped[slots] <- source_first == target_other | 5L - source_first == target_other
  complemented <- rep.int(FALSE, length(rows))
  complemented[slots] <- source_first != target_first & source_first != target_other
  offsets <- c(0, cumsum(vapply(sources, function(s) length(s$ids), integer(1))))
  for (i in seq_along(sources)) {
    selected <- !is.na(rows) & rows > offsets[i] & rows <= offsets[i + 1L]
    for (label in c("inversion", "strand complementation")) {
      mask <- if (label == "inversion") swapped else complemented
      count <- sum(selected & mask)
      if (count) .info(sprintf("%s: %d matched row(s) requiring %s", sources[[i]]$context, count, label))
    }
  }
  list(source_rows = rows, is_allele_swapped = swapped)
}
