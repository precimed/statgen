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
