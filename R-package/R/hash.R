.allele_hash_p <- 2147483647
.allele_hash_base1 <- 257
.allele_hash_base2 <- 263
.allele_hash_max_chars <- 150L
.hash64_shift <- bit64::as.integer64("4294967296")
.single_base_hash64 <- bit64::as.integer64(c(
  "1387274436937",
  "1395864371531",
  "1413044240719",
  "1468878815580"
))
names(.single_base_hash64) <- c("A", "C", "G", "T")

.allele_hash64 <- function(alleles) {
  alleles <- as.character(alleles)
  if (!length(alleles)) {
    return(bit64::integer64(0L))
  }
  if (anyNA(alleles)) {
    stop("Alleles must be non-missing character values", call. = FALSE)
  }

  lengths <- nchar(alleles, type = "chars", allowNA = FALSE)
  long <- lengths > .allele_hash_max_chars
  if (any(long)) {
    warning("Allele length exceeds 150 characters; hashing uses first 150 characters", call. = FALSE)
    alleles[long] <- substr(alleles[long], 1L, .allele_hash_max_chars)
    lengths <- pmin(lengths, .allele_hash_max_chars)
  }

  out <- bit64::integer64(length(alleles))
  single_base <- lengths == 1L
  if (any(single_base)) {
    single_idx <- which(single_base)
    assigned <- rep(FALSE, length(single_idx))
    single_alleles <- alleles[single_idx]
    for (base in names(.single_base_hash64)) {
      matches <- single_alleles == base
      if (any(matches)) {
        out[single_idx[matches]] <- .single_base_hash64[[base]]
        assigned[matches] <- TRUE
      }
    }
    if (any(!assigned)) {
      other_idx <- single_idx[!assigned]
      out[other_idx] <- .allele_hash64_utf8(alleles[other_idx])
    }
  }
  if (all(single_base)) {
    return(out)
  }

  multi_idx <- which(!single_base)
  out[multi_idx] <- .allele_hash64_utf8(alleles[multi_idx])
  out
}

.allele_hash64_utf8 <- function(alleles) {
  h1 <- rep.int(1, length(alleles))
  h2 <- rep.int(1, length(alleles))
  bytes <- lapply(alleles, charToRaw)
  byte_lengths <- lengths(bytes)
  max_len <- max(byte_lengths)
  if (max_len > 0L) {
    for (j in seq_len(max_len)) {
      active <- byte_lengths >= j
      byte <- vapply(bytes[active], function(x) as.integer(x[[j]]), integer(1))
      x <- byte + 1
      h1[active] <- (h1[active] * .allele_hash_base1 + x) %% .allele_hash_p
      h2[active] <- (h2[active] * .allele_hash_base2 + x) %% .allele_hash_p
    }
  }
  bit64::as.integer64(h1) * .hash64_shift + bit64::as.integer64(h2)
}

.is_single_base_hash <- function(hash) {
  hash %in% .single_base_hash64
}

.is_single_nucleotide_variant_hash <- function(a1_hash, a2_hash) {
  .is_single_base_hash(a1_hash) & .is_single_base_hash(a2_hash)
}

.is_strand_ambiguous_hash <- function(a1_hash, a2_hash) {
  h <- .single_base_hash64
  (a1_hash == h[["A"]] & a2_hash == h[["T"]]) |
    (a1_hash == h[["T"]] & a2_hash == h[["A"]]) |
    (a1_hash == h[["C"]] & a2_hash == h[["G"]]) |
    (a1_hash == h[["G"]] & a2_hash == h[["C"]])
}
