.allele_hash_p <- 2147483647
.allele_hash_base1 <- 257
.allele_hash_base2 <- 263
.allele_hash_max_chars <- 1000000L
# At most 8192 terms per allele reduction, independent of batch memory.
.allele_hash_block_terms <- 8192L
# Raw batches use at most twice this target (1 MiB). Allow about 64 bytes
# per expanded base for raw data, products, indices and reduction buffers.
.allele_hash_batch_bytes <- 512L * 1024L
.hash64_shift <- bit64::as.integer64("4294967296")
.single_base_hash64 <- bit64::as.integer64(c(
  "1387274436937",
  "1395864371531",
  "1413044240719",
  "1468878815580"
))
names(.single_base_hash64) <- c("A", "C", "G", "T")

.warn_allele_cap <- function(a1, a2, path, rows) {
  lengths <- list(nchar(a1, type = "bytes"), nchar(a2, type = "bytes"))
  affected <- lapply(lengths, function(x) !is.na(x) & x > .allele_hash_max_chars)
  count <- sum(vapply(affected, sum, numeric(1)))
  if (count) {
    first <- vapply(affected, function(x) if (any(x)) which(x)[[1]] else Inf, numeric(1))
    col <- which.min(first)
    i <- first[[col]]
    warning(sprintf("%s: %d allele(s) exceed hash cap %d; first affected row %d, column a%d, length %d; hashing uses first %d bases",
                    path, count, .allele_hash_max_chars, rows[[i]], col, lengths[[col]][[i]], .allele_hash_max_chars), call. = FALSE)
  }
}

.allele_hash_batches <- function(lengths) {
  # Aggregate capped bytes can exceed R's integer range even though each
  # individual allele length and bounded batch index fit in an integer.
  ends <- cumsum(as.double(lengths))
  # Native partitioning is linear; findInterval would rescan its sorted input
  # on every call on supported R versions. Small batches use at most twice
  # the byte target; each allele reduction still has at most 8192 terms.
  limit <- min(.allele_hash_block_terms, .allele_hash_batch_bytes)
  long <- lengths > limit
  split <- c(TRUE, diff(floor((ends - 1) / .allele_hash_batch_bytes)) != 0) |
    long | c(FALSE, head(long, -1L))
  batch_starts <- which(split)
  batch_ends <- c(batch_starts[-1L] - 1L, length(lengths))
  list(starts = batch_starts, ends = batch_ends)
}

.allele_hash_powers <- function(length) {
  powers <- matrix(1, nrow = length + 1L, ncol = 2L)
  factors <- bit64::as.integer64(c(.allele_hash_base1, .allele_hash_base2))
  modulus <- bit64::as.integer64(.allele_hash_p)
  used <- 1L
  # PERF: logarithmic dispatch, linear native work. Residue products stay
  # below 2^62; convert only the reduced, exactly representable values to double.
  while (used < nrow(powers)) {
    count <- min(used, nrow(powers) - used)
    rows <- seq_len(count)
    for (lane in 1:2) {
      values <- bit64::as.integer64(powers[rows, lane])
      powers[used + rows, lane] <- as.double((values * factors[[lane]]) %% modulus)
    }
    factors <- (factors * factors) %% modulus
    used <- used + count
  }
  powers
}

.allele_hash64 <- function(alleles) {
  alleles <- as.character(alleles)
  if (!length(alleles)) return(bit64::integer64(0L))
  if (anyNA(alleles)) stop("Alleles must be non-missing character values", call. = FALSE)
  base_idx <- match(alleles, names(.single_base_hash64))
  single <- !is.na(base_idx)
  idx <- which(!single)
  if (!length(idx)) return(unname(.single_base_hash64)[base_idx])
  # Plain double vectors support in-place batch writes. integer64 subassignment
  # copies the whole output, so assemble the exact integer64 result only once.
  h1 <- numeric(length(alleles))
  h2 <- numeric(length(alleles))
  x <- c(66, 68, 72, 85)[base_idx[single]]  # ASCII byte + 1 for A/C/G/T.
  h1[single] <- .allele_hash_base1 + x
  h2[single] <- .allele_hash_base2 + x
  alleles <- alleles[idx]
  lengths <- nchar(alleles, type = "bytes")
  lengths <- pmin(lengths, .allele_hash_max_chars)
  block_terms <- min(.allele_hash_block_terms, .allele_hash_batch_bytes)
  max_len <- max(lengths)
  # Build once per operation, lazily after the single-base fast path.
  powers <- .allele_hash_powers(max_len)
  batches <- .allele_hash_batches(lengths)
  batch_starts <- batches$starts
  batch_ends <- batches$ends
  # PERF: bounded raw batches and rare long-allele blocks, never SNP-position scans.
  for (batch_i in seq_along(batch_starts)) {
    start <- batch_starts[[batch_i]]
    stop <- batch_ends[[batch_i]]
    lens <- lengths[start:stop]
    if (lens[[1]] > block_terms) {
      len <- lens[[1]]
      lanes <- powers[len + 1L, ]
      for (offset in seq.int(0L, len - 1L, by = block_terms)) {
        last <- min(len, offset + block_terms)
        x <- as.integer(charToRaw(substr(alleles[[start]], offset + 1L, last))) + 1
        positions <- len - seq.int(offset + 1L, last) + 1L
        for (lane in 1:2) {
          lanes[[lane]] <- (lanes[[lane]] + sum(x * powers[positions, lane])) %% .allele_hash_p
        }
      }
      h1[idx[[start]]] <- lanes[[1]]
      h2[idx[[start]]] <- lanes[[2]]
    } else {
      batch <- substr(alleles[start:stop], 1L, .allele_hash_max_chars)
      x <- as.integer(charToRaw(paste0(batch, collapse = ""))) + 1
      starts <- cumsum(lens) - lens
      groups <- rep.int(seq_along(lens), lens)
      positions <- rep.int(starts + lens, lens) - seq_along(x) + 1L
      lanes <- matrix(0, nrow = length(lens), ncol = 2L)
      for (lane in 1:2) {
        sums <- rowsum(x * powers[positions, lane], groups, reorder = FALSE)[, 1L]
        lanes[, lane] <- (sums + powers[lens + 1L, lane]) %% .allele_hash_p
      }
      h1[idx[start:stop]] <- lanes[, 1L]
      h2[idx[start:stop]] <- lanes[, 2L]
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
