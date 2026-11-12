.genotype_cache_schema <- "genotype_cache/0.2"
.bed_missing_int <- -1L
.haploid_modes <- c("raw", "ploidy_scaled")

num_sample <- function(x, ...) UseMethod("num_sample")
fid <- function(x, ...) UseMethod("fid")
iid <- function(x, ...) UseMethod("iid")
father_id <- function(x, ...) UseMethod("father_id")
mother_id <- function(x, ...) UseMethod("mother_id")
sex <- function(x, ...) UseMethod("sex")
is_male <- function(x, ...) UseMethod("is_male")
is_female <- function(x, ...) UseMethod("is_female")
ploidy_male <- function(x, ...) UseMethod("ploidy_male")
ploidy_female <- function(x, ...) UseMethod("ploidy_female")
source_layout <- function(x, ...) UseMethod("source_layout")
source_row0 <- function(x, ...) UseMethod("source_row0")
is_subject_present <- function(x, shard, ...) UseMethod("is_subject_present")
fetch_genotypes_int8 <- function(x, snp_indices, ...) UseMethod("fetch_genotypes_int8")
fetch_genotypes <- function(x, snp_indices, haploid_mode = NULL, ...) UseMethod("fetch_genotypes")

load_genotype <- function(bfile_prefix, reference, allow_indel_swaps = FALSE, match_by = "genomic_position") {
  allow_indel_swaps <- .validate_matching_options(match_by, allow_indel_swaps)
  prefix <- .validate_path_scalar(bfile_prefix, "bfile_prefix")
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  ref_shards <- shards(reference)
  if (!length(ref_shards)) {
    stop("load_genotype requires a reference with at least one shard", call. = FALSE)
  }
  ref_labels <- vapply(ref_shards, function(s) s$label, character(1))
  .validate_requested_shards(ref_labels, ref_labels, "load_genotype reference")

  if (grepl("@", prefix, fixed = TRUE)) {
    layout <- "sharded"
    sources <- stats::setNames(vector("list", length(ref_labels)), ref_labels)
    for (label in ref_labels) {
      sources[[label]] <- .load_genotype_source_record(gsub("@", label, prefix, fixed = TRUE), label)
    }
  } else {
    layout <- "non_sharded"
    shared <- .load_genotype_source_record(prefix, NULL)
    if ("X" %in% ref_labels && any(shared$bim$chr == "X") && !file.exists(paste0(prefix, ".ploidy"))) {
      warning(
        sprintf("%s: chrX genotype source has no .ploidy sidecar; defaulting chrX rows to male/female ploidy (1, 2)", prefix),
        call. = FALSE
      )
    }
    sources <- stats::setNames(rep(list(shared), length(ref_labels)), ref_labels)
  }

  subject_maps <- stats::setNames(vector("list", length(ref_labels)), ref_labels)
  if (identical(layout, "non_sharded")) {
    panel_fam <- sources[[1]]$fam
    for (label in ref_labels) {
      if (!any(sources[[label]]$bim$chr == label)) {
        warning(
          sprintf("%s: no source BIM rows for requested reference shard %s; marking shard absent", prefix, sQuote(label)),
          call. = FALSE
        )
      }
      subject_maps[[label]] <- .all_subjects_present(sources[[label]]$source_num_sample)
    }
  } else {
    autosomal <- ref_labels[ref_labels != "X"]
    if (length(autosomal)) {
      panel_fam <- sources[[autosomal[[1]]]]$fam
      if (length(autosomal) > 1L) {
        for (label in autosomal[-1L]) {
          if (!.fam_equal(sources[[label]]$fam, panel_fam)) {
            stop(sprintf("load_genotype: autosomal FAM mismatch for shard %s", sQuote(label)), call. = FALSE)
          }
        }
      }
      for (label in autosomal) {
        subject_maps[[label]] <- .all_subjects_present(sources[[label]]$source_num_sample)
      }
      if ("X" %in% ref_labels) {
        subject_maps[["X"]] <- .map_chrx_subjects(sources[["X"]]$fam, panel_fam, "load_genotype")
      }
    } else {
      panel_fam <- sources[["X"]]$fam
      subject_maps[["X"]] <- .all_subjects_present(sources[["X"]]$source_num_sample)
    }
  }

  id_matches <- vector("list", length(ref_shards))
  if (match_by == "id") {
    records <- if (layout == "non_sharded") list(shared) else unname(sources)
    inputs <- lapply(records, function(record) list(
      context = record$bim_path, ids = record$bim$snp,
      a1 = record$bim$a1_hash64, a2 = record$bim$a2_hash64, chr = record$bim$chr
    ))
    id_match <- .match_id_panel(reference, inputs, chromosome_scope = TRUE)
    source_offsets <- if (layout == "non_sharded") rep.int(0L, length(ref_shards)) else
      c(0, utils::head(cumsum(vapply(records, function(s) nrow(s$bim), integer(1))), -1L))
    offsets <- shard_offsets(reference)
    for (i in seq_along(ref_shards)) {
      ix <- seq.int(offsets$start0[i] + 1L, offsets$stop0[i])
      local <- id_match$source_rows[ix]
      local[!is.na(local)] <- local[!is.na(local)] - source_offsets[i]
      id_matches[[i]] <- list(source_rows = local, is_allele_swapped = id_match$is_allele_swapped[ix])
    }
  }

  out <- vector("list", length(ref_shards))
  for (i in seq_along(ref_shards)) {
    label <- ref_shards[[i]]$label
    out[[i]] <- .build_genotype_shard(ref_shards[[i]], sources[[label]], subject_maps[[label]], allow_indel_swaps, id_matches[[i]])
  }
  .new_genotype_panel(
    out,
    fid = panel_fam$fid,
    iid = panel_fam$iid,
    father_id = panel_fam$father_id,
    mother_id = panel_fam$mother_id,
    sex = panel_fam$sex,
    source_layout = layout
  )
}

save_genotype_cache <- function(panel, path) {
  if (!inherits(panel, "GenotypePanel")) {
    stop("panel must be a GenotypePanel", call. = FALSE)
  }
  path <- .validate_path_scalar(path, "path")
  offsets <- shard_offsets(panel)
  metadata <- list(
    schema = .genotype_cache_schema,
    n_shards = length(panel$shards),
    shard_labels = vapply(panel$shards, function(s) s$label, character(1)),
    shard_checksums = vapply(panel$shards, function(s) s$reference_checksum, character(1)),
    shard_start0 = offsets$start0,
    shard_stop0 = offsets$stop0,
    source_layout = panel$source_layout,
    bed_paths = vapply(panel$shards, function(s) s$bed_path, character(1)),
    bed_file_sizes = vapply(panel$shards, function(s) s$bed_file_size, numeric(1)),
    source_num_snp = vapply(panel$shards, function(s) s$source_num_snp, integer(1)),
    source_num_sample = vapply(panel$shards, function(s) s$source_num_sample, integer(1)),
    num_sample = panel$num_sample
  )
  payload <- list(
    metadata = metadata,
    is_present = is_present(panel),
    is_allele_swapped = is_allele_swapped(panel),
    ploidy_male = ploidy_male(panel),
    ploidy_female = ploidy_female(panel),
    source_row0 = source_row0(panel),
    subject_present = do.call(cbind, lapply(panel$shards, function(s) s$subject_present)),
    source_subject_row0 = do.call(cbind, lapply(panel$shards, function(s) s$source_subject_row0)),
    fid = fid(panel),
    iid = iid(panel),
    father_id = father_id(panel),
    mother_id = mother_id(panel),
    sex = sex(panel),
    is_male = is_male(panel),
    is_female = is_female(panel)
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, path)
  invisible(NULL)
}

load_genotype_cache <- function(path, shards = NULL, bed_path = NULL) {
  path <- .validate_path_scalar(path, "path")
  payload <- readRDS(path)
  .validate_genotype_cache_payload(payload)
  meta <- payload$metadata
  swap_flags <- if (meta$schema == .genotype_cache_schema) payload$is_allele_swapped else rep.int(FALSE, length(payload$is_present))
  selected <- .validate_requested_shards(shards, meta$shard_labels, "load_genotype_cache")
  label_to_index <- stats::setNames(seq_along(meta$shard_labels), meta$shard_labels)

  indices <- unname(label_to_index[selected])
  paths <- .resolve_bed_paths(meta, indices, bed_path)
  info <- file.info(paths)
  bad <- is.na(info$size) | info$isdir | info$size != meta$bed_file_sizes[indices]
  bad[is.na(bad)] <- TRUE
  if (any(bad)) {
    msg <- paste("Missing or incorrectly sized BED files:", paste(unique(paths[bad]), collapse = ", "))
    if (!is.null(bed_path)) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }

  out <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    label <- selected[[i]]
    idx <- label_to_index[[label]]
    start <- meta$shard_start0[[idx]] + 1L
    stop <- meta$shard_stop0[[idx]]
    out[[i]] <- .new_genotype_shard(
      label = label,
      bed_path = paths[[i]],
      bed_file_size = meta$bed_file_sizes[[idx]],
      source_num_snp = meta$source_num_snp[[idx]],
      source_num_sample = meta$source_num_sample[[idx]],
      source_row0 = payload$source_row0[start:stop],
      subject_present = payload$subject_present[, idx],
      source_subject_row0 = payload$source_subject_row0[, idx],
      is_present = payload$is_present[start:stop],
      is_allele_swapped = swap_flags[start:stop],
      ploidy_male = payload$ploidy_male[start:stop],
      ploidy_female = payload$ploidy_female[start:stop],
      reference_checksum = meta$shard_checksums[[idx]]
    )
  }
  .new_genotype_panel(
    out,
    fid = payload$fid,
    iid = payload$iid,
    father_id = payload$father_id,
    mother_id = payload$mother_id,
    sex = payload$sex,
    source_layout = meta$source_layout
  )
}

num_snp.GenotypeShard <- function(x, ...) x$num_snp
num_snp.GenotypePanel <- function(x, ...) x$num_snp
num_sample.GenotypePanel <- function(x, ...) x$num_sample
shards.GenotypePanel <- function(x, ...) x$shards
shard_offsets.GenotypePanel <- function(x, ...) x$shard_offsets
fid.GenotypePanel <- function(x, ...) x$fid
iid.GenotypePanel <- function(x, ...) x$iid
father_id.GenotypePanel <- function(x, ...) x$father_id
mother_id.GenotypePanel <- function(x, ...) x$mother_id
sex.GenotypePanel <- function(x, ...) x$sex
is_male.GenotypePanel <- function(x, ...) x$sex == 1L
is_female.GenotypePanel <- function(x, ...) x$sex == 2L
is_allele_swapped <- function(x, ...) UseMethod("is_allele_swapped")
is_allele_swapped.GenotypeShard <- function(x, ...) x$is_allele_swapped
is_allele_swapped.GenotypePanel <- function(x, ...) unlist(lapply(x$shards, is_allele_swapped), use.names = FALSE)
is_present.GenotypeShard <- function(x, ...) x$is_present
is_present.GenotypePanel <- function(x, ...) unlist(lapply(x$shards, is_present), use.names = FALSE)
ploidy_male.GenotypeShard <- function(x, ...) x$ploidy_male
ploidy_male.GenotypePanel <- function(x, ...) unlist(lapply(x$shards, ploidy_male), use.names = FALSE)
ploidy_female.GenotypeShard <- function(x, ...) x$ploidy_female
ploidy_female.GenotypePanel <- function(x, ...) unlist(lapply(x$shards, ploidy_female), use.names = FALSE)
source_layout.GenotypePanel <- function(x, ...) x$source_layout
source_row0.GenotypeShard <- function(x, ...) x$source_row0
source_row0.GenotypePanel <- function(x, ...) unlist(lapply(x$shards, source_row0), use.names = FALSE)

is_subject_present.GenotypePanel <- function(x, shard, ...) {
  label <- .validate_char_scalar(shard, "shard")
  for (s in x$shards) {
    if (identical(s$label, label)) {
      return(s$subject_present)
    }
  }
  stop(sprintf("GenotypePanel.is_subject_present: requested shard %s is not loaded", sQuote(label)), call. = FALSE)
}

select_shards.GenotypePanel <- function(x, shards, ...) {
  available <- vapply(x$shards, function(s) s$label, character(1))
  selected <- .validate_requested_shards(shards, available, "GenotypePanel.select_shards")
  by_label <- stats::setNames(x$shards, available)
  .new_genotype_panel(
    unname(by_label[selected]),
    fid = x$fid,
    iid = x$iid,
    father_id = x$father_id,
    mother_id = x$mother_id,
    sex = x$sex,
    source_layout = x$source_layout
  )
}

fetch_genotypes_int8.GenotypePanel <- function(x, snp_indices, ...) {
  if (length(list(...))) stop("Unexpected fetch arguments", call. = FALSE)
  indices <- .normalize_genotype_indices(snp_indices, x$num_snp)
  out <- matrix(.bed_missing_int, nrow = x$num_sample, ncol = length(indices))
  if (!length(indices)) {
    storage.mode(out) <- "integer"
    return(out)
  }
  addressed <- .addressed_genotype_shards(x, indices)
  .validate_requested_genotype_present(x, addressed)

  for (entry in addressed) {
    shard <- entry$shard
    effective_bed <- shard$bed_path
    .validate_bed_file(effective_bed, shard$source_num_sample, shard$source_num_snp)
    source_rows <- shard$source_row0[entry$local_indices]
    source_geno <- .read_bed_rows_int8_codes(effective_bed, source_rows, shard$source_num_sample)
    swapped_columns <- which(shard$is_allele_swapped[entry$local_indices])
    if (length(swapped_columns)) {
      calls <- source_geno[, swapped_columns, drop = FALSE]
      flip_lut <- c(-1L, 2L, 1L, 0L)
      source_geno[, swapped_columns] <- matrix(flip_lut[calls + 2L], nrow = nrow(calls), ncol = ncol(calls))
    }
    panel_rows <- which(shard$subject_present)
    if (length(panel_rows)) {
      source_subject_rows <- shard$source_subject_row0[panel_rows] + 1L
      out[panel_rows, entry$cols] <- source_geno[source_subject_rows, , drop = FALSE]
    }
  }
  storage.mode(out) <- "integer"
  out
}

fetch_genotypes.GenotypePanel <- function(x, snp_indices, haploid_mode = NULL, ...) {
  if (length(list(...))) stop("Unexpected fetch arguments", call. = FALSE)
  mode <- if (is.null(haploid_mode)) "raw" else as.character(haploid_mode)
  if (length(mode) != 1L || is.na(mode) || !(mode %in% .haploid_modes)) {
    stop("GenotypePanel.fetch_genotypes: haploid_mode must be 'raw' or 'ploidy_scaled'", call. = FALSE)
  }
  indices <- .normalize_genotype_indices(snp_indices, x$num_snp)
  gi <- fetch_genotypes_int8(x, indices)
  out <- matrix(as.numeric(gi), nrow = nrow(gi), ncol = ncol(gi))
  out[gi == .bed_missing_int] <- NaN
  if (identical(mode, "ploidy_scaled")) {
    out <- .apply_ploidy_scaled(x, out, indices)
  }
  out
}

save_cache.GenotypePanel <- function(x, path, ...) save_genotype_cache(x, path)

print.GenotypePanel <- function(x, ...) {
  labels <- vapply(x$shards, function(s) s$label, character(1))
  cat(sprintf(
    "<GenotypePanel: %d SNPs x %d samples across %d shard(s): %s>\n",
    x$num_snp, x$num_sample, length(labels), paste(labels, collapse = ", ")
  ))
  invisible(x)
}

.new_genotype_shard <- function(label, bed_path, bed_file_size, source_num_snp,
                                source_num_sample, source_row0, subject_present,
                                source_subject_row0, is_present, is_allele_swapped,
                                ploidy_male, ploidy_female, reference_checksum) {
  label <- as.character(label)
  source_num_snp <- as.integer(source_num_snp)
  source_num_sample <- as.integer(source_num_sample)
  bed_file_size <- as.numeric(bed_file_size)
  if (is.na(source_num_snp) || is.na(source_num_sample) || source_num_snp < 0L || source_num_sample < 0L) {
    stop("source_num_snp and source_num_sample must be non-negative", call. = FALSE)
  }
  expected <- .expected_bed_size(source_num_sample, source_num_snp)
  if (!identical(as.numeric(expected), bed_file_size)) {
    stop(sprintf("Invalid genotype metadata for shard %s: bed_file_size mismatch", label), call. = FALSE)
  }

  source_row0 <- as.integer(source_row0)
  is_present <- as.logical(is_present)
  ploidy_male <- as.numeric(ploidy_male)
  ploidy_female <- as.numeric(ploidy_female)
  .validate_equal_lengths(
    list(source_row0 = source_row0, is_present = is_present, is_allele_swapped = is_allele_swapped, ploidy_male = ploidy_male, ploidy_female = ploidy_female),
    sprintf("GenotypeShard %s: SNP-axis vector lengths must match", label)
  )
  if (!identical(is_present, source_row0 >= 0L)) {
    stop(sprintf("GenotypeShard %s: is_present must match source_row0 >= 0", label), call. = FALSE)
  }
  bad_source <- source_row0 < -1L | source_row0 >= source_num_snp
  if (any(bad_source)) {
    stop(sprintf("GenotypeShard %s: source_row0 out of source BIM bounds", label), call. = FALSE)
  }
  if (!is.logical(is_allele_swapped) || anyNA(is_allele_swapped) || any(is_allele_swapped & !is_present)) {
    stop(sprintf("GenotypeShard %s: is_allele_swapped must be logical and requires is_present", label), call. = FALSE)
  }
  absent <- !is_present
  if (any(!is.nan(ploidy_male[absent])) || any(!is.nan(ploidy_female[absent]))) {
    stop(sprintf("GenotypeShard %s: absent SNPs must have NaN ploidy", label), call. = FALSE)
  }
  present_ploidy <- c(ploidy_male[is_present], ploidy_female[is_present])
  if (length(present_ploidy) && any(!(present_ploidy %in% c(0, 1, 2)))) {
    stop(sprintf("GenotypeShard %s: present ploidy values must be 0, 1, or 2", label), call. = FALSE)
  }

  subject_present <- as.logical(subject_present)
  source_subject_row0 <- as.integer(source_subject_row0)
  .validate_equal_lengths(
    list(subject_present = subject_present, source_subject_row0 = source_subject_row0),
    sprintf("GenotypeShard %s: sample-axis vector lengths must match", label)
  )
  if (!identical(subject_present, source_subject_row0 >= 0L)) {
    stop(sprintf("GenotypeShard %s: subject_present must match source_subject_row0 >= 0", label), call. = FALSE)
  }
  bad_subject <- source_subject_row0 < -1L | source_subject_row0 >= source_num_sample
  if (any(bad_subject)) {
    stop(sprintf("GenotypeShard %s: source_subject_row0 out of source FAM bounds", label), call. = FALSE)
  }

  structure(
    list(
      chr = label,
      label = label,
      bed_path = .validate_path_scalar(bed_path, "bed_path"),
      bed_file_size = bed_file_size,
      source_num_snp = source_num_snp,
      source_num_sample = source_num_sample,
      source_row0 = source_row0,
      subject_present = subject_present,
      source_subject_row0 = source_subject_row0,
      is_present = is_present,
      is_allele_swapped = is_allele_swapped,
      ploidy_male = ploidy_male,
      ploidy_female = ploidy_female,
      reference_checksum = as.character(reference_checksum),
      num_snp = length(is_present)
    ),
    class = "GenotypeShard"
  )
}

.new_genotype_panel <- function(shard_list, fid, iid, father_id, mother_id, sex, source_layout) {
  if (!length(shard_list)) {
    stop("GenotypePanel requires at least one shard", call. = FALSE)
  }
  if (!(source_layout %in% c("non_sharded", "sharded"))) {
    stop("source_layout must be 'non_sharded' or 'sharded'", call. = FALSE)
  }
  fid <- as.character(fid)
  iid <- as.character(iid)
  father_id <- as.character(father_id)
  mother_id <- as.character(mother_id)
  sex <- as.integer(sex)
  .validate_equal_lengths(
    list(fid = fid, iid = iid, father_id = father_id, mother_id = mother_id, sex = sex),
    "GenotypePanel FAM vector lengths must match"
  )
  if (any(!(sex %in% c(0L, 1L, 2L)))) {
    stop("GenotypePanel sex values must be 0, 1, or 2", call. = FALSE)
  }
  num_sample <- length(fid)
  pos <- 0L
  offsets <- data.frame(
    shard_label = character(length(shard_list)),
    start0 = integer(length(shard_list)),
    stop0 = integer(length(shard_list)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(shard_list)) {
    if (length(shard_list[[i]]$subject_present) != num_sample) {
      stop(sprintf("Genotype shard %s: subject_present length does not match panel sample axis", shard_list[[i]]$label), call. = FALSE)
    }
    offsets$shard_label[[i]] <- shard_list[[i]]$label
    offsets$start0[[i]] <- pos
    pos <- pos + shard_list[[i]]$num_snp
    offsets$stop0[[i]] <- pos
  }
  structure(
    list(
      shards = shard_list,
      fid = fid,
      iid = iid,
      father_id = father_id,
      mother_id = mother_id,
      sex = sex,
      source_layout = source_layout,
      num_sample = num_sample,
      num_snp = pos,
      shard_offsets = offsets
    ),
    class = "GenotypePanel"
  )
}

.load_genotype_source_record <- function(prefix, label) {
  paths <- .require_genotype_source_paths(prefix, label)
  bim <- .parse_bim(paths$bim)
  source_num_snp <- nrow(bim)
  bim$line <- seq_len(source_num_snp)
  bim$source_row0 <- if (source_num_snp) seq.int(0L, source_num_snp - 1L) else integer()
  if (!is.null(label) && any(bim$chr != label)) {
    idx <- which(bim$chr != label)[[1]]
    stop(
      sprintf("%s:%d: sharded genotype BIM for %s contains chr %s", paths$bim, idx, sQuote(label), sQuote(bim$chr[[idx]])),
      call. = FALSE
    )
  }
  fam <- .parse_fam(paths$fam)
  if (identical(label, "X") && is.null(paths$ploidy)) {
    warning(
      sprintf("%s: chrX genotype source has no .ploidy sidecar; defaulting chrX rows to male/female ploidy (1, 2)", prefix),
      call. = FALSE
    )
  }
  ploidy <- .parse_ploidy(paths$ploidy, source_num_snp, bim$chr)
  bed_size <- .validate_bed_file(paths$bed, nrow(fam), source_num_snp)
  list(
    bed_path = paths$bed,
    bim_path = paths$bim,
    bed_file_size = bed_size,
    source_num_snp = source_num_snp,
    source_num_sample = nrow(fam),
    bim = bim,
    fam = fam,
    ploidy_male = ploidy$male,
    ploidy_female = ploidy$female
  )
}

.require_genotype_source_paths <- function(prefix, label) {
  out <- list(
    bed = paste0(prefix, ".bed"),
    bim = paste0(prefix, ".bim"),
    fam = paste0(prefix, ".fam")
  )
  for (path in out) {
    if (!file.exists(path)) {
      where <- if (is.null(label)) "" else sprintf(" for requested shard %s", sQuote(label))
      stop(sprintf("Missing genotype source%s: %s", where, path), call. = FALSE)
    }
  }
  ploidy <- paste0(prefix, ".ploidy")
  out$ploidy <- if (file.exists(ploidy)) ploidy else NULL
  out
}

.all_subjects_present <- function(source_num_sample) {
  list(
    subject_present = rep.int(TRUE, source_num_sample),
    source_subject_row0 = if (source_num_sample) seq.int(0L, source_num_sample - 1L) else integer()
  )
}

.fam_equal <- function(lhs, rhs) {
  identical(lhs[c("fid", "iid", "father_id", "mother_id", "sex")], rhs[c("fid", "iid", "father_id", "mother_id", "sex")])
}

.map_chrx_subjects <- function(chrx_fam, panel_fam, where) {
  if (.fam_equal(chrx_fam, panel_fam)) {
    return(.all_subjects_present(nrow(panel_fam)))
  }
  panel_key <- paste(panel_fam$fid, panel_fam$iid, sep = "\r")
  src_key <- paste(chrx_fam$fid, chrx_fam$iid, sep = "\r")
  panel_row <- match(src_key, panel_key)
  missing <- is.na(panel_row)
  if (any(missing)) {
    idx <- which(missing)[[1]]
    stop(sprintf("%s: chrX FAM subject (%s, %s) is absent from the autosomal sample axis", where, sQuote(chrx_fam$fid[[idx]]), sQuote(chrx_fam$iid[[idx]])), call. = FALSE)
  }
  for (field in c("father_id", "mother_id", "sex")) {
    bad <- chrx_fam[[field]] != panel_fam[[field]][panel_row]
    if (any(bad)) {
      idx <- which(bad)[[1]]
      stop(sprintf("%s: chrX FAM metadata mismatch for subject (%s, %s)", where, sQuote(chrx_fam$fid[[idx]]), sQuote(chrx_fam$iid[[idx]])), call. = FALSE)
    }
  }
  subject_present <- rep.int(FALSE, nrow(panel_fam))
  source_subject_row0 <- rep.int(-1L, nrow(panel_fam))
  subject_present[panel_row] <- TRUE
  source_subject_row0[panel_row] <- seq.int(0L, nrow(chrx_fam) - 1L)
  list(subject_present = subject_present, source_subject_row0 = source_subject_row0)
}

.build_genotype_shard <- function(ref_shard, source, subject_map, allow_indel_swaps, id_match = NULL) {
  if (!is.null(id_match)) {
    source_bim <- source$bim
    local_match <- id_match$source_rows
    is_allele_swapped <- id_match$is_allele_swapped
  } else {
    is_allele_swapped <- rep.int(FALSE, ref_shard$num_snp)
    source_bim <- source$bim[source$bim$chr == ref_shard$label, , drop = FALSE]
    if (!nrow(source_bim)) {
      local_match <- rep.int(NA_integer_, ref_shard$num_snp)
    } else {
      match <- .match_variant_keys(
        bp(ref_shard), a1_hash64(ref_shard), a2_hash64(ref_shard),
        source_bim$bp, source_bim$a1_hash64, source_bim$a2_hash64,
        ref_shard$label, "genotype", allow_indel_swaps
      )
      local_match <- match$source_rows
      is_allele_swapped <- match$is_allele_swapped
      .report_variant_swaps(source$bim_path, ref_shard$label, sum(is_allele_swapped),
                            match$num_strand_ambiguous, match$num_disabled_indels)
    }
  }

  is_present <- !is.na(local_match)
  source_row0_vec <- rep.int(-1L, ref_shard$num_snp)
  ploidy_male_vec <- rep.int(NaN, ref_shard$num_snp)
  ploidy_female_vec <- rep.int(NaN, ref_shard$num_snp)
  if (any(is_present)) {
    source_rows <- source_bim$source_row0[local_match[is_present]]
    source_row0_vec[is_present] <- source_rows
    ploidy_male_vec[is_present] <- source$ploidy_male[source_rows + 1L]
    ploidy_female_vec[is_present] <- source$ploidy_female[source_rows + 1L]
  }

  .new_genotype_shard(
    label = ref_shard$label,
    bed_path = source$bed_path,
    bed_file_size = source$bed_file_size,
    source_num_snp = source$source_num_snp,
    source_num_sample = source$source_num_sample,
    source_row0 = source_row0_vec,
    subject_present = subject_map$subject_present,
    source_subject_row0 = subject_map$source_subject_row0,
    is_present = is_present,
    is_allele_swapped = is_allele_swapped,
    ploidy_male = ploidy_male_vec,
    ploidy_female = ploidy_female_vec,
    reference_checksum = ref_shard$checksum
  )
}

.normalize_genotype_indices <- function(snp_indices, num_snp_value) {
  if (is.null(snp_indices)) {
    stop("GenotypePanel.fetch_genotypes: snp_indices must be integer indices", call. = FALSE)
  }
  if (!length(snp_indices)) {
    return(integer())
  }
  arr <- as.numeric(snp_indices)
  if (length(arr) != length(snp_indices) || any(is.na(arr) | !is.finite(arr) | floor(arr) != arr)) {
    stop("GenotypePanel.fetch_genotypes: snp_indices must be integer indices", call. = FALSE)
  }
  arr <- as.integer(arr)
  bad <- arr < 1L | arr > num_snp_value
  if (any(bad)) {
    stop(
      sprintf("GenotypePanel.fetch_genotypes: SNP index %d is out of bounds for num_snp=%d", arr[which(bad)[[1]]], num_snp_value),
      call. = FALSE
    )
  }
  arr
}

.addressed_genotype_shards <- function(panel, indices) {
  out <- list()
  offsets <- shard_offsets(panel)
  for (i in seq_along(panel$shards)) {
    start <- offsets$start0[[i]] + 1L
    stop <- offsets$stop0[[i]]
    cols <- which(indices >= start & indices <= stop)
    if (length(cols)) {
      out[[length(out) + 1L]] <- list(
        shard = panel$shards[[i]],
        cols = cols,
        local_indices = indices[cols] - offsets$start0[[i]]
      )
    }
  }
  out
}

.resolve_bed_paths <- function(meta, indices, bed_path) {
  if (is.null(bed_path)) {
    return(meta$bed_paths[indices])
  }
  labels <- meta$shard_labels[indices]
  override <- .validate_path_scalar(bed_path, "bed_path")
  if (grepl("@", override, fixed = TRUE)) {
    if (identical(meta$source_layout, "non_sharded")) {
      stop("load_genotype_cache: @ override incompatible with non-sharded panel metadata", call. = FALSE)
    }
    return(vapply(labels, function(label) gsub("@", label, override, fixed = TRUE), character(1)))
  }
  if (identical(meta$source_layout, "sharded") && length(indices)) {
    nsnp <- unique(meta$source_num_snp[indices])
    nsample <- unique(meta$source_num_sample[indices])
    if (length(nsnp) != 1L || length(nsample) != 1L) {
      stop("load_genotype_cache: flat bed_path override requires equal source_num_snp and source_num_sample across selected shards", call. = FALSE)
    }
  }
  rep(override, length(indices))
}

.validate_requested_genotype_present <- function(panel, addressed) {
  start_by_label <- stats::setNames(shard_offsets(panel)$start0, shard_offsets(panel)$shard_label)
  for (entry in addressed) {
    present <- entry$shard$is_present[entry$local_indices]
    if (!all(present)) {
      j <- which(!present)[[1]]
      panel_index <- start_by_label[[entry$shard$label]] + entry$local_indices[[j]]
      stop(
        sprintf("GenotypePanel.fetch_genotypes: requested SNP %d in shard %s is not present in the genotype source", panel_index, sQuote(entry$shard$label)),
        call. = FALSE
      )
    }
  }
}

.apply_ploidy_scaled <- function(panel, geno, indices) {
  if (!length(indices)) {
    return(geno)
  }
  unknown <- panel$sex == 0L
  male <- panel$sex == 1L
  female <- panel$sex == 2L
  addressed <- .addressed_genotype_shards(panel, indices)
  for (entry in addressed) {
    shard <- entry$shard
    pm <- shard$ploidy_male[entry$local_indices]
    pf <- shard$ploidy_female[entry$local_indices]
    differs <- pm != pf
    unknown_present <- unknown & shard$subject_present
    if (any(unknown_present) && any(differs) && any(is.finite(geno[unknown_present, entry$cols[differs], drop = FALSE]))) {
      stop("GenotypePanel.fetch_genotypes: haploid_mode='ploidy_scaled' requires known FAM sex when male and female ploidy differ", call. = FALSE)
    }
    if (any(male)) {
      geno[male, entry$cols] <- sweep(geno[male, entry$cols, drop = FALSE], 2L, pm / 2, `*`)
    }
    if (any(female)) {
      geno[female, entry$cols] <- sweep(geno[female, entry$cols, drop = FALSE], 2L, pf / 2, `*`)
    }
    same_ploidy <- !differs
    if (any(unknown) && any(same_ploidy)) {
      cols <- entry$cols[same_ploidy]
      geno[unknown, cols] <- sweep(geno[unknown, cols, drop = FALSE], 2L, pf[same_ploidy] / 2, `*`)
    }
  }
  geno
}

.validate_genotype_cache_payload <- function(payload) {
  schema <- payload$metadata$schema
  if (length(schema) != 1L || is.na(schema) || !(schema %in% c(.genotype_cache_schema, "genotype_cache/0.1"))) {
    stop("Unsupported genotype cache schema", call. = FALSE)
  }
  cache <- .validate_cache_payload_metadata(
    payload,
    schema,
    "genotype",
    extra_fields = c("bed_paths", "bed_file_sizes", "source_num_snp", "source_num_sample")
  )
  meta <- cache$meta
  n_shards <- cache$n_shards
  total <- cache$total
  .validate_requested_shards(meta$shard_labels, meta$shard_labels, "load_genotype_cache")
  if (!(meta$source_layout %in% c("non_sharded", "sharded"))) {
    stop("Invalid genotype cache: source_layout must be 'non_sharded' or 'sharded'", call. = FALSE)
  }
  for (field in c("is_present", "ploidy_male", "ploidy_female", "source_row0", "subject_present", "source_subject_row0", "fid", "iid", "father_id", "mother_id", "sex", "is_male", "is_female")) {
    if (is.null(payload[[field]])) {
      stop(sprintf("Invalid genotype cache: missing field %s", sQuote(field)), call. = FALSE)
    }
  }
  if (length(payload$is_present) != total || length(payload$source_row0) != total) {
    stop("Invalid genotype cache: SNP-axis vector lengths mismatch", call. = FALSE)
  }
  if (length(payload$ploidy_male) != total || length(payload$ploidy_female) != total) {
    stop("Invalid genotype cache: ploidy vector lengths mismatch", call. = FALSE)
  }
  if (schema == .genotype_cache_schema) {
    flags <- payload$is_allele_swapped
    if (!is.logical(flags) || !is.null(dim(flags)) || length(flags) != total || anyNA(flags)) {
      stop("Invalid genotype cache: is_allele_swapped must be a logical SNP-axis vector", call. = FALSE)
    }
    if (any(flags & !payload$is_present)) {
      stop("Invalid genotype cache: is_allele_swapped requires is_present", call. = FALSE)
    }
  }
  num_sample_value <- as.integer(meta$num_sample)
  if (is.na(num_sample_value) || num_sample_value < 0L) {
    stop("Invalid genotype cache: num_sample must be non-negative", call. = FALSE)
  }
  for (field in c("fid", "iid", "father_id", "mother_id", "sex", "is_male", "is_female")) {
    if (length(payload[[field]]) != num_sample_value) {
      stop(sprintf("Invalid genotype cache: %s length mismatch", field), call. = FALSE)
    }
  }
  sex_vec <- as.integer(payload$sex)
  if (any(!(sex_vec %in% c(0L, 1L, 2L)))) {
    stop("Invalid genotype cache: sex values must be 0, 1, or 2", call. = FALSE)
  }
  if (!identical(as.logical(payload$is_male), sex_vec == 1L)) {
    stop("Invalid genotype cache: is_male does not match sex", call. = FALSE)
  }
  if (!identical(as.logical(payload$is_female), sex_vec == 2L)) {
    stop("Invalid genotype cache: is_female does not match sex", call. = FALSE)
  }
  if (!identical(as.integer(dim(payload$subject_present)), c(num_sample_value, n_shards))) {
    stop("Invalid genotype cache: subject_present shape mismatch", call. = FALSE)
  }
  if (!identical(as.integer(dim(payload$source_subject_row0)), c(num_sample_value, n_shards))) {
    stop("Invalid genotype cache: source_subject_row0 shape mismatch", call. = FALSE)
  }
  for (i in seq_len(n_shards)) {
    expected <- .expected_bed_size(meta$source_num_sample[[i]], meta$source_num_snp[[i]])
    if (!identical(as.numeric(expected), as.numeric(meta$bed_file_sizes[[i]]))) {
      stop(sprintf("Invalid genotype cache: bed_file_size mismatch for shard %s", sQuote(meta$shard_labels[[i]])), call. = FALSE)
    }
  }
  invisible(NULL)
}
