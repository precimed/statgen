.statgen_state <- new.env(parent = emptyenv())
.statgen_state$verbosity <- "info"

set_verbosity <- function(level) {
  if (!is.character(level) || length(level) != 1L || is.na(level)) {
    stop("level must be a character scalar", call. = FALSE)
  }
  if (!(level %in% c("quiet", "info"))) {
    stop("level must be one of 'quiet' or 'info'", call. = FALSE)
  }
  .statgen_state$verbosity <- level
  invisible(level)
}

get_verbosity <- function() {
  .statgen_state$verbosity
}

.info <- function(message) {
  if (get_verbosity() == "info") cat(message, "\n", sep = "", file = stderr())
  invisible(NULL)
}

.validate_allow_indel_swaps <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("allow_indel_swaps must be a logical scalar", call. = FALSE)
  }
  value
}

.report_variant_swaps <- function(path, label, swaps, strand_ambiguous, disabled_indels) {
  context <- sprintf("%s: shard %s: ", path, label)
  if (swaps > 0L) .info(sprintf("%s%d successful allele swap(s)", context, swaps))
  if (strand_ambiguous > 0L) .info(sprintf("%s%d strand-ambiguous A/T or C/G swap(s)", context, strand_ambiguous))
  if (disabled_indels > 0L) .info(sprintf("%s%d indel row(s) unmatched because indel swaps are disabled", context, disabled_indels))
  invisible(NULL)
}

.validate_matching_options <- function(match_by, allow_indel_swaps) {
  allow_indel_swaps <- .validate_allow_indel_swaps(allow_indel_swaps)
  if (!is.character(match_by) || length(match_by) != 1L || is.na(match_by) ||
      !(match_by %in% c("genomic_position", "id"))) {
    stop('match_by must be "genomic_position" or "id"', call. = FALSE)
  }
  if (match_by == "id" && allow_indel_swaps) {
    stop('allow_indel_swaps=TRUE is incompatible with match_by="id"', call. = FALSE)
  }
  allow_indel_swaps
}
