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

.validate_allow_indel_swaps <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("allow_indel_swaps must be a logical scalar", call. = FALSE)
  }
  value
}

.report_variant_swaps <- function(path, label, swaps, strand_ambiguous, disabled_indels) {
  if (get_verbosity() != "info") return(invisible(NULL))
  context <- sprintf("%s: shard %s: ", path, label)
  if (swaps > 0L) cat(sprintf("%s%d successful allele swap(s)\n", context, swaps), file = stderr())
  if (strand_ambiguous > 0L) cat(sprintf("%s%d strand-ambiguous A/T or C/G swap(s)\n", context, strand_ambiguous), file = stderr())
  if (disabled_indels > 0L) cat(sprintf("%s%d indel row(s) unmatched because indel swaps are disabled\n", context, disabled_indels), file = stderr())
  invisible(NULL)
}
