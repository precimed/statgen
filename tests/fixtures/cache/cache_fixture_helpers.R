statgen_repo_root <- function() {
  normalizePath(file.path(statgen_cache_dir, "..", "..", ".."), mustWork = TRUE)
}

statgen_source_r_package <- function(repo_root) {
  r_files <- list.files(file.path(repo_root, "R-package", "R"), pattern = "[.]R$", full.names = TRUE)
  for (path in r_files) {
    source(path)
  }
  set_verbosity("quiet")
}

statgen_cache_fixture_reference <- function(repo_root) {
  load_reference(file.path("tests", "fixtures", "reference", "sharded", "@.bim"))
}

statgen_cache_fixture_annotations <- function(repo_root, reference) {
  fixture_dir <- file.path("tests", "fixtures")
  binary <- load_annotations(
    file.path(fixture_dir, "annotations", c("anno1.bed", "anno2.bed")),
    reference
  )
  continuous <- load_annotation(
    file.path(fixture_dir, "annotations", "continuous.annot"),
    reference,
    has_header = TRUE,
    value_columns = c("score", "weight"),
    annotation_metadata_path = file.path(fixture_dir, "annotations", "continuous.meta")
  )
  union_annotations(binary, continuous)
}

statgen_generate_r_cache_fixture <- function(object_name) {
  repo_root <- statgen_repo_root()
  cache_dir <- file.path(repo_root, "tests", "fixtures", "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  statgen_source_r_package(repo_root)
  reference <- statgen_cache_fixture_reference(repo_root)
  fixture_dir <- file.path("tests", "fixtures")

  if (identical(object_name, "reference")) {
    save_reference_cache(reference, file.path(cache_dir, "reference_r_reference_cache_0_1.rds"))
  } else if (identical(object_name, "sumstats")) {
    panel <- load_sumstats(file.path(fixture_dir, "sumstats", "traits_complete.tsv.gz"), reference)
    save_sumstats_cache(panel, file.path(cache_dir, "sumstats_r_sumstats_cache_0_1.rds"))
  } else if (identical(object_name, "annotations")) {
    save_annotations_cache(
      statgen_cache_fixture_annotations(repo_root, reference),
      file.path(cache_dir, "annotations_r_annotations_cache_0_2.rds")
    )
  } else if (identical(object_name, "genotype")) {
    root <- file.path("tests", "fixtures", "allele_swaps")
    reference <- load_reference(file.path(root, "reference.bim"))
    panel <- load_genotype(file.path(root, "genotype", "@"), reference, allow_indel_swaps = TRUE)
    save_genotype_cache(panel, file.path(cache_dir, "genotype_r_genotype_cache_0_2.rds"))
  } else {
    stop(sprintf("unknown cache fixture object: %s", object_name), call. = FALSE)
  }
}
