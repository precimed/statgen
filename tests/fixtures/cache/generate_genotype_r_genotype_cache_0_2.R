#!/usr/bin/env Rscript
# Run from the repository root; uses shared swap inputs with indel swaps enabled.
statgen_cache_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
source(file.path(statgen_cache_dir, "cache_fixture_helpers.R"))
statgen_generate_r_cache_fixture("genotype")
