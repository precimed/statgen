\name{statgen-package}
\alias{statgen-package}
\docType{package}
\title{Statistical Genetics Data Objects and Loaders}
\description{
CRAN-compatible R runtime for statgen statistical genetics data objects.
}
\details{
\pkg{statgen} provides R objects and loaders for reference-aligned statistical
genetics data: reference variants, genotypes, linkage disequilibrium (LD),
annotations, and GWAS summary statistics. The package is designed for analysis
workflows where genome build, chromosome names, and DNA strand have already
been harmonized upstream.

The central object is a \code{ReferencePanel}, which defines the ordered SNP
axis used by the other panel types. \code{LDPanel}, \code{GenotypePanel},
\code{AnnotationPanel}, and \code{Sumstats} objects expose data in those
reference coordinates. See \link{allele_matching} for how loading adjusts
allele orientation.

Typical workflows start by loading or caching a reference panel with
\code{\link{load_reference}}, then loading aligned data with functions such as
\code{\link{load_ld}}, \code{\link{load_genotype}},
\code{\link{load_annotations}}, or \code{\link{load_sumstats}}. Panel objects
provide accessors for aligned vectors and matrices, \code{\link{select_shards}}
for canonical chromosome subsetting, and cache helpers for repeated analyses.

LD distributions are directory-based runtime artifacts discovered through
\code{ld_manifest.json}; they are loaded with \code{\link{load_ld}} rather than
as individual shard files. R reads the Python \code{.npz} LD shard format
directly after \code{\link{prepare_ld_npz_for_r}} adds an R reference-cache
sidecar. R-native object caches use RDS files, while portable source inputs
such as PLINK BIM/BED/FAM, BED annotations, and summary-statistics TSV files
remain external inputs.
}
\seealso{
\link{allele_matching},
\code{\link{load_reference}}, \code{\link{load_ld}},
\code{\link{load_genotype}}, \code{\link{load_annotations}},
\code{\link{load_sumstats}}, \code{\link{version}},
\code{\link{set_verbosity}}
}
