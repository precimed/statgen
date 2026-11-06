\name{sumstats}
\alias{load_sumstats}
\alias{create_sumstats}
\alias{save_sumstats_cache}
\alias{load_sumstats_cache}
\alias{logpvec}
\alias{zvec}
\alias{nvec}
\alias{is_present}
\alias{beta_vec}
\alias{se_vec}
\alias{eaf_vec}
\alias{info_vec}
\title{Sumstats objects}
\description{
Load, create, inspect, and cache reference-aligned GWAS summary statistics.
}
\usage{
load_sumstats(path, reference, allow_indel_swaps = FALSE)
create_sumstats(reference, p, z = NULL, n = NULL,
  beta = NULL, se = NULL, eaf = NULL, info = NULL)
save_sumstats_cache(sumstats, path)
load_sumstats_cache(path, shards = NULL)

logpvec(x, ...)
zvec(x, ...)
nvec(x, ...)
is_present(x, ...)
beta_vec(x, ...)
se_vec(x, ...)
eaf_vec(x, ...)
info_vec(x, ...)
}
\arguments{
\item{path}{Character scalar path to a summary-statistics TSV or RDS cache.}
\item{reference}{A \code{ReferencePanel}.}
\item{allow_indel_swaps}{Logical scalar enabling indel swaps. Default is
\code{FALSE}; see \link{allele_matching}.}
\item{sumstats, x}{A \code{Sumstats} object.}
\item{p, z, n, beta, se, eaf, info}{Aligned
numeric vectors. Optional vectors may be \code{NULL}. Missing \code{p}
(\code{NA} or \code{NaN}) sets \code{logpvec} and present optional vectors
to \code{NaN} at those positions. Other \code{p} values must be finite in
\code{[0, 1]}; zero becomes \code{Inf} in \code{logpvec}.}
\item{shards}{Optional character vector of shard labels in canonical order.}
\item{...}{Reserved for S3 methods.}
}
\value{
\code{load_sumstats()}, \code{create_sumstats()}, and
\code{load_sumstats_cache()} return \code{Sumstats} S3 objects. Accessors
return ordinary R vectors, with \code{NULL} for absent optional fields.
}
\details{
Loading prefers exact matches, swaps SNVs automatically, and uses
\code{allow_indel_swaps} for indels. Effects and EAF describe the reference
panel's \code{a1} (ALT) allele. \code{create_sumstats()} takes vectors already
in that order and orientation; it performs no allele matching or swapping.
The TSV \code{p} column is required; empty fields, \code{NA}, and \code{NaN}
are missing values. Rows with missing \code{p} are validated, then dropped
before matching. At \code{info} verbosity, loading reports the number dropped.
Missing reference variants are represented by
\code{NaN} in \code{logpvec(x)} and by \code{is_present(x) == FALSE}.
Gzipped TSV input requires the suggested \pkg{R.utils} package so
\code{data.table::fread()} can read compressed files portably.
}
\seealso{\link{allele_matching}, \code{load_reference}, \code{load_genotype}, \code{set_verbosity}}
\examples{
ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen")
sum_path <- system.file("extdata", "traits_complete.tsv.gz", package = "statgen")
ref <- load_reference(ref_path)
s <- load_sumstats(sum_path, ref)
head(logpvec(s))

cache <- tempfile(fileext = ".rds")
save_sumstats_cache(s, cache)
s_cached <- load_sumstats_cache(cache)
head(logpvec(s_cached))
}
