\name{genotype}
\alias{load_genotype}
\alias{save_genotype_cache}
\alias{load_genotype_cache}
\alias{num_sample}
\alias{fid}
\alias{iid}
\alias{father_id}
\alias{mother_id}
\alias{sex}
\alias{is_male}
\alias{is_female}
\alias{ploidy_male}
\alias{ploidy_female}
\alias{source_layout}
\alias{source_row0}
\alias{is_allele_swapped}
\alias{is_subject_present}
\alias{fetch_genotypes_int8}
\alias{fetch_genotypes}
\title{GenotypePanel objects}
\description{
Load PLINK bfile-backed genotype metadata, inspect reference-aligned genotype
panels, cache metadata, and fetch hardcall genotype slices from BED files.
}
\usage{
load_genotype(bfile_prefix, reference, allow_indel_swaps = FALSE)
save_genotype_cache(panel, path)
load_genotype_cache(path, shards = NULL)

num_sample(x, ...)
fid(x, ...)
iid(x, ...)
father_id(x, ...)
mother_id(x, ...)
sex(x, ...)
is_male(x, ...)
is_female(x, ...)
ploidy_male(x, ...)
ploidy_female(x, ...)
source_layout(x, ...)
source_row0(x, ...)
is_allele_swapped(x, ...)
is_subject_present(x, shard, ...)

fetch_genotypes_int8(x, snp_indices, bed_path = NULL, ...)
fetch_genotypes(x, snp_indices, bed_path = NULL,
  haploid_mode = NULL, ...)
}
\arguments{
\item{bfile_prefix}{Character scalar PLINK bfile prefix, without suffixes.
Use \code{@} as a shard-label placeholder for sharded input.}
\item{reference}{A \code{ReferencePanel}.}
\item{allow_indel_swaps}{Logical scalar enabling swaps for allele pairs with
one single-base and one longer allele. SNV swaps are always enabled.}
\item{panel, x}{A \code{GenotypePanel}.}
\item{path}{Character scalar path to an RDS genotype metadata cache.}
\item{shards}{Optional character vector of shard labels in canonical order.}
\item{shard}{A loaded shard label for subject-presence lookup.}
\item{snp_indices}{One-based panel-global SNP indices.}
\item{bed_path}{Optional BED path override. A sharded override may contain
\code{@}.}
\item{haploid_mode}{Optional output mode, \code{"raw"} or
\code{"ploidy_scaled"}. The default is \code{"raw"}.}
\item{...}{Reserved for S3 methods.}
}
\value{
\code{load_genotype()} and \code{load_genotype_cache()} return
\code{GenotypePanel} S3 objects. Metadata accessors return ordinary R vectors.
\code{fetch_genotypes_int8()} returns an integer matrix with samples as rows
and requested SNPs as columns, counting reference \code{a1} and using
\code{-1} for missing hardcalls.
\code{fetch_genotypes()} returns a numeric matrix with missing calls as
\code{NaN}.
}
\details{
Genotype panels are reference-aligned metadata handles for PLINK BED payloads.
Hardcalls are decoded on demand; caches store only metadata needed to locate
and map BED rows. Exact matches take priority over allele swaps.
\code{is_allele_swapped()} returns a logical vector on the full reference axis;
exact matches and absent variants are false. Swapped nonmissing calls are
flipped before ploidy scaling.
}
\examples{
ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen")
geno_bed <- system.file("extdata", "genotype_1.bed", package = "statgen")
geno_prefix <- sub("\\\\.bed$", "", geno_bed)
ref <- load_reference(ref_path)
g <- load_genotype(geno_prefix, ref)
geno <- fetch_genotypes(g, c(1, 2))
dim(geno)
head(geno)

cache <- tempfile(fileext = ".rds")
save_genotype_cache(g, cache)
g_cached <- load_genotype_cache(cache)
num_sample(g_cached)
sum(is_present(g_cached))
}
