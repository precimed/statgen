\name{reference}
\alias{load_reference}
\alias{save_reference_cache}
\alias{load_reference_cache}
\alias{num_snp}
\alias{shards}
\alias{chr}
\alias{snp}
\alias{bp}
\alias{a1}
\alias{a2}
\alias{a1_hash64}
\alias{a2_hash64}
\alias{is_single_nucleotide_variant}
\alias{is_strand_ambiguous}
\alias{shard_offsets}
\alias{select_shards}
\alias{is_object_compatible}
\alias{validate_checksums}
\alias{save_cache}
\title{ReferencePanel objects}
\description{
Load PLINK BIM reference panels as statgen \code{ReferencePanel} S3 objects,
inspect aligned reference vectors, and save or load R-native reference caches.
}
\usage{
# Constructors and cache I/O
load_reference(path, shards = NULL)
load_reference_cache(path, shards = NULL)
save_reference_cache(panel, path)

# Accessors
num_snp(x, ...)
shards(x, ...)
chr(x, ...)
snp(x, ...)
bp(x, ...)
a1(x, ...)
a2(x, ...)
a1_hash64(x, ...)
a2_hash64(x, ...)
is_single_nucleotide_variant(x, ...)
is_strand_ambiguous(x, ...)
shard_offsets(x, ...)

# Operations
select_shards(x, shards, ...)
is_object_compatible(reference, object, ...)
validate_checksums(x, ...)
save_cache(x, path, ...)
}
\arguments{
\item{path}{Character scalar path. For \code{load_reference()}, this is a BIM
file or sharded BIM template containing \code{@}. For
\code{load_reference_cache()}, this is an RDS reference cache. For cache
writers, this is the output cache path.}
\item{shards}{Optional character vector of shard labels in canonical order.}
\item{panel, reference, x, object}{Reference panel or compatible statgen object.}
\item{...}{Reserved for S3 methods.}
}
\details{
Reference panels are represented as \code{ReferencePanel} S3 objects backed by
an ordered list of reference shards. Reference BIM \code{a1} is ALT and
\code{a2} is REF; other loading functions use this panel's allele orientation.

Most exported functions in this API are S3 generics operating on
\code{ReferencePanel} objects and other reference-aligned statgen objects.
Users should access panel data through accessor functions such as \code{bp(x)}
and \code{a1(x)} rather than directly accessing internal fields.

\code{is_object_compatible(reference, object)} dispatches on the first argument,
which should be a \code{ReferencePanel}. Other reference-aligned panel classes
participate by providing \code{shards()} methods that return their ordered shard
objects.

\code{shards(x)} returns the ordered shard objects. \code{save_cache()} is the
cross-object cache-writing generic; \code{save_reference_cache()} is the
explicit reference-panel cache writer.
}
\value{
\code{load_reference()} and \code{load_reference_cache()} return a
\code{ReferencePanel} S3 object. Accessors return ordinary R vectors or data
frames; allele hashes are \code{bit64::integer64} vectors.
}
\seealso{\link{allele_matching}}
\examples{
path <- system.file("extdata", "reference_chr1.bim", package = "statgen")
ref <- load_reference(path)
num_snp(ref)
head(snp(ref))
sum(is_strand_ambiguous(ref))

cache <- tempfile(fileext = ".rds")
save_reference_cache(ref, cache)
ref_cached <- load_reference_cache(cache)
num_snp(ref_cached)

save_cache(ref, cache)
}
