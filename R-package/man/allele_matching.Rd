\name{allele_matching}
\alias{allele_matching}
\title{Allele Matching and Orientation}
\description{
How loading makes summary statistics and genotype counts consistent with your
reference panel.
}
\section{Preparing your data}{
Prepare all files so they use the same genome build, chromosome names, and
DNA strand. \pkg{statgen} does not repair strand differences, normalize alleles,
or check alleles against a genome FASTA file. Alleles must be uppercase DNA
strings.
}
\section{Which allele do results describe?}{
Your reference panel is a table of variants that determines the order of
results and which allele they describe. The allele columns have different
meanings depending on the file:
\itemize{
\item In the reference panel's BIM file, \code{a1} (column 5) is the alternative
allele (ALT), and \code{a2} (column 6) is the genome reference allele (REF).
\item In a genotype BIM file, \code{a1} is the allele counted by genotype calls,
and \code{a2} is the other allele.
\item In summary statistics, \code{a1} (or \code{EffectAllele}) is the allele
described by the effect size, Z score, and effect allele frequency (EAF).
\code{a2} (or \code{OtherAllele}) is the other allele.
}
Loaded effects, EAF, and genotype counts describe the reference panel's
\code{a1} (ALT) allele. \code{create_sumstats()} expects vectors already in the
reference panel's order and allele orientation; it does not match or swap
alleles.
}
\section{Matching variants}{
\code{load_sumstats()} and \code{load_genotype()} identify variants by
chromosome, position, and alleles; marker IDs do not affect matching. They
prefer a row whose allele order already matches the reference panel, even
when its optional statistics or calls are missing. Otherwise, they can use
a row with reversed alleles. Sumstats rows with missing \code{p} are validated
and dropped before matching, so a valid swapped row can match when an exact
row has missing \code{p}. Each remaining row is used at most once.

Single nucleotide variants (SNVs) are swapped automatically, including A/T
and C/G pairs, assuming your files already use the same strand. Swapping
insertions and deletions is off by default. Enable it with
\code{allow_indel_swaps = TRUE} on either loading function only when you know
the reversed alleles describe the same biological variant. For example:
\preformatted{
trait <- load_sumstats(path, reference, allow_indel_swaps = TRUE)
genotype <- load_genotype(bfile_prefix, reference, allow_indel_swaps = TRUE)
}
This option applies to pairs with exactly one single-base and one longer
allele. Pairs where both alleles have multiple bases must already match the
reference panel's allele order.

\code{load_reference()} rejects repeated variants and SNPs listed with both
allele orders, such as A/G and G/A at the same position. Indels with reversed
allele orders can be listed as separate entries.

When multiple rows in your data compete to match one reference variant,
loading fails. A row with matching allele order takes priority over rows with
reversed alleles. Duplicate summary-statistics rows that cannot match any
reference variant are ignored; genotype loading may still reject duplicate
BIM entries.
}
\section{Adjusted values and swap information}{
When a match requires a swap, Z and beta change sign and EAF becomes
\code{1 - EAF}. Other statistics (p, n, se, and info) and missing values are
preserved. Genotype hardcalls are flipped as \code{2 - x} before any ploidy
scaling. Missing calls and calls for absent subjects remain missing.

\code{is_allele_swapped(genotype)} returns a logical vector in reference panel
order. It is true for variants matched by swapping alleles and false for
exact matches and absent variants. A swapped variant can still have missing
genotype calls. Summary-statistics objects expose adjusted values without a
swap mask.
}
\section{Loading messages}{
At the default \code{"info"} verbosity, loading reports nonzero counts of
swapped variants, with an A/T or C/G subset, by input file and chromosome
shard. It also reports unused indel rows that could match an otherwise
unmatched reference variant if indel swaps were enabled. Reversed rows
superseded by exact matches are not reported.

Use \code{set_verbosity("quiet")} to hide progress messages and loading
summaries. Warnings remain visible.
}
\section{Genotype caches and BED overrides}{
New genotype caches use \code{genotype_cache/0.2} and preserve swap information.
Historical \code{/0.1} caches remain readable with their original mappings
and calls, treating all stored variants as unswapped. Loading an old cache
does not find additional matches; reload the original files to obtain them.
Saving a loaded old cache writes \code{/0.2} while retaining its mappings.
Summary-statistics caches keep their existing schema and adjusted values.

The indel option applies when loading your data and is not stored in caches.
Cache loads and chromosome subsets retain the stored allele orientation. To
change the option, load the original files again and save a new cache.

Genotype caches contain metadata, so BED files are still needed when fetching
calls. If you move a BED file, supply its new location with \code{bed_path}.
A replacement must preserve the original variant and sample layout and
allele order, including any chromosome-specific sample subsets. Stored swap
information still applies to the replacement BED data.
}
\seealso{
\code{\link{load_reference}}, \code{\link{load_sumstats}},
\code{\link{create_sumstats}}, \code{\link{load_genotype}},
\code{\link{is_allele_swapped}}, \code{\link{fetch_genotypes}},
\code{\link{load_genotype_cache}}, \code{\link{set_verbosity}}
}
