\name{allele_matching}
\alias{allele_matching}
\title{Allele Matching and Orientation}
\description{
How loading makes summary statistics and genotype counts consistent with your
reference panel.
}
\section{Preparing your data}{
For genomic-position matching, prepare files with the same genome build,
chromosome names, and DNA strand. ID matching allows different builds and
strands for eligible SNVs when IDs identify the same biological variants.
\pkg{statgen} does not verify those ID assignments, normalize chromosome names
or alleles, perform liftover, or check a genome FASTA. Alleles must be uppercase
DNA strings.
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
\section{Genomic-position matching}{
\code{load_sumstats()} and \code{load_genotype()} default to
\code{match_by = "genomic_position"}, using chromosome, position, and alleles;
marker IDs do not affect this mode. They
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
\section{ID matching}{
Set \code{match_by = "id"} to match marker IDs and alleles instead of positions:
\preformatted{
trait <- load_sumstats(path, reference, match_by = "id")
genotype <- load_genotype(bfile_prefix, reference, match_by = "id")
}
Both matching modes require allele matching. ID mode accepts distinct
single-base A/C/G/T alleles, excluding palindromic A/T and C/G pairs. It accepts
allele exchange, strand complementation, or both: A/G, G/A, T/C, and C/T form
one equivalence group. A/C and its equivalents form a different group, so one
ID may identify a variant in each group. Multi-base alleles are excluded even
when their order matches exactly. \code{allow_indel_swaps = TRUE} cannot be
combined with ID matching.

IDs must contain 1--32 printable non-whitespace ASCII characters. Empty,
missing, overlength, or invalid-character IDs and the placeholders \code{.},
\code{NA}, and \code{NaN} are excluded; \code{0} is an ordinary ID. Comparison is exact and case-sensitive:
\code{rs123} and \code{RS123} do not match. There is no coordinate fallback.
Eligibility applies to IDs returned by the existing readers, including their
quoting and whitespace handling. Original reference strings and row order
are preserved.

Summary statistics require \code{SNP}, the allele columns, and \code{P};
coordinate columns are optional and ignored. Rows with missing P are validated,
then dropped before matching. All remaining source rows are searched across
the supplied reference panel. Genotype loading retains chromosome-based file,
sample, and ploidy handling, searching only the corresponding chromosome.
Thus source rows with the same ID/allele key on two chromosomes are ambiguous
for sumstats, but only the same-chromosome row competes for a genotype slot.

Repeated reference ID/allele keys are assessed across the entire supplied panel;
all slots in those groups are excluded with a warning. Selecting a reference
subset before loading defines that assessment scope. Multiple eligible source
rows for a remaining slot are an error, even if values agree or genotype calls
are all missing. There is no exact-match priority in ID mode. Unmatched or
excluded slots remain missing.
}
\section{Adjusted values and swap information}{
When a match requires allele-count/effect inversion, Z and beta change sign and EAF becomes
\code{1 - EAF}. Other statistics (p, n, se, and info) and missing values are
preserved. Genotype hardcalls are flipped as \code{2 - x} before any ploidy
scaling. Missing calls and calls for absent subjects remain missing.

\code{is_allele_swapped(genotype)} returns a logical vector in reference panel
order. It is true when counts are inverted and false for absent variants.
Strand complementation alone leaves values unchanged and the flag false;
complementation plus allele exchange requires inversion. A swapped variant can still have missing
genotype calls. Summary-statistics objects expose adjusted values without a
swap mask.
}
\section{Loading messages}{
In genomic-position mode at the default \code{"info"} verbosity, loading reports
nonzero counts of swapped variants, with an A/T or C/G subset, by input file and chromosome
shard. It also reports unused indel rows that could match an otherwise
unmatched reference variant if indel swaps were enabled. Reversed rows
superseded by exact matches are not reported.

In ID mode, loading reports nonzero counts of excluded palindromic and
multi-base variants, and matched rows requiring inversion or complementation.
Unusable IDs and ambiguous reference keys produce warnings.

Use \code{set_verbosity("quiet")} to hide progress messages and loading
summaries. Warnings remain visible.
}
\section{Long alleles and historical reference caches}{
Source loaders validate and preserve complete uppercase DNA sequences. Matching
hashes use the first 1,000,000 bases of each allele. Longer alleles produce one
warning per input file. Matching and hash-key duplicate checks cannot distinguish
alleles that share this prefix but differ afterwards.

Older reference caches, including those inside older LD distributions, may use
only the first 150 bases. They may fail to match alleles longer than 150 bases
in newly loaded sources, despite unchanged reference checksums. Cache loading
retains stored hashes and mappings. Rebuild the reference or its cache from
canonical sources to use the new hashes. When rebuilding an LD distribution's
reference cache, preserve reference row order and allele orientation.
}
\section{Genotype caches and BED overrides}{
New genotype caches use \code{genotype_cache/0.2} and preserve swap information.
Historical \code{/0.1} caches remain readable with their original mappings
and calls, treating all stored variants as unswapped. Loading an old cache
does not find additional matches; reload the original files to obtain them.
Saving a loaded old cache writes \code{/0.2} while retaining its mappings.
Summary-statistics caches keep their existing schema and adjusted values.

The \code{allow_indel_swaps} and \code{match_by} options apply when loading
source data and are not stored in caches. Cache loads and chromosome subsets
retain the stored values, mappings, and allele orientation. To change matching
options, load the original files again and save a new cache.

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
