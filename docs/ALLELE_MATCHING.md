# Allele Matching

See the [matching specification](../spec/reference.md#source-to-reference-matching)
for the formal contract.

## Preparing your data

For genomic-position matching, prepare files with the same genome build,
chromosome names, and DNA strand. ID matching supports different genome builds
and strands for eligible SNVs when IDs identify the same biological variants;
loading does not verify those ID assignments. Genotype matching still requires
the same chromosome labels. Alleles must be uppercase DNA strings. Loading
does not normalize chromosome labels or alleles, perform liftover, or check a
genome FASTA file.

## Which allele do results describe?

Your reference panel is a table of variants that determines the order of
results and which allele they describe. In its BIM file, column 5 (`a1`) is the
alternative allele (ALT), and column 6 (`a2`) is the genome reference allele
(REF). In genotype files, `a1` is the counted allele and `a2` is the other
allele. In summary statistics, `a1` (or `EffectAllele`) is the effect allele,
EAF gives its frequency, and `a2` (or `OtherAllele`) is the other allele.

Loaded effects, EAF, and genotype counts describe the reference panel's `a1`
(ALT) allele. `create_sumstats` expects vectors already in the reference panel's
order and allele orientation; it performs no matching or swaps.
Missing `p` values in these vectors mark absent variants: `logpvec` and all
present optional vectors are `NaN` at those positions.

## Matching variants

`load_sumstats` and `load_genotype` accept `match_by="genomic_position"`
(default) or `match_by="id"`. This selects the position or ID component of
the matching key; both modes also require allele matching. Results preserve
the supplied reference panel's data and row order.

### Genomic-position matching

The default mode identifies variants by chromosome, position, and alleles;
marker IDs do not affect this mode. The loaders prefer a row whose allele
order already matches the reference panel, even when its optional statistics
or calls are missing. Otherwise, they can use a row with reversed alleles.
Sumstats rows with missing `p` are validated and dropped before matching, so
a valid swapped row can match when an exact row has missing `p`.
Each remaining row in your file is used at most once:

- SNV swaps are automatic, including A/T and C/G under the same-strand assumption.
- `allow_indel_swaps` defaults to false. Enabling it permits swaps with exactly
  one single-base allele and one longer allele. Users must establish that
  these reversed indels represent the same biological variant.
- Two multi-base alleles match exactly but never through swap fallback.

Enable indel swaps on each source load when appropriate:

```python
trait = load_sumstats(path, reference, allow_indel_swaps=True)
```

```matlab
trait = statgen.load_sumstats(path, reference, 'allow_indel_swaps', true);
```

```r
trait <- load_sumstats(path, reference, allow_indel_swaps = TRUE)
```

The same option applies to `load_genotype`.

`load_reference` rejects repeated variants and SNPs listed with both allele
orders, such as A/G and G/A at the same position. Indels with reversed allele
orders can still be listed as separate entries.

When multiple rows in your data compete to match one reference variant,
loading fails. A row with matching allele order takes priority over rows with
reversed alleles. Duplicate summary-statistics rows that cannot match any
reference variant are ignored; genotype loading may still reject duplicate
BIM entries.

### ID matching

Use ID mode when marker IDs identify the same biological variants across your
inputs, even when their positions or strands differ:

```python
trait = load_sumstats(path, reference, match_by="id")
genotype = load_genotype(bfile_prefix, reference, match_by="id")
```

```matlab
trait = statgen.load_sumstats(path, reference, 'match_by', 'id');
genotype = statgen.load_genotype(bfile_prefix, reference, 'match_by', 'id');
```

```r
trait <- load_sumstats(path, reference, match_by = "id")
genotype <- load_genotype(bfile_prefix, reference, match_by = "id")
```

Only pairs of distinct single-base A/C/G/T alleles are eligible, excluding
palindromic A/T and C/G pairs. Allele exchange, strand complementation, and both
are accepted equally. For example, A/G, G/A, T/C, and C/T belong to one allele
equivalence group; A/C and its equivalents belong to another. The same ID can
match a variant in each group. Multi-base variants are excluded even if their
alleles match exactly. `allow_indel_swaps=true` is incompatible with ID mode.

Eligible IDs have 1–32 printable, non-whitespace ASCII characters, excluding
the case-sensitive placeholders `.`, `NA`, and `NaN`. Empty, missing,
overlength, or invalid-character IDs are also excluded. `0`, `na`, and `nan`
are ordinary IDs. Matching is exact and case-sensitive: `rs123` and `RS123`
are different IDs and do not match. There is no coordinate fallback.
Eligibility applies to values returned by the existing readers; their quoting
and whitespace handling remain unchanged, so unusual IDs can be read
differently across runtimes.

The loaders search different source scopes:

| Loader | ID-mode source scope |
| --- | --- |
| `load_sumstats` | All source rows across the supplied reference panel; coordinates are ignored entirely. `SNP`, the allele columns, and `P` are required. |
| `load_genotype` | Only rows on the corresponding reference chromosome; source-file routing, samples, and source ploidy are preserved. |

For example, if a reference slot on chromosome 1 has ID `rs123` and alleles
A/G, source rows `rs123:A:G` on chromosome 1 and `rs123:T:C` on chromosome 2
compete for that sumstats slot and cause an ambiguity error. Only the
chromosome 1 row competes when loading genotypes.

Reference ambiguity is assessed across the entire supplied panel in both
loaders. Repeated ID/allele-group keys exclude every reference slot in the
group and produce a warning. Selecting a reference subset before loading
defines the assessment scope. Multiple eligible source rows for a remaining
slot are an error, even when values agree or genotype calls are entirely
missing. There is no exact-match priority in ID mode. Source duplicates that
cannot fill any remaining slot need no ambiguity check. Sumstats rows with
missing `p` are validated and dropped before eligibility and candidate counting.
Unmatched and excluded reference slots remain missing.

## Adjusted values and swap information

When a match requires allele-count/effect inversion, Z and beta change sign
and EAF becomes `1 - EAF`. The results therefore describe the reference
panel's `a1` (ALT) allele. Other statistics and missing values are preserved. Genotype hardcalls
flip as `2 - x` before any ploidy scaling. Missing calls and calls for absent
subjects stay missing.
`genotype.is_allele_swapped` in Python/MATLAB, or
`is_allele_swapped(genotype)` in R, is a full reference-axis logical mask,
true when counts require inversion and false for absent variants. Strand
complementation alone leaves values unchanged and the flag false; relative to
reference A/G, source T/C needs no inversion, while G/A or C/T does. An inverted
variant can still have missing genotype calls. Sumstats expose aligned values
without a swap mask.

## Loading messages

In genomic-position mode at the default `info` verbosity, loading reports
nonzero counts of swapped variants, with an A/T or C/G subset, by input file
and chromosome shard. It also reports unused indel rows that could match an
otherwise unmatched reference variant if indel swaps were enabled. Reversed rows superseded by exact matches
are not reported.

In ID mode, loading reports nonzero source and reference counts of excluded
palindromic and multi-base variants, plus matched rows requiring inversion or
strand complementation. Rows requiring both contribute to both match counts.
Unusable IDs and ambiguous reference keys produce warnings.

Use `set_verbosity("quiet")` to hide progress messages and loading summaries;
warnings remain visible. MATLAB uses `statgen.set_verbosity('quiet')`.

## Long alleles and historical reference caches

Source loaders validate and preserve complete uppercase DNA sequences. Matching
hashes use the first 1,000,000 bases of each allele. Longer alleles produce one
warning per input file; alleles sharing this prefix cannot be distinguished by
matching or hash-key duplicate checks, even if their suffixes differ.

Older reference caches, including those inside older LD distributions, may use
only the first 150 bases. They may fail to match alleles longer than 150 bases
in newly loaded sources, even though the reference checksum is unchanged.
Cache loading retains stored hashes and mappings. Rebuild the reference or its
cache from canonical sources to use the new hashes. For an LD distribution,
keep its reference row order and allele orientation when rebuilding its
reference cache.

## Genotype caches and BED overrides

New genotype caches use `genotype_cache/0.2` and retain swap flags. Existing
`/0.1` caches remain readable with all-false flags, preserving their stored
mappings and calls. Reload the original files to obtain additional matches
from the new matching policy; loading an old cache does not rematch variants.
Saving that loaded cache writes `/0.2` while retaining its mappings.

`allow_indel_swaps` and `match_by` apply when loading source data and are not
stored in caches. Cache loads and shard subsets retain the stored values,
mappings, and orientation without rematching. Reload sources to change the
matching mode. Cache formats and reference checksums are unchanged.

Genotype caches contain metadata, so BED files are still needed to fetch
calls. If you move a BED file, supply its new location with `bed_path`.
A replacement must preserve the original variant and sample layout and allele
order, including any chromosome-specific sample subsets. Stored swap flags
still apply to the replacement BED data.
