# Allele Matching

See the [matching specification](../spec/reference.md#source-to-reference-matching)
for the formal contract.

## Preparing your data

Prepare all files with the same genome build, chromosome names, and DNA
strand before loading them. Alleles must be uppercase DNA strings. Loading
does not repair strand differences, normalize alleles, or check them against
a genome FASTA file.

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

## Matching variants

`load_sumstats` and `load_genotype` identify variants by chromosome, position,
and alleles; marker IDs do not affect matching. They prefer a row whose allele
order already matches the reference panel, even when its optional statistics
or calls are missing. Otherwise, they can use a row with reversed alleles.
Each row in your file is used at most once:

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

## Adjusted values and swap information

When a match requires a swap, Z and beta change sign and EAF becomes
`1 - EAF`. The results therefore describe the reference panel's `a1` (ALT)
allele. Other statistics and missing values are preserved. Genotype hardcalls
flip as `2 - x` before any ploidy scaling. Missing calls and calls for absent
subjects stay missing.
`genotype.is_allele_swapped` in Python/MATLAB, or
`is_allele_swapped(genotype)` in R, is a full reference-axis logical mask,
false for exact matches and absent variants. A swapped variant can still have
missing genotype calls. Sumstats expose aligned values without a swap mask.

## Loading messages

At the default `info` verbosity, loading reports nonzero counts of swapped
variants, with an A/T or C/G subset, by input file and chromosome shard. It also
reports unused indel rows that could match an otherwise unmatched reference
variant if indel swaps were enabled. Reversed rows superseded by exact matches
are not reported. Use
`set_verbosity("quiet")` to hide progress messages and loading summaries;
warnings remain visible. MATLAB uses `statgen.set_verbosity('quiet')`.

## Genotype caches and BED overrides

New genotype caches use `genotype_cache/0.2` and retain swap flags. Existing
`/0.1` caches remain readable with all-false flags, preserving their stored
mappings and calls. Reload the original files to obtain additional matches
from the new matching policy; loading an old cache does not rematch variants.
Saving that loaded cache writes `/0.2` while retaining its mappings.

The indel option belongs to loading your data and is not stored in caches.
Cache loads and shard subsets retain stored orientation. Sumstats caches keep
their existing schema and aligned values.

Genotype caches contain metadata, so BED files are still needed to fetch
calls. If you move a BED file, supply its new location with `bed_path`.
A replacement must preserve the original variant and sample layout and allele
order, including any chromosome-specific sample subsets. Stored swap flags
still apply to the replacement BED data.
