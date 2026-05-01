# GenotypeShard and GenotypePanel

## Represents

Genotype data consists of per-sample diploid allele counts for a set of
variants. Counts are oriented to `a1` under the allele contract in
[conventions.md](conventions.md). PLINK 1 bfile format is the canonical input.
Other formats (PGEN, BGEN, VCF, dosage) are out of scope and should be
converted by project-specific preprocessing before entering `statgen`.

## Disk representation

Canonical disk format is a PLINK bfile triplet:

```text
PREFIX.bed
PREFIX.bim
PREFIX.fam
```

A sharded panel stores one triplet per chromosome:

```text
genotypes/chrN.bed
genotypes/chrN.bim
genotypes/chrN.fam
```

A non-sharded single PLINK triplet combined across chromosomes is also an
acceptable disk representation.

Converted dense genotype matrices are caches only.

## In-memory objects

A `GenotypeShard` is defined by:

- `chr`: chromosome, including X;
- `bed_path`, `bim_path`, `fam_path`;
- `num_snp` and `num_sample`;
- BIM rows with `chr`, `snp`, `cm`, `bp`, `a1`, `a2`;
- FAM sample IDs in file order;
- allele contract: `.bed` genotypes are interpreted as counts/dosages of `a1`.

A `GenotypePanel` is an ordered collection of `GenotypeShard` objects, usually
one shard per chromosome.

## Representation

Implementations may use language-native containers. The object retains paths
and BIM/FAM metadata. Dense genotype matrix access is deferred.

## API

Genotype implementation is deferred for the current implementation phase. The
contract below records the intended object boundary, but exact matrix access,
missing-genotype handling, sample/variant slicing, and reference-alignment
validation remain to be specified before implementation.

```text
load_genotype(path_or_paths) -> GenotypePanel

GenotypePanel.num_snp -> int
GenotypePanel.num_sample -> int
GenotypePanel.bim -> BIM metadata table
GenotypePanel.fam -> FAM/sample metadata table
```

Expected behavior:

- `path_or_paths` accepts a single bfile prefix (non-sharded) or a list of
  prefixes (sharded).
- loads BIM and FAM metadata; `.bed` genotype data is read on demand.
- no reference argument; no caching.
- Metadata accessors are read-only and may be implemented before genotype
  matrix access. Dense genotype matrix access remains deferred.
