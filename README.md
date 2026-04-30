# statgen

`statgen` is the working design area for a paired Python and MATLAB/Octave
package for statistical genetics data objects and basic inference utilities.

The canonical specification starts at [spec/SPEC.md](spec/SPEC.md). The spec is
object-centered, with separate files for:

- reference shards and panels;
- genotype shards and panels;
- LD shards and panels;
- annotation shards and panels;
- GWAS summary statistics;
- implementation and testing strategy.

`statgen` consumes post-harmonization inputs. Variant build, `genomatch` NCBI
contig naming, and allele orientation are expected to be resolved upstream by
`genomatch` or an equivalent pipeline. Loaders preserve contig labels and
allele fields exactly; they do not normalize aliases such as `chr1` to `1` or
swap alleles.

## Object overview

| Object | Meaning | Portable disk representation | MATLAB/Octave representation |
| --- | --- | --- | --- |
| `ReferenceShard` / `ReferencePanel` | One chromosome SNP table / ordered multi-chromosome SNP table | External input: sharded (one `.bim` per chromosome, `@` placeholder) or non-sharded (single `.bim`); auto-split by `chr` column on load | Struct with SNP columns; panel struct with ordered shards and optional offsets |
| `GenotypeShard` / `GenotypePanel` | One chromosome genotype shard / ordered genotype panel | External input: sharded (one PLINK `.bed/.bim/.fam` per chromosome) or single non-sharded bfile | Struct with paths and BIM/FAM metadata; `.bed` read on demand |
| `LDShard` / `LDPanel` | One chromosome sparse LD / ordered LD panel | Built by `statgen_build_ld.py`: sharded or non-sharded to match the reference; upper-triangle directory with `metadata.json`, `ld_idx1.i32`, `ld_idx2.i32`, `ld_r.f32`, `mafvec.f32`; tied to a specific reference (checksum in `metadata.json`) | Sparse matrix `LD_r` of size `num_snp x num_snp`; tied to the paired reference; optional `.mat` cache |
| `AnnotationShard` / `AnnotationPanel` | SNP annotation matrices painted from BED intervals onto a reference | External BED files (one per annotation); painted `annomat` is in-memory only; tied to the paired reference | Dense matrix `annomat` of size `num_snp x num_annot`; optional `.mat` cache |
| `SumstatsShard` / `Sumstats` | One sumstats shard / ordered shards for one trait/source | External input: single `.tsv.gz` per trait; rows matched to reference by `chr:bp:a1:a2` on load | Shard struct of aligned vectors; `Sumstats` struct with ordered shard cells; tied to the paired reference; optional aligned `.mat` cache |

## API provenance

| Object | `load` | `save` | `save_cache` / `load_cache` |
| --- | --- | --- | --- |
| Reference | ✓ from `.bim` | — external | ✓ |
| Genotype | ✓ from bfile | — external | — |
| LD | ✓ from binary triplets | — build script | ✓ sparse `.mat` |
| Annotations | ✓ paint from BED | — external | ✓ painted `.mat` |
| Sumstats | ✓ from `.tsv.gz` | — external | ✓ aligned `.mat` |

Shared preprocessing scripts may use external tools or Python utility
functions. Scripts should not depend on MATLAB. MATLAB functionality should be
validated with Octave through `pytest`-driven tests.

Planned implementation layout:

```text
python
matlab
script
tests
spec
```
