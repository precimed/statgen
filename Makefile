R_PACKAGE_VERSION := $(shell awk '/^Version:/ {print $$2}' R-package/DESCRIPTION)
R_PACKAGE_TARBALL := /tmp/statgen_$(R_PACKAGE_VERSION).tar.gz
R_DOCS_DIR := /tmp/statgen-docs
R_DOCS_LIB := /tmp/statgen-r-lib

.PHONY: install fixtures prepare-r-fixtures test test-python test-octave test-matlab test-r r-manual r-vignettes r-docs rcmd-check

install:
	pip install -e python/

fixtures:
	python tests/fixtures/generate.py
	$(MAKE) prepare-r-fixtures

prepare-r-fixtures:
	mkdir -p R-package/inst/extdata
	cp tests/fixtures/reference/sharded/1.bim R-package/inst/extdata/reference_chr1.bim
	cp tests/fixtures/reference/sharded/X.bim R-package/inst/extdata/reference_chrX.bim
	cp tests/fixtures/sumstats/traits.tsv.gz R-package/inst/extdata/traits.tsv.gz
	cp tests/fixtures/sumstats/traits_complete.tsv.gz R-package/inst/extdata/traits_complete.tsv.gz
	cp tests/fixtures/annotations/anno1.bed R-package/inst/extdata/anno1.bed
	cp tests/fixtures/annotations/anno2.bed R-package/inst/extdata/anno2.bed
	cp tests/fixtures/annotations/grouped.annot R-package/inst/extdata/grouped.annot
	cp tests/fixtures/genotype/sharded/1.bed R-package/inst/extdata/genotype_1.bed
	cp tests/fixtures/genotype/sharded/1.bim R-package/inst/extdata/genotype_1.bim
	cp tests/fixtures/genotype/sharded/1.fam R-package/inst/extdata/genotype_1.fam
	cp tests/fixtures/genotype/sharded/X.bed R-package/inst/extdata/genotype_X.bed
	cp tests/fixtures/genotype/sharded/X.bim R-package/inst/extdata/genotype_X.bim
	cp tests/fixtures/genotype/sharded/X.fam R-package/inst/extdata/genotype_X.fam
	cp tests/fixtures/genotype/sharded/X.ploidy R-package/inst/extdata/genotype_X.ploidy
	mkdir -p R-package/inst/extdata/allele_swaps
	cp -r tests/fixtures/allele_swaps/. R-package/inst/extdata/allele_swaps/
	$(RM) -r R-package/inst/extdata/ld
	mkdir -p R-package/inst/extdata/ld
	cp tests/fixtures/ld/python/ld_manifest.json R-package/inst/extdata/ld/ld_manifest.json
	cp tests/fixtures/ld/python/reference_cache.npz R-package/inst/extdata/ld/reference_cache.npz
	cp tests/fixtures/ld/python/reference_chr1.bim R-package/inst/extdata/ld/reference_chr1.bim
	cp tests/fixtures/ld/python/reference_chrX.bim R-package/inst/extdata/ld/reference_chrX.bim
	cp tests/fixtures/ld/python/ld_chr1.npz R-package/inst/extdata/ld/ld_chr1.npz
	cp tests/fixtures/ld/python/ld_chrX_female.npz R-package/inst/extdata/ld/ld_chrX_female.npz
	cp tests/fixtures/ld/python/ld_chrX_male.npz R-package/inst/extdata/ld/ld_chrX_male.npz
	cp tests/fixtures/ld/python/ld_chrX_combined.npz R-package/inst/extdata/ld/ld_chrX_combined.npz
	Rscript -e "local({ for (f in list.files('R-package/R', pattern = '[.]R$$', full.names = TRUE)) source(f); prepare_ld_npz_for_r('R-package/inst/extdata/ld') })"

test: prepare-r-fixtures
	pytest tests/

test-python: prepare-r-fixtures
	pytest tests/ -m "not octave"

test-octave:
	pytest tests/ -m octave

test-matlab:
	STATGEN_MATLAB=1 pytest tests/ -m octave

test-r: prepare-r-fixtures
	pytest tests/ -m r

r-manual:
	$(RM) /tmp/statgen-manual.pdf
	R CMD Rd2pdf R-package --output=/tmp/statgen-manual.pdf

r-vignettes: prepare-r-fixtures
	mkdir -p $(R_DOCS_DIR) $(R_DOCS_LIB)
	R CMD INSTALL --library=$(R_DOCS_LIB) R-package
	Rscript -e ".libPaths(c('$(R_DOCS_LIB)', .libPaths())); rmarkdown::render('R-package/vignettes/statgen.Rmd', output_dir = '$(R_DOCS_DIR)')"

r-docs: r-vignettes
	$(RM) $(R_DOCS_DIR)/statgen-manual.pdf
	R CMD Rd2pdf R-package --output=$(R_DOCS_DIR)/statgen-manual.pdf

rcmd-check: prepare-r-fixtures
	cd /tmp && R CMD build $(CURDIR)/R-package
	R CMD check --as-cran --output=/tmp $(R_PACKAGE_TARBALL)
