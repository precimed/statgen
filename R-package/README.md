# statgen R Package

Statistical genetics data objects and loaders for R. This package is part of the larger `statgen` repository, which also includes Python and MATLAB/Octave implementations.

## Installation

### From GitHub (development)

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install from the R-package directory
devtools::install()

# Or from the repository root
devtools::install("R-package")
```

### From CRAN (release)

```r
install.packages("statgen")
```

## Documentation

- [Main statgen repository](https://github.com/precimed/statgen)
- Package vignette: `vignette("statgen")`

## Allele Alignment

Loaders match reversed SNVs automatically; indel swaps are opt-in.
See `?allele_matching` and `vignette("statgen")` for details and examples.

## Annotation Loading

`load_annotations()` loads one or more 3-column BED files as binary annotation
columns. `load_annotation()` loads one BED-like source file and can paint
selected numeric value columns as continuous annotations or expand one
`group_column` into multiple binary annotations. Both return an
`AnnotationPanel` whose `annomat()` accessor is a sparse numeric
`Matrix::dgCMatrix`; `is_binary()` and `annotation_metadata()` expose per-column
metadata.
