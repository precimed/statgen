# statgen 0.3.4

* `load_annotation()` gains a `group_column` argument to create one binary annotation per distinct group value in a source column.

# statgen 0.3.3

* New `load_annotation()` function loads a single BED-like annotation file, supporting multiple value columns, optional headers, annotation name overrides, and sidecar metadata files.
* `AnnotationPanel` now exposes `is_binary` and `annotation_metadata` accessors.
* Annotation matrices are now sparse numeric matrices (previously binary-only).

# statgen 0.3.2

* Initial CRAN submission.
