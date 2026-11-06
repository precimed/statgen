"""R phase 2 acceptance tests: sumstats and annotations."""

import gzip
import json
import math

import numpy as np
import pytest

from statgen._utils import allele_hash64
from statgen.annotations import load_annotations
from statgen.reference import load_reference
from statgen.sumstats import load_sumstats
from tests.conftest import FIXTURES_DIR, R_PACKAGE_DIR, run_rscript, skipif_no_rscript

SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
SUMSTATS = FIXTURES_DIR / "sumstats/traits.tsv.gz"
ANNOTATIONS = [FIXTURES_DIR / "annotations/anno1.bed", FIXTURES_DIR / "annotations/anno2.bed"]


def _source_phase2_script(expr: str) -> str:
    files = [
        "utils.R",
        "verbosity.R",
        "hash.R",
        "variant_match.R",
        "bfile_utils.R",
        "reference.R",
        "sumstats.R",
        "annotation_painting.R",
        "annotations.R",
    ]
    sources = " ".join(
        f"source({json.dumps(str(R_PACKAGE_DIR / 'R' / filename))});"
        for filename in files
    )
    return f"{sources} {expr}"


def _parse_float_csv(text: str) -> list[float]:
    out = []
    for token in text.split(","):
        if token == "NaN":
            out.append(math.nan)
        elif token == "Inf":
            out.append(math.inf)
        else:
            out.append(float(token))
    return out


@pytest.mark.r
@skipif_no_rscript
def test_r_allele_hash64_matches_python_known_values():
    alleles = ["A", "C", "G", "T", "N", "-", "ACGT", "ATCGGCTA"]
    expected = ",".join(str(int(x)) for x in allele_hash64(alleles))
    result = run_rscript(
        _source_phase2_script(
            "h <- .allele_hash64(c('A', 'C', 'G', 'T', 'N', '-', 'ACGT', 'ATCGGCTA')); "
            "cat(paste(as.character(h), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == expected


@pytest.mark.r
@skipif_no_rscript
def test_r_allele_hash64_warns_and_truncates_long_alleles():
    expected = str(int(allele_hash64(["A" * 150])[0]))
    result = run_rscript(
        _source_phase2_script(
            "h <- suppressWarnings(.allele_hash64(paste(rep('A', 151), collapse=''))); "
            "cat(as.character(h), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == expected


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_matches_python_fixture_alignment():
    reference = load_reference(SHARDED_REF)
    with pytest.warns(RuntimeWarning):
        py = load_sumstats(SUMSTATS, reference)

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(SUMSTATS))}, ref)); "
            "cat(num_snp(s), '\\n'); "
            "cat(paste(vapply(shards(s), function(x) x$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(format(logpvec(s), digits=17), collapse=','), '\\n'); "
            "cat(paste(format(zvec(s), digits=17), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_present(s)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "8"
    assert lines[1] == "1,X"
    np.testing.assert_allclose(_parse_float_csv(lines[2]), py.logpvec, equal_nan=True)
    np.testing.assert_allclose(_parse_float_csv(lines[3]), py.zvec, equal_nan=True)
    assert lines[4] == ",".join(str(int(x)) for x in py.is_present)


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_ignores_extra_columns(tmp_path):
    path = tmp_path / "extra_cols.tsv.gz"
    with gzip.open(path, "wt") as f:
        f.write(
            "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tN\tStudyN\tDirection\n"
            "1\t100\trs1\tA\tG\t0.01\t2.5\t1000\t999\t+\n"
            "X\t100\trsx1\tA\tG\t0.003\t3.0\t500\t499\t-\n"
        )

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(path))}, ref)); "
            "cat(num_snp(s), '\\n'); "
            "cat(paste(format(zvec(s)[c(1, 6)], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(nvec(s)[c(1, 6)], digits=17), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "8"
    np.testing.assert_allclose(_parse_float_csv(lines[1]), np.array([2.5, 3.0]))
    np.testing.assert_allclose(_parse_float_csv(lines[2]), np.array([1000, 500]))


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_cache_subset_and_create(tmp_path):
    cache = tmp_path / "sumstats.rds"
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(SUMSTATS))}, ref)); "
            f"save_sumstats_cache(s, {json.dumps(str(cache))}); "
            f"loaded <- suppressWarnings(load_sumstats_cache({json.dumps(str(cache))}, shards = 'X')); "
            "created <- create_sumstats(ref, p = rep(0.5, num_snp(ref)), z = rep(1, num_snp(ref)), n = rep(100, num_snp(ref)), beta = rep(0.1, num_snp(ref)), se = rep(0.2, num_snp(ref)), eaf = rep(0.3, num_snp(ref)), info = rep(0.9, num_snp(ref))); "
            "cat(num_snp(loaded), '\\n'); "
            "cat(paste(format(logpvec(loaded), digits=17), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_present(created)), collapse=','), '\\n'); "
            "cat(paste(format(nvec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(beta_vec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(se_vec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(eaf_vec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(info_vec(created)[1:3], digits=17), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "3"
    np.testing.assert_allclose(
        _parse_float_csv(lines[1]),
        np.array([-math.log10(0.003), -math.log10(0.6), math.nan]),
        equal_nan=True,
    )
    assert lines[2] == "1,1,1,1,1,1,1,1"
    np.testing.assert_allclose(_parse_float_csv(lines[3]), np.array([100, 100, 100]))
    np.testing.assert_allclose(_parse_float_csv(lines[4]), np.array([0.1, 0.1, 0.1]))
    np.testing.assert_allclose(_parse_float_csv(lines[5]), np.array([0.2, 0.2, 0.2]))
    np.testing.assert_allclose(_parse_float_csv(lines[6]), np.array([0.3, 0.3, 0.3]))
    np.testing.assert_allclose(_parse_float_csv(lines[7]), np.array([0.9, 0.9, 0.9]))


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_constructor_rejects_inconsistent_optional_fields():
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            "s <- create_sumstats(ref, p = rep(0.5, num_snp(ref)), z = rep(1, num_snp(ref))); "
            "bad <- s; bad$shards[[2]]$zvec <- NULL; "
            "ok <- FALSE; tryCatch(.new_sumstats_panel(bad$shards), error = function(e) ok <<- grepl('consistent zvec presence', e$message)); "
            "ok2 <- FALSE; tryCatch(zvec(bad), error = function(e) ok2 <<- grepl('inconsistent zvec presence', e$message)); "
            "stopifnot(ok, ok2)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_annotations_match_python_fixture_matrix_and_cache(tmp_path):
    reference = load_reference(SHARDED_REF)
    py = load_annotations(ANNOTATIONS, reference)
    expected_col_major = py.annomat.toarray().astype(int).ravel(order="F")
    cache = tmp_path / "annotations.rds"

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"a <- load_annotations(c({json.dumps(str(ANNOTATIONS[0]))}, {json.dumps(str(ANNOTATIONS[1]))}), ref); "
            f"save_annotations_cache(a, {json.dumps(str(cache))}); "
            f"loaded <- load_annotations_cache({json.dumps(str(cache))}); "
            "selected <- select_annotations(a, c('anno2', 'anno1')); "
            "cat(num_snp(a), '\\n'); "
            "cat(num_annot(a), '\\n'); "
            "cat(paste(annonames(a), collapse=','), '\\n'); "
            "cat(class(annomat(a))[[1]], '\\n'); "
            "cat(paste(as.integer(is_binary(a)), collapse=','), '\\n'); "
            "cat(paste(grepl('source_file', annotation_metadata(a)), collapse=','), '\\n'); "
            "cat(paste(as.integer(as.vector(as.matrix(annomat(loaded)))), collapse=','), '\\n'); "
            "cat(paste(annonames(selected), collapse=','), '\\n'); "
            "cat(paste(colnames(annomat(selected)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[:6] == ["8", "2", "anno1,anno2", "dgCMatrix", "1,1", "TRUE,TRUE"]
    # R stores matrices in column-major order; this flattening matches the
    # Python expected_col_major vector built from the dense annotation matrix.
    assert lines[6] == ",".join(str(x) for x in expected_col_major.tolist())
    assert lines[7] == "anno2,anno1"
    assert lines[8] == "anno2,anno1"


@pytest.mark.r
@skipif_no_rscript
def test_r_load_annotation_continuous_headered_and_cache_metadata(tmp_path):
    annot = tmp_path / "wide.annot"
    annot.write_text(
        "chrom\tstart0\tend0\tscore\tweight\n"
        "1\t99\t200\t0.5\t10\n"
        "1\t299\t400\t1.5\t20\n"
        "X\t99\t301\t2.5\t30\n"
    )
    sidecar = tmp_path / "wide.meta"
    sidecar.write_text("chrom meta\nstart meta\nend meta\nscore meta\nweight meta", encoding="utf-8")
    cache = tmp_path / "annotations.rds"

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"a <- load_annotation({json.dumps(str(annot))}, ref, has_header = TRUE, value_columns = c('weight', 'score'), annotation_metadata_path = {json.dumps(str(sidecar))}); "
            f"save_annotations_cache(a, {json.dumps(str(cache))}); "
            f"b <- load_annotations_cache({json.dumps(str(cache))}); "
            "M <- as.matrix(annomat(b)); "
            "cat(paste(annonames(b), collapse=','), '\\n'); "
            "cat(paste(annotation_metadata(b), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_binary(b)), collapse=','), '\\n'); "
            "cat(class(annomat(b))[[1]], '\\n'); "
            "cat(paste(format(as.vector(M), digits=17), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "weight,score"
    assert lines[1] == "weight meta,score meta"
    assert lines[2] == "0,0"
    assert lines[3] == "dgCMatrix"
    np.testing.assert_allclose(
        _parse_float_csv(lines[4]),
        np.array(
            [
                10.0, 10.0, 20.0, 20.0, 0.0, 30.0, 30.0, 30.0,
                0.5, 0.5, 1.5, 1.5, 0.0, 2.5, 2.5, 2.5,
            ]
        ),
    )


@pytest.mark.r
@skipif_no_rscript
def test_r_load_annotation_group_column_fixture():
    grouped = FIXTURES_DIR / "annotations/grouped.annot"
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"a <- load_annotation({json.dumps(str(grouped))}, ref, has_header = TRUE, group_column = 'group'); "
            f"b <- load_annotation({json.dumps(str(grouped))}, ref, has_header = TRUE, group_column = 4L); "
            "M <- as.matrix(annomat(a)); "
            "meta <- lapply(annotation_metadata(a), jsonlite::fromJSON); "
            "meta_indexed <- lapply(annotation_metadata(b), jsonlite::fromJSON); "
            "cat(paste(annonames(a), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_binary(a)), collapse=','), '\\n'); "
            "cat(paste(as.integer(as.vector(M)), collapse=','), '\\n'); "
            "cat(paste(vapply(meta, function(x) x$group_value, character(1)), collapse=','), '\\n'); "
            "cat(paste(vapply(meta, function(x) x$num_source_intervals, integer(1)), collapse=','), '\\n'); "
            "cat(paste(annonames(b), collapse=','), '\\n'); "
            "cat(paste(vapply(meta_indexed, function(x) as.integer(x$group_column), integer(1)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "coding,regulatory"
    assert lines[1] == "1,1"
    assert lines[2] == "1,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1"
    assert lines[3] == "coding,regulatory"
    assert lines[4] == "2,2"
    assert lines[5] == "coding,regulatory"
    assert lines[6] == "4,4"


@pytest.mark.r
@skipif_no_rscript
def test_r_load_annotation_group_column_overlap_union(tmp_path):
    grouped = tmp_path / "grouped_overlap.annot"
    grouped.write_text(
        "chrom\tstart0\tend0\tgroup\n"
        "1\t99\t250\tcoding\n"
        "1\t150\t401\tcoding\n"
        "X\t99\t301\tregulatory\n"
    )
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"a <- load_annotation({json.dumps(str(grouped))}, ref, has_header = TRUE, group_column = 'group'); "
            "M <- as.matrix(annomat(a)); "
            "meta <- lapply(annotation_metadata(a), jsonlite::fromJSON); "
            "cat(paste(annonames(a), collapse=','), '\\n'); "
            "cat(paste(as.integer(as.vector(M)), collapse=','), '\\n'); "
            "cat(paste(vapply(meta, function(x) x$num_source_intervals, integer(1)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "coding,regulatory"
    assert lines[1] == "1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1"
    assert lines[2] == "2,1"


@pytest.mark.r
@skipif_no_rscript
def test_r_load_annotation_binary_and_batch_sidecar_metadata(tmp_path):
    bed = tmp_path / "raw.bed"
    bed.write_text("1\t99\t200\n")
    binary_meta = tmp_path / "raw.meta"
    binary_meta.write_text("stored\nmetadata\n", encoding="utf-8")
    batch_meta = tmp_path / "anno1.meta"
    batch_meta.write_text("batch\nmetadata\n", encoding="utf-8")

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"single <- load_annotation({json.dumps(str(bed))}, ref, annotation_names = 'renamed', annotation_metadata_path = {json.dumps(str(binary_meta))}); "
            f"batch <- load_annotations(c({json.dumps(str(ANNOTATIONS[0]))}, {json.dumps(str(ANNOTATIONS[1]))}), ref, annotation_metadata_paths = c({json.dumps(str(batch_meta))}, NA_character_)); "
            "cat(paste(annonames(single), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_binary(single)), collapse=','), '\\n'); "
            "cat(as.character(identical(annotation_metadata(single), 'stored\\nmetadata\\n')), '\\n'); "
            "cat(paste(as.integer(as.vector(as.matrix(annomat(single)))), collapse=','), '\\n'); "
            "cat(as.character(identical(annotation_metadata(batch)[[1]], 'batch\\nmetadata\\n')), '\\n'); "
            "cat(as.character(grepl('anno2.bed', annotation_metadata(batch)[[2]], fixed = TRUE)), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.rstrip() for line in result.stdout.splitlines()]
    assert lines[0] == "renamed"
    assert lines[1] == "1"
    assert lines[2] == "TRUE"
    assert lines[3] == "1,1,0,0,0,0,0,0"
    assert lines[4] == "TRUE"
    assert lines[5] == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_load_annotation_validation_and_old_cache_upgrade(tmp_path):
    overlap = tmp_path / "overlap.annot"
    overlap.write_text("1\t99\t200\t1\n1\t150\t250\t2\n")
    nonfinite = tmp_path / "nonfinite.annot"
    nonfinite.write_text("1\t99\t200\tNaN\n")
    wide = tmp_path / "wide.annot"
    wide.write_text("1\t99\t200\t1\t2\n")
    old_cache = tmp_path / "old_annotations.rds"

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            "a <- load_annotations(c('tests/fixtures/annotations/anno1.bed', 'tests/fixtures/annotations/anno2.bed'), ref); "
            f"payload <- list(metadata = list(schema = 'annotations_cache/0.1', n_shards = length(shards(a)), shard_labels = vapply(shards(a), function(s) s$label, character(1)), shard_checksums = vapply(shards(a), function(s) s$reference_checksum, character(1)), shard_start0 = shard_offsets(a)$start0, shard_stop0 = shard_offsets(a)$stop0), annomat = as(annomat(a), 'lgCMatrix'), annonames = annonames(a)); "
            f"saveRDS(payload, {json.dumps(str(old_cache))}); "
            f"loaded <- load_annotations_cache({json.dumps(str(old_cache))}); "
            f"ok1 <- FALSE; tryCatch(load_annotation({json.dumps(str(overlap))}, ref), error = function(e) ok1 <<- grepl('overlap', e$message)); "
            f"ok2 <- FALSE; tryCatch(load_annotation({json.dumps(str(nonfinite))}, ref), error = function(e) ok2 <<- grepl('finite numeric', e$message)); "
            f"ok3 <- FALSE; tryCatch(load_annotation({json.dumps(str(wide))}, ref), error = function(e) ok3 <<- grepl('requires explicit value_columns', e$message)); "
            f"ok4 <- FALSE; tryCatch(load_annotation({json.dumps(str(wide))}, ref, value_columns = 4), error = function(e) ok4 <<- grepl('requires annotation_names', e$message)); "
            f"empty_group <- {json.dumps(str(tmp_path / 'empty_group.annot'))}; "
            "writeLines(c('chrom\\tstart0\\tend0\\tgroup', '1\\t99\\t200\\t'), empty_group); "
            "ok5 <- FALSE; tryCatch(load_annotation(empty_group, ref, has_header = TRUE, group_column = 'group'), error = function(e) ok5 <<- grepl('group_column values must be non-empty', e$message)); "
            f"grouped <- {json.dumps(str(tmp_path / 'grouped.annot'))}; "
            "writeLines(c('chrom\\tstart0\\tend0\\tgroup\\tscore', '1\\t99\\t200\\tcoding\\t1'), grouped); "
            "ok6 <- FALSE; tryCatch(load_annotation(grouped, ref, has_header = TRUE, group_column = 'group', value_columns = 'score'), error = function(e) ok6 <<- grepl('mutually exclusive', e$message)); "
            "ok7 <- FALSE; tryCatch(load_annotation(grouped, ref, has_header = TRUE, group_column = 'group', annotation_names = 'coding'), error = function(e) ok7 <<- grepl('annotation_names is invalid', e$message)); "
            "ok8 <- FALSE; tryCatch(load_annotation(grouped, ref, has_header = TRUE, group_column = 'group', annotation_metadata = 'meta'), error = function(e) ok8 <<- grepl('annotation_metadata is invalid', e$message)); "
            "ok9 <- FALSE; tryCatch(load_annotation(grouped, ref, has_header = TRUE, group_column = 'group', annotation_metadata_path = grouped), error = function(e) ok9 <<- grepl('annotation_metadata_path is invalid', e$message)); "
            "cat(as.character(ok1), as.character(ok2), as.character(ok3), as.character(ok4), as.character(ok5), as.character(ok6), as.character(ok7), as.character(ok8), as.character(ok9), '\\n'); "
            "cat(class(annomat(loaded))[[1]], '\\n'); "
            "cat(paste(as.integer(is_binary(loaded)), collapse=','), '\\n'); "
            "cat(paste(annotation_metadata(loaded), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE", "dgCMatrix", "1,1", ","]


@pytest.mark.r
@skipif_no_rscript
def test_r_annotation_union_uses_reference_compatible_panels():
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            "a <- create_annotation(ref, rep(c(1, 0), length.out = num_snp(ref)), 'a'); "
            "b <- create_annotation(ref, rep(c(0, 1), length.out = num_snp(ref)), 'b'); "
            "u <- union_annotations(a, b); "
            "bad <- b; bad$shards[[1]]$reference_checksum <- 'bad'; "
            "warned <- FALSE; ok <- FALSE; withCallingHandlers(tryCatch(union_annotations(a, bad), error = function(e) ok <<- grepl('compatible reference alignment', e$message)), warning = function(w) warned <<- TRUE); "
            "msg <- NULL; tryCatch(union_annotations(a, bad), error = function(e) msg <<- e$message); "
            "cat(paste(annonames(u), collapse=','), '\\n'); "
            "cat(paste(colnames(annomat(u)), collapse=','), '\\n'); "
            "cat(as.character(ok), '\\n'); "
            "cat(as.character(warned), '\\n'); "
            "cat(as.character(grepl('reference_checksum mismatch', msg)), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["a,b", "a,b", "TRUE", "FALSE", "TRUE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_annotations_allow_leading_comments_but_reject_late_comments(tmp_path):
    leading = tmp_path / "leading.bed"
    leading.write_text("# header\n\n1\t99\t200\n")
    late = tmp_path / "late.bed"
    late.write_text("1\t99\t200\n# late\n")

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}, shards = '1'); "
            f"a <- load_annotations({json.dumps(str(leading))}, ref); "
            f"ok <- FALSE; tryCatch(load_annotations({json.dumps(str(late))}, ref), error = function(e) ok <<- grepl('blank or comment line after BED data row', e$message)); "
            "cat(num_snp(a), num_annot(a), '\\n'); "
            "cat(paste(as.integer(as.vector(as.matrix(annomat(a)))), collapse=','), '\\n'); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["5 1", "1,1,0,0,0", "TRUE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_phase2_rejects_malformed_inputs(tmp_path):
    bad_sumstats = tmp_path / "bad.tsv"
    bad_sumstats.write_text("chr\tbp\ta1\ta2\tp\n1\t100\tA\tG\t1.5\n")
    bad_nan_sumstats = tmp_path / "bad_nan.tsv"
    bad_nan_sumstats.write_text("chr\tbp\ta1\ta2\tp\n1\t100\tA\tG\tNaN\n")
    bad_allele_sumstats = tmp_path / "bad_allele.tsv"
    bad_allele_sumstats.write_text("chr\tbp\ta1\ta2\tp\n1\t100\ta\tG\t0.1\n")
    same_allele_sumstats = tmp_path / "same_allele.tsv"
    same_allele_sumstats.write_text("chr\tbp\ta1\ta2\tp\n1\t100\tA\tA\t0.1\n")
    ignored_chr_sumstats = tmp_path / "ignored_chr.tsv"
    ignored_chr_sumstats.write_text("chr\tbp\ta1\ta2\tp\nY\t1\tN\tN\tNaN\nMT\t1\tN\tN\tNaN\n1\t100\tA\tG\t0.1\n")
    bad_bed = tmp_path / "bad.bed"
    bad_bed.write_text("1\t100\n")
    bad_cache = tmp_path / "bad_sumstats_cache.rds"
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"saveRDS(list(metadata = list(schema = 'sumstats_cache/0.1', n_shards = 0L, shard_labels = character(), shard_checksums = character(), shard_start0 = numeric(), shard_stop0 = numeric()), logpvec = numeric()), {json.dumps(str(bad_cache))}); "
            f"ok1 <- FALSE; tryCatch(load_sumstats({json.dumps(str(bad_sumstats))}, ref), error = function(e) ok1 <<- grepl('p must be finite numeric', e$message)); "
            f"missing <- suppressWarnings(load_sumstats({json.dumps(str(bad_nan_sumstats))}, ref)); ok1_nan <- all(is.nan(logpvec(missing))) && !any(is_present(missing)); "
            f"ok1_allele <- FALSE; tryCatch(load_sumstats({json.dumps(str(bad_allele_sumstats))}, ref), error = function(e) ok1_allele <<- grepl('uppercase DNA bases', e$message)); "
            f"ok1_same <- FALSE; tryCatch(load_sumstats({json.dumps(str(same_allele_sumstats))}, ref), error = function(e) ok1_same <<- grepl('a1 and a2 must differ', e$message)); "
            f"ignored <- load_sumstats({json.dumps(str(ignored_chr_sumstats))}, ref); ok1_ignored <- sum(is_present(ignored)) == 1L; "
            f"ok2 <- FALSE; tryCatch(load_annotations({json.dumps(str(bad_bed))}, ref), error = function(e) ok2 <<- grepl('at least 3', e$message)); "
            "ok3 <- FALSE; tryCatch(load_annotations('', ref), error = function(e) ok3 <<- grepl('bed_paths must be a non-empty character vector', e$message)); "
            "ok4 <- FALSE; tryCatch(load_annotations(NA_character_, ref), error = function(e) ok4 <<- grepl('bed_paths must be a non-empty character vector', e$message)); "
            f"ok5 <- FALSE; tryCatch(load_sumstats_cache({json.dumps(str(bad_cache))}), error = function(e) ok5 <<- grepl('n_shards must be at least 1', e$message)); "
            "stopifnot(ok1, ok1_nan, ok1_allele, ok1_same, ok1_ignored, ok2, ok3, ok4, ok5)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_swapped_allele_alignment_and_duplicate_key_error(tmp_path):
    swapped = tmp_path / "swapped.tsv"
    swapped.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\t0.01\n"
        "1\t200\tT\tC\t0.02\n"
    )
    duplicate = tmp_path / "duplicate.tsv"
    duplicate.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\t0.01\n"
        "1\t100\tA\tG\t0.02\n"
    )
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}, shards = '1'); "
            "msg <- NULL; "
            f"s <- withCallingHandlers(load_sumstats({json.dumps(str(swapped))}, ref), warning = function(w) {{ msg <<- paste(msg, conditionMessage(w), sep='|'); invokeRestart('muffleWarning') }}); "
            "stopifnot(!grepl('would match the reference', paste(msg, collapse='|'))); "
            "stopifnot(identical(as.integer(is_present(s)), c(1L, 1L, 0L, 0L, 0L))); "
            f"ok <- FALSE; tryCatch(load_sumstats({json.dumps(str(duplicate))}, ref), error = function(e) ok <<- grepl('duplicate sumstats/reference matching key', e$message)); "
            "stopifnot(ok)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_variant_matching_preserves_integer64_hashes_above_double_precision(tmp_path):
    ref_path = tmp_path / "high_hash.bim"
    ref_path.write_text(
        "1\trsHigh1\t0\t100\tACGT\tCCCCCCCC\n"
        "1\trsHigh2\t0\t200\tTGCA\tGGGGGGGG\n"
    )
    sumstats_path = tmp_path / "high_hash.tsv"
    sumstats_path.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tACGT\tCCCCCCCC\t0.01\n"
        "1\t200\tTGCA\tGGGGGGGG\t0.2\n"
    )
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(ref_path))}); "
            "limit <- bit64::as.integer64('9007199254740992'); "
            "stopifnot(all(a1_hash64(ref) > limit), all(a2_hash64(ref) > limit)); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(sumstats_path))}, ref)); "
            "h <- .allele_hash64(c('ACGT', 'TGCA', 'CCCCCCCC', 'GGGGGGGG')); "
            "direct <- .match_variant_keys(c(100L, 200L), h[1:2], h[3:4], c(200L, 100L), h[c(2, 1)], h[c(4, 3)], '1', 'sumstats'); "
            "dup_ok <- FALSE; tryCatch(.match_variant_keys(c(100L), h[1], h[3], c(100L, 100L), h[c(1, 1)], h[c(3, 3)], '1', 'variant'), error = function(e) dup_ok <<- grepl('duplicate variant/reference matching key', e$message)); "
            "cat(paste(as.integer(is_present(s)), collapse=','), '\\n'); "
            "cat(paste(format(logpvec(s), digits=17), collapse=','), '\\n'); "
            "cat(paste(direct$source_rows, collapse=','), '\\n'); "
            "cat(as.character(dup_ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "1,1"
    np.testing.assert_allclose(_parse_float_csv(lines[1]), np.array([2.0, -math.log10(0.2)]))
    assert lines[2] == "2,1"
    assert lines[3] == "TRUE"


# Hashes of these alleles are >= 0x7FF0000000000000, so their raw bits read as
# NaN when an integer64 is treated as a double. Two such values are needed to
# expose unstable grouping under base::order().
NAN_RANGE_A2_X = "GTAGAA"
NAN_RANGE_A2_Y = "GTAGAC"


@pytest.mark.r
@skipif_no_rscript
def test_r_variant_keys_sort_hashes_in_double_nan_bit_range(tmp_path):
    x, y = NAN_RANGE_A2_X, NAN_RANGE_A2_Y
    ref_path = tmp_path / "ref.bim"
    ref_path.write_text(f"1\trs1\t0\t100\tT\t{x}\n1\trs2\t0\t100\tT\t{y}\n")
    dup_bim = tmp_path / "dup.bim"
    dup_bim.write_text(
        f"1\trs1\t0\t100\tT\t{x}\n1\trs2\t0\t100\tT\t{y}\n1\trs3\t0\t100\tT\t{x}\n"
    )
    exact = tmp_path / "exact.tsv"
    exact.write_text(f"chr\tbp\ta1\ta2\tp\n1\t100\tT\t{y}\t0.001\n1\t100\tT\t{x}\t0.01\n")
    dup_sumstats = tmp_path / "dup.tsv"
    dup_sumstats.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        f"1\t100\tT\t{x}\t0.01\n1\t100\tT\t{y}\t0.1\n1\t100\tT\t{x}\t0.5\n"
    )
    result = run_rscript(
        _source_phase2_script(
            f"h <- .allele_hash64(c({json.dumps(x)}, {json.dumps(y)}, 'T')); "
            "stopifnot(all(h[1:2] >= bit64::as.integer64('9218868437227405312'))); "
            "fails <- function(f) tryCatch({ f(); FALSE }, error = function(e) grepl('duplicate', e$message)); "
            "stopifnot(fails(function() .check_unique_variant_keys(rep(1L, 3L), h[c(3, 3, 3)], h[c(1, 2, 1)], 't'))); "
            "stopifnot(fails(function() .check_unique_variant_keys(rep(1L, 3L), h[c(1, 2, 1)], h[c(3, 3, 3)], 't'))); "
            f"stopifnot(fails(function() load_reference({json.dumps(str(dup_bim))}))); "
            f"ref <- load_reference({json.dumps(str(ref_path))}); "
            f"stopifnot(fails(function() suppressWarnings(load_sumstats({json.dumps(str(dup_sumstats))}, ref)))); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(exact))}, ref)); "
            "stopifnot(all.equal(logpvec(s), c(2, 3)))"
        )
    )
    assert result.returncode == 0, result.stderr
