"""R phase 3 acceptance tests: genotype panels."""

import json

import numpy as np
import pytest

from statgen.genotype import load_genotype
from statgen.reference import load_reference
from tests.conftest import (
    FIXTURES_DIR,
    GENOTYPE_CHR1_CALLS,
    R_PACKAGE_DIR,
    copy_genotype_shard_files,
    copy_sharded_genotype,
    run_rscript,
    skipif_no_rscript,
    write_plink_bed_calls,
)

REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
REF_NONSHARDED = FIXTURES_DIR / "reference/nonsharded/all.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


def _source_phase3_script(expr: str) -> str:
    files = [
        "utils.R",
        "verbosity.R",
        "hash.R",
        "variant_match.R",
        "bfile_utils.R",
        "reference.R",
        "genotype.R",
    ]
    sources = " ".join(
        f"source({json.dumps(str(R_PACKAGE_DIR / 'R' / filename))});"
        for filename in files
    )
    return f"{sources} {expr}"


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_sharded_metadata_matches_python_fixture():
    reference = load_reference(REF_SHARDED)
    py = load_genotype(G_SHARDED, reference)

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}); "
            f"g <- load_genotype({json.dumps(str(G_SHARDED))}, ref); "
            "cat(source_layout(g), '\\n'); "
            "cat(num_snp(g), num_sample(g), '\\n'); "
            "cat(paste(vapply(shards(g), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(vapply(shards(g), function(s) s$chr, character(1)), collapse=','), '\\n'); "
            "cat(paste(source_row0(g), collapse=','), '\\n'); "
            "cat(paste(format(ploidy_male(g), digits=17), collapse=','), '\\n'); "
            "cat(paste(format(ploidy_female(g), digits=17), collapse=','), '\\n'); "
            "cat(paste(fid(g), collapse=','), '\\n'); "
            "cat(paste(iid(g), collapse=','), '\\n'); "
            "cat(paste(sex(g), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_male(g)), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_subject_present(g, 'X')), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == py.source_layout
    assert lines[1] == f"{py.num_snp} {py.num_sample}"
    assert lines[2] == ",".join(s.label for s in py.shards)
    assert lines[3] == ",".join(s.chr for s in py.shards)
    assert lines[4] == ",".join(str(int(x)) for x in py.source_row0)
    np.testing.assert_allclose([float(x) for x in lines[5].split(",")], py.ploidy_male, equal_nan=True)
    np.testing.assert_allclose([float(x) for x in lines[6].split(",")], py.ploidy_female, equal_nan=True)
    assert lines[7] == ",".join(py.fid)
    assert lines[8] == ",".join(py.iid)
    assert lines[9] == ",".join(str(int(x)) for x in py.sex)
    assert lines[10] == ",".join(str(int(x)) for x in py.is_male)
    assert lines[11] == ",".join(str(int(x)) for x in py.is_subject_present("X"))


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_nonsharded_metadata_and_cache_subset(tmp_path):
    cache = tmp_path / "genotype.rds"
    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_NONSHARDED))}); "
            f"g <- load_genotype({json.dumps(str(G_NONSHARDED))}, ref); "
            f"save_genotype_cache(g, {json.dumps(str(cache))}); "
            f"h <- load_genotype_cache({json.dumps(str(cache))}, shards = 'X'); "
            "cat(source_layout(g), '\\n'); "
            "cat(paste(vapply(shards(g), function(s) s$bed_path, character(1)), collapse=','), '\\n'); "
            "cat(paste(source_row0(g), collapse=','), '\\n'); "
            "cat(source_layout(h), num_snp(h), num_sample(h), '\\n'); "
            "cat(paste(vapply(shards(h), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(source_row0(h), collapse=','), '\\n'); "
            "cat(paste(fid(h), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "non_sharded"
    assert lines[1] == f"{FIXTURES_DIR / 'genotype/nonsharded/all.bed'},{FIXTURES_DIR / 'genotype/nonsharded/all.bed'}"
    assert lines[2] == "0,1,2,3,4,5,6,7"
    assert lines[3] == "non_sharded 3 4"
    assert lines[4] == "X"
    assert lines[5] == "5,6,7"
    assert lines[6] == "FAM1,FAM1,FAM2,FAM2"


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_select_shards_and_reference_compatibility():
    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}); "
            f"g <- load_genotype({json.dumps(str(G_SHARDED))}, ref); "
            "gx <- select_shards(g, 'X'); "
            "bad <- g; bad$shards[[1]]$reference_checksum <- 'bad'; "
            "cat(num_snp(gx), num_sample(gx), source_layout(gx), '\\n'); "
            "cat(paste(vapply(shards(gx), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(as.character(is_object_compatible(ref, g)), '\\n'); "
            "cat(as.character(suppressWarnings(is_object_compatible(ref, bad))), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert [line.strip() for line in result.stdout.strip().splitlines()] == [
        "3 4 sharded",
        "X",
        "TRUE",
        "FALSE",
    ]


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_nonsharded_absent_chrx_does_not_warn_about_ploidy(tmp_path):
    for suffix in (".bim", ".fam", ".bed"):
        (tmp_path / f"all{suffix}").write_bytes((FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes())

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}); "
            "messages <- character(); "
            "withCallingHandlers("
            f"g <- load_genotype({json.dumps(str(tmp_path / 'all'))}, ref), "
            "warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart('muffleWarning') }); "
            "cat(paste(messages, collapse='|'), '\\n'); "
            "cat(paste(as.integer(is_present(g)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.splitlines()]
    assert "chrX genotype source has no .ploidy" not in lines[0]
    assert "no source BIM rows for requested reference shard" in lines[0]
    assert lines[1] == "1,1,1,1,1,0,0,0"


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_fetch_rejects_absent_snp_and_bad_nonsharded_override(tmp_path):
    for suffix in (".bim", ".fam", ".bed"):
        (tmp_path / f"all{suffix}").write_bytes((FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes())

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}); "
            f"g_absent <- suppressWarnings(load_genotype({json.dumps(str(tmp_path / 'all'))}, ref)); "
            f"g_flat <- load_genotype({json.dumps(str(G_NONSHARDED))}, load_reference({json.dumps(str(REF_NONSHARDED))})); "
            "ok1 <- FALSE; tryCatch(fetch_genotypes_int8(g_absent, 6), error = function(e) ok1 <<- grepl('not present', e$message)); "
            f"ok2 <- FALSE; tryCatch(fetch_genotypes_int8(g_flat, 1, bed_path = {json.dumps(str(FIXTURES_DIR / 'genotype/sharded/@.bed'))}), error = function(e) ok2 <<- grepl('@ override incompatible', e$message)); "
            "cat(as.character(ok1), as.character(ok2), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_fetch_decodes_bed_order_missing_and_ploidy(tmp_path):
    for suffix in (".bim", ".fam"):
        (tmp_path / f"chr1{suffix}").write_bytes((FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes())
    write_plink_bed_calls(tmp_path / "chr1.bed", GENOTYPE_CHR1_CALLS)

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}, shards = '1'); "
            f"g <- load_genotype({json.dumps(str(tmp_path / 'chr1'))}, ref); "
            "gi <- fetch_genotypes_int8(g, c(2, 1, 2, 5)); "
            "gf <- fetch_genotypes(g, c(2, 1, 2, 5)); "
            "gf[is.nan(gf)] <- -9; "
            "cat(paste(as.integer(gi), collapse=','), '\\n'); "
            "cat(paste(as.numeric(gf), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "0,1,2,-1,2,-1,1,0,0,1,2,-1,-1,2,0,1"
    assert lines[1] == "0,1,2,-9,2,-9,1,0,0,1,2,-9,-9,2,0,1"


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_chrx_subset_and_ploidy_scaled(tmp_path):
    copy_sharded_genotype(tmp_path)
    (tmp_path / "X.fam").write_text(
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t1\t-9\n"
    )
    write_plink_bed_calls(
        tmp_path / "X.bed",
        np.array([[0, 1], [2, -1], [1, 0]], dtype=np.int8),
    )

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}); "
            f"g <- load_genotype({json.dumps(str(tmp_path / '@'))}, ref); "
            "gi <- fetch_genotypes_int8(g, 6); "
            "gi_empty <- fetch_genotypes_int8(g, integer(0)); "
            "scaled <- fetch_genotypes(g, 6, haploid_mode = 'ploidy_scaled'); "
            "scaled[is.nan(scaled)] <- -9; "
            "cat(paste(as.integer(is_subject_present(g, 'X')), collapse=','), '\\n'); "
            "cat(paste(shards(g)[[2]]$source_subject_row0, collapse=','), '\\n'); "
            "cat(paste(as.integer(gi), collapse=','), '\\n'); "
            "cat(paste(dim(gi_empty), collapse='x'), typeof(gi_empty), '\\n'); "
            "cat(paste(as.numeric(scaled), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["1,0,1,0", "1,-1,0,-1", "1,-1,0,-1", "4x0 integer", "0.5,-9,0,-9"]


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_chrx_only_reference_and_unknown_sex_scaled_error(tmp_path):
    copy_genotype_shard_files(tmp_path, "X")
    (tmp_path / "X.fam").write_text(
        "FAM1\tIND1\t0\t0\t0\t-9\n"
        "FAM1\tIND2\t0\t0\t2\t-9\n"
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM2\tIND4\t0\t0\t2\t-9\n"
    )

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}, shards = 'X'); "
            f"g <- load_genotype({json.dumps(str(tmp_path / 'X'))}, ref); "
            "ok <- FALSE; tryCatch(fetch_genotypes(g, 1, haploid_mode = 'ploidy_scaled'), error = function(e) ok <<- grepl('requires known FAM sex', e$message)); "
            "cat(source_layout(g), num_snp(g), num_sample(g), '\\n'); "
            "cat(paste(as.integer(is_subject_present(g, 'X')), collapse=','), '\\n'); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert [line.strip() for line in result.stdout.strip().splitlines()] == [
        "non_sharded 3 4",
        "1,1,1,1",
        "TRUE",
    ]


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_ploidy_scaled_unknown_sex_scales_only_equal_ploidy_snps():
    result = run_rscript(
        _source_phase3_script(
            "shard <- structure(list(label = 'X', num_snp = 2L, "
            "subject_present = TRUE, ploidy_male = c(0, 2), ploidy_female = c(2, 2)), "
            "class = 'GenotypeShard'); "
            "panel <- structure(list(sex = 0L, shards = list(shard), num_snp = 2L, "
            "shard_offsets = data.frame(shard_label = 'X', start0 = 0L, stop0 = 2L)), "
            "class = 'GenotypePanel'); "
            "out <- .apply_ploidy_scaled(panel, matrix(c(Inf, 2), nrow = 1L), c(1L, 2L)); "
            "ok_error <- FALSE; "
            "tryCatch(.apply_ploidy_scaled(panel, matrix(c(1, 2), nrow = 1L), c(1L, 2L)), "
            "error = function(e) ok_error <<- grepl('requires known FAM sex', e$message)); "
            "cat(as.character(is.infinite(out[1, 1])), out[1, 2], as.character(ok_error), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE 2 TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_cache_load_is_lazy_and_override_fetches(tmp_path):
    copy_sharded_genotype(tmp_path)
    cache = tmp_path / "genotype.rds"
    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}); "
            f"g <- load_genotype({json.dumps(str(tmp_path / '@'))}, ref); "
            f"save_genotype_cache(g, {json.dumps(str(cache))}); "
            f"invisible(file.remove({json.dumps(str(tmp_path / '1.bed'))}, {json.dumps(str(tmp_path / 'X.bed'))})); "
            f"h <- load_genotype_cache({json.dumps(str(cache))}); "
            "ok1 <- FALSE; tryCatch(fetch_genotypes_int8(h, 1), error = function(e) ok1 <<- grepl('File not found', e$message)); "
            f"gi <- fetch_genotypes_int8(h, c(1, 6), bed_path = {json.dumps(str(FIXTURES_DIR / 'genotype/sharded/@.bed'))}); "
            "cat(num_snp(h), '\\n'); "
            "cat(as.character(ok1), '\\n'); "
            "cat(paste(as.integer(gi), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["8", "TRUE", "2,2,2,2,2,2,2,2"]


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_validation_failures(tmp_path):
    copy_genotype_shard_files(tmp_path, "1")
    (tmp_path / "1.fam").write_text("FAM1\tIND1\t0\t0\t3\t-9\n")
    bad_cache = tmp_path / "bad_genotype.rds"
    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}, shards = '1'); "
            f"ok1 <- FALSE; tryCatch(load_genotype({json.dumps(str(tmp_path / '1'))}, ref), error = function(e) ok1 <<- grepl('FAM sex', e$message)); "
            f"saveRDS(list(metadata = list(schema = 'genotype_cache/0.1', n_shards = 0L)), {json.dumps(str(bad_cache))}); "
            f"ok2 <- FALSE; tryCatch(load_genotype_cache({json.dumps(str(bad_cache))}), error = function(e) ok2 <<- grepl('n_shards', e$message)); "
            "stopifnot(ok1, ok2)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_fam_empty_field_message_and_zero_parent_placeholders(tmp_path):
    copy_genotype_shard_files(tmp_path, "1")
    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}, shards = '1'); "
            f"g <- load_genotype({json.dumps(str(tmp_path / '1'))}, ref); "
            "cat(paste(father_id(g), collapse=','), '\\n'); "
            f"bad <- {json.dumps(str(tmp_path / 'bad.fam'))}; "
            "writeLines('FAM1 IND1 0 1 -9', bad); "
            f"invisible(file.copy({json.dumps(str(tmp_path / '1.bim'))}, {json.dumps(str(tmp_path / 'bad.bim'))}, overwrite = TRUE)); "
            f"invisible(file.copy({json.dumps(str(tmp_path / '1.bed'))}, {json.dumps(str(tmp_path / 'bad.bed'))}, overwrite = TRUE)); "
            "ok <- FALSE; "
            f"tryCatch(load_genotype({json.dumps(str(tmp_path / 'bad'))}, ref), "
            "error = function(e) ok <<- grepl('expected 6 whitespace-delimited columns', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["0,0,0,0", "TRUE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_duplicate_fid_iid_rejected(tmp_path):
    copy_genotype_shard_files(tmp_path, "1")
    fam = tmp_path / "1.fam"
    first = fam.read_text().splitlines()[0]
    fam.write_text(fam.read_text() + first + "\n")

    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(REF_SHARDED))}, shards = '1'); "
            f"ok <- FALSE; tryCatch(load_genotype({json.dumps(str(tmp_path / '1'))}, ref), error = function(e) ok <<- grepl('duplicate FAM subject pair', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_bed_helpers_use_double_arithmetic_beyond_32bit_offsets(tmp_path):
    sparse_bed = tmp_path / "large_sparse.bed"
    result = run_rscript(
        _source_phase3_script(
            "num_sample <- 500000; "
            "source_row0 <- 17000; "
            "bytes_per_snp <- floor((num_sample + 3) / 4); "
            "expected_size <- .expected_bed_size(num_sample, 20000); "
            f"con <- file({json.dumps(str(sparse_bed))}, 'w+b'); "
            "writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con); "
            "invisible(seek(con, where = 3 + source_row0 * bytes_per_snp, origin = 'start', rw = 'write')); "
            "writeBin(as.raw(rep(0, bytes_per_snp)), con); "
            "close(con); "
            f"g <- .read_bed_rows_int8_codes({json.dumps(str(sparse_bed))}, source_row0, num_sample); "
            "cat(format(expected_size, scientific = FALSE), '\\n'); "
            "cat(dim(g), '\\n'); "
            "cat(g[1, 1], g[num_sample, 1], '\\n')"
        ),
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["2500000003", "500000 1", "2 2"]


@pytest.mark.r
@skipif_no_rscript
def test_r_genotype_matches_hashes_in_double_nan_bit_range(tmp_path):
    # Both a2 hashes are >= 0x7FF0000000000000; see test_r_phase2.py.
    bim = "1\trs1\t0\t100\tT\tGTAGAA\n1\trs2\t0\t100\tT\tGTAGAC\n"
    (tmp_path / "ref.bim").write_text(bim)
    (tmp_path / "g.bim").write_text(bim)
    (tmp_path / "g.fam").write_text("f1 s1 0 0 1 -9\nf2 s2 0 0 2 -9\n")
    (tmp_path / "g.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00, 0x00]))
    result = run_rscript(
        _source_phase3_script(
            f"ref <- load_reference({json.dumps(str(tmp_path / 'ref.bim'))}); "
            f"g <- load_genotype({json.dumps(str(tmp_path / 'g'))}, ref); "
            "stopifnot(identical(is_present(g), c(TRUE, TRUE)))"
        )
    )
    assert result.returncode == 0, result.stderr
