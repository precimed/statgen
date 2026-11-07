"""R phase 1 acceptance tests: reference panel loading and cache behavior."""

import hashlib
import json

import pytest

from tests.conftest import FIXTURES_DIR, R_PACKAGE_DIR, run_rscript, skipif_no_rscript
from statgen.reference import load_reference

SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
NONSHARDED = FIXTURES_DIR / "reference/nonsharded/all.bim"

CHR1_ROWS = [
    ("1", "rs1001", 0, 100, "A", "G"),
    ("1", "rs1002", 0, 200, "C", "T"),
    ("1", "rs1003", 0, 300, "A", "C"),
    ("1", "rs1004", 0, 400, "G", "A"),
    ("1", "rs1005", 0, 500, "T", "C"),
]
CHRX_ROWS = [
    ("X", "rsX001", 0, 100, "A", "G"),
    ("X", "rsX002", 0, 200, "C", "T"),
    ("X", "rsX003", 0, 300, "A", "C"),
]


def _checksum(rows):
    text = "".join(f"{r[0]}:{r[3]}:{r[4]}:{r[5]}\n" for r in rows)
    return hashlib.md5(text.encode()).hexdigest()


CHR1_CHECKSUM = _checksum(CHR1_ROWS)
CHRX_CHECKSUM = _checksum(CHRX_ROWS)


def _source_reference_script(expr: str) -> str:
    files = [
        "utils.R",
        "hash.R",
        "variant_match.R",
        "bfile_utils.R",
        "reference.R",
    ]
    sources = " ".join(
        f"source({json.dumps(str(R_PACKAGE_DIR / 'R' / filename))});"
        for filename in files
    )
    return f"{sources} {expr}"


@pytest.mark.r
@skipif_no_rscript
def test_r_bim_reader_preserves_character_fields(tmp_path):
    path = tmp_path / "characters.bim"
    rows = [
        ("1", "000123", "0.00", "000100", "AC", "G"),
        ("1", "1e6", "1e-3", "200", "T", "CA"),
        ("X", "NA", "0", "300", "G", "A"),
    ]
    path.write_text("".join("\t".join(row) + "\n" for row in rows))
    result = run_rscript(_source_reference_script(
        f'df <- .parse_bim({json.dumps(str(path))}); '
        'stopifnot(identical(df$chr, c("1", "1", "X")), '
        'identical(df$snp, c("000123", "1e6", "NA")), '
        'identical(df$cm, c("0.00", "1e-3", "0")), '
        'identical(df$bp, c(100L, 200L, 300L)), '
        'identical(df$a1, c("AC", "T", "G")), '
        'identical(df$a2, c("G", "CA", "A")))'
    ))
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
@pytest.mark.parametrize("row", ["1 rs1 0 100 A", "1 rs1 0 100 A G extra"])
def test_r_bim_reader_requires_six_columns(tmp_path, row):
    path = tmp_path / "columns.bim"
    path.write_text((row + "\n") * 3)
    result = run_rscript(_source_reference_script(
        f'failed <- tryCatch({{.parse_bim({json.dumps(str(path))}); FALSE}}, '
        'error=function(e) grepl("expected 6 whitespace-delimited columns", conditionMessage(e), fixed=TRUE)); '
        'stopifnot(failed)'
    ))
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_bim_reader_rejects_short_row_mid_file(tmp_path):
    # A malformed row must fail the load rather than silently truncate the file.
    path = tmp_path / "malformed.bim"
    rows = [f"1\trs{i}\t0\t{i}\tA\tG" for i in range(1, 31)]
    rows.insert(15, "1\tshort\t0\t31\tA")
    path.write_text("\n".join(rows) + "\n")
    result = run_rscript(_source_reference_script(
        f'failed <- tryCatch({{.parse_bim({json.dumps(str(path))}); FALSE}}, error=function(e) TRUE); '
        'stopifnot(failed)'
    ))
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_load_reference_sharded_matches_fixture_contract():
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(SHARDED))}); "
            "cat(num_snp(ref), '\\n'); "
            "cat(paste(vapply(shards(ref), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(vapply(shards(ref), function(s) s$checksum, character(1)), collapse=','), '\\n'); "
            "cat(paste(bp(ref), collapse=','), '\\n'); "
            "cat(paste(chr(ref), collapse=','), '\\n'); "
            "cat(as.character(inherits(a1_hash64(ref), 'integer64')), '\\n'); "
            "cat(as.character(validate_checksums(ref)), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == [
        "8",
        "1,X",
        f"{CHR1_CHECKSUM},{CHRX_CHECKSUM}",
        "100,200,300,400,500,100,200,300",
        "1,1,1,1,1,X,X,X",
        "TRUE",
        "TRUE",
    ]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_hash_checksum_and_offset_parity_with_python():
    py_ref = load_reference(SHARDED)
    py_a1_hash = ",".join(str(int(x)) for x in py_ref.a1_hash64)
    py_a2_hash = ",".join(str(int(x)) for x in py_ref.a2_hash64)
    py_checksums = ",".join(s.checksum for s in py_ref.shards)
    py_offsets = ";".join(
        f"{row['shard_label']}:{row['start0']}:{row['stop0']}"
        for row in py_ref.shard_offsets
    )

    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(SHARDED))}); "
            "off <- shard_offsets(ref); "
            "cat(paste(as.character(a1_hash64(ref)), collapse=','), '\\n'); "
            "cat(paste(as.character(a2_hash64(ref)), collapse=','), '\\n'); "
            "cat(paste(vapply(shards(ref), function(s) s$checksum, character(1)), collapse=','), '\\n'); "
            "cat(paste(paste(off$shard_label, off$start0, off$stop0, sep=':'), collapse=';'), '\\n'); "
            "cat(typeof(off$start0), typeof(off$stop0), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == [py_a1_hash, py_a2_hash, py_checksums, py_offsets, "integer integer"]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_nonsharded_split_and_cache_subset(tmp_path):
    cache = tmp_path / "ref_cache.rds"
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(NONSHARDED))}); "
            f"save_reference_cache(ref, {json.dumps(str(cache))}); "
            f"loaded <- load_reference_cache({json.dumps(str(cache))}, shards = 'X'); "
            "cat(num_snp(ref), '\\n'); "
            "cat(paste(vapply(shards(ref), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(num_snp(loaded), '\\n'); "
            "cat(paste(snp(loaded), collapse=','), '\\n'); "
            "cat(as.character(validate_checksums(loaded)), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["8", "1,X", "3", "rsX001,rsX002,rsX003", "TRUE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_cache_rejects_zero_shards(tmp_path):
    cache = tmp_path / "zero_shards.rds"
    result = run_rscript(
        _source_reference_script(
            "payload <- list("
            "metadata = list(schema = .reference_cache_schema, n_shards = 0L, "
            "shard_labels = character(), shard_checksums = character(), "
            "shard_start0 = integer(), shard_stop0 = integer()), "
            "bp = integer(), snp_text_by_shard = character(), "
            "a1_text_by_shard = character(), a2_text_by_shard = character(), "
            "a1_hash64 = bit64::as.integer64(integer()), "
            "a2_hash64 = bit64::as.integer64(integer())); "
            f"saveRDS(payload, {json.dumps(str(cache))}); "
            "ok <- FALSE; "
            f"tryCatch(load_reference_cache({json.dumps(str(cache))}), "
            "error = function(e) ok <<- grepl('n_shards must be at least 1', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_select_shards_source_loaded_panel():
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(SHARDED))}); "
            "sub <- select_shards(ref, 'X'); "
            "off <- shard_offsets(sub); "
            "cat(num_snp(sub), '\\n'); "
            "cat(paste(vapply(shards(sub), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(bp(sub), collapse=','), '\\n'); "
            "cat(paste(paste(off$shard_label, off$start0, off$stop0, sep=':'), collapse=';'), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["3", "X", "100,200,300", "X:0:3"]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_is_object_compatible():
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(SHARDED))}); "
            f"same <- load_reference({json.dumps(str(SHARDED))}); "
            "sub <- select_shards(ref, 'X'); "
            "cat(as.character(is_object_compatible(ref, same)), '\\n'); "
            "cat(as.character(suppressWarnings(is_object_compatible(ref, sub))), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["TRUE", "FALSE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_is_object_compatible_warns_with_label_pair():
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(SHARDED))}); "
            "bad <- ref; bad$shards[[2]]$label <- '2'; "
            "msg <- NULL; "
            "ok <- withCallingHandlers(is_object_compatible(ref, bad), warning = function(w) { msg <<- conditionMessage(w); invokeRestart('muffleWarning') }); "
            "stopifnot(!ok); "
            "cat(msg, '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert "reference X, object 2" in result.stdout


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_is_object_compatible_wrong_dispatch_object_returns_false():
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(SHARDED))}); "
            "msg <- NULL; "
            "ok <- withCallingHandlers(is_object_compatible(list(), ref), warning = function(w) { msg <<- conditionMessage(w); invokeRestart('muffleWarning') }); "
            "cat(as.character(ok), '\\n'); "
            "cat(msg, '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == [
        "FALSE",
        "statgen: is_object_compatible: reference must be a ReferencePanel",
    ]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_rejects_bad_bim_inputs(tmp_path):
    bad = tmp_path / "bad.bim"
    bad.write_text("1\trs1\t0\t100\tA\tA\n")
    result = run_rscript(
        _source_reference_script(
            f"ok <- FALSE; tryCatch(load_reference({json.dumps(str(bad))}), error = function(e) ok <<- grepl('a1 and a2 must differ', e$message)); "
            "stopifnot(ok)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
@pytest.mark.parametrize(
    ("contents", "pattern"),
    [
        (
            "1\trs1\t0\t200\tA\tG\n"
            "1\trs2\t0\t100\tC\tT\n",
            "chr_rank, bp",
        ),
        (
            "X\trsx\t0\t100\tA\tG\n"
            "1\trs1\t0\t200\tC\tT\n",
            "chr_rank, bp",
        ),
        (
            "1\trs1\t0\t100\tA\tG\n"
            "1\trs2\t0\t100\tA\tG\n",
            "duplicate",
        ),
    ],
)
def test_r_reference_rejects_structural_bim_violations(tmp_path, contents, pattern):
    bad = tmp_path / "bad.bim"
    bad.write_text(contents)
    result = run_rscript(
        _source_reference_script(
            f"ok <- FALSE; tryCatch(load_reference({json.dumps(str(bad))}), error = function(e) ok <<- grepl({json.dumps(pattern)}, e$message)); "
            "stopifnot(ok)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_ignores_y_mt_rows(tmp_path):
    path = tmp_path / "with_y_mt.bim"
    path.write_text(
        "1\trs1\t0\t100\tA\tG\n"
        "Y\trsY\t0\t150\tA\tC\n"
        "MT\trsMT\t0\t180\tG\tA\n"
        "X\trsX\t0\t200\tC\tT\n"
    )
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(path))}); "
            "cat(num_snp(ref), '\\n'); "
            "cat(paste(vapply(shards(ref), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(chr(ref), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["2", "1,X", "1,X"]


@pytest.mark.r
@skipif_no_rscript
def test_r_reference_variant_type_accessors(tmp_path):
    path = tmp_path / "variant_types.bim"
    path.write_text(
        "1\trs1\t0\t100\tA\tT\n"
        "1\trs2\t0\t200\tT\tA\n"
        "1\trs3\t0\t300\tC\tG\n"
        "1\trs4\t0\t400\tA\tC\n"
        "1\trs5\t0\t500\tAC\tG\n"
        "1\trs6\t0\t600\tA\tAT\n"
    )
    result = run_rscript(
        _source_reference_script(
            f"ref <- load_reference({json.dumps(str(path))}); "
            "cat(paste(as.integer(is_single_nucleotide_variant(ref)), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_strand_ambiguous(ref)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["1,1,1,1,0,0", "1,1,1,0,0,0"]
