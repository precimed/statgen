"""Cache BED relocation and availability across R and Octave."""
import json
import shutil

import pytest

from tests.conftest import FIXTURES_DIR, run_octave, run_rscript, skipif_no_octave, skipif_no_rscript
from tests.test_r_phase3 import _source_phase3_script


@pytest.mark.parametrize("sharded", [False, True])
@pytest.mark.parametrize("runtime", ["r", "octave"])
def test_cache_paths(tmp_path, sharded, runtime):
    if runtime == "r" and skipif_no_rscript.args[0]:
        pytest.skip("Rscript unavailable")
    if runtime == "octave" and skipif_no_octave.args[0]:
        pytest.skip("Octave unavailable")
    ref = FIXTURES_DIR / "reference/sharded/@.bim"
    source = FIXTURES_DIR / ("genotype/sharded" if sharded else "genotype/nonsharded")
    copy = tmp_path / "source"
    shutil.copytree(source, copy)
    prefix = copy / ("@" if sharded else "all")
    override = source / ("@.bed" if sharded else "all.bed")
    missing = tmp_path / ("missing@.bed" if sharded else "missing.bed")
    bad = tmp_path / "bad.bed"
    bad.write_bytes(b"bad")
    bad_header = tmp_path / "bad_header.bed"
    sample_bed = source / ("X.bed" if sharded else "all.bed")
    payload = sample_bed.read_bytes()
    bad_header.write_bytes(b"\x00\x00\x01" + payload[3:])
    cache = tmp_path / ("cache.rds" if runtime == "r" else "cache.mat")
    saved = tmp_path / ("saved.rds" if runtime == "r" else "saved.mat")
    beds = [copy / "1.bed", copy / "X.bed"] if sharded else [copy / "all.bed"]
    if runtime == "r":
        q = lambda x: json.dumps(str(x))
        delete = "; ".join(f"unlink({q(p)})" for p in beds)
        expr = (
            f"ref <- load_reference({q(ref)}); g <- load_genotype({q(prefix)}, ref); "
            f"save_genotype_cache(g, {q(cache)}); {delete}; "
            f"nwarn <- 0L; h <- withCallingHandlers(load_genotype_cache({q(cache)}), "
            "warning=function(w) { nwarn <<- nwarn + 1L; invokeRestart('muffleWarning') }); "
            "stopifnot(nwarn == 1L, num_snp(h) == num_snp(g)); "
            "fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error=identity), 'error'); "
            f"stopifnot(fails(fetch_genotypes_int8(h, 1)), fails(load_genotype_cache({q(cache)}, bed_path={q(missing)})), "
            f"fails(load_genotype_cache({q(cache)}, shards='X', bed_path={q(bad)}))); "
            f"writeBin(as.raw(c(1,2,3)), {q(beds[0])}); nwarn <- 0L; "
            f"withCallingHandlers(load_genotype_cache({q(cache)}), warning=function(w) {{ nwarn <<- nwarn+1L; invokeRestart('muffleWarning') }}); stopifnot(nwarn == 1L); "
            f"header <- load_genotype_cache({q(cache)}, shards='X', bed_path={q(bad_header)}); stopifnot(fails(fetch_genotypes_int8(header,1))); "
            f"h <- load_genotype_cache({q(cache)}, bed_path={q(override)}); "
            f"save_genotype_cache(h, {q(saved)}); h <- load_genotype_cache({q(saved)}); "
            "stopifnot(identical(fetch_genotypes_int8(h,c(1,6)), matrix(2L,4,2))); "
            f"x <- load_genotype_cache({q(cache)}, shards='X', bed_path={q(override)}); "
            "stopifnot(identical(fetch_genotypes_int8(x,1), matrix(2L,4,1)));"
        )
        result = run_rscript(_source_phase3_script(expr))
    else:
        delete = " ".join(f"delete('{p}');" for p in beds)
        script = (
            f"ref = statgen.load_reference('{ref}'); g = statgen.load_genotype('{prefix}', ref); "
            f"g.save_cache('{cache}'); {delete} lastwarn(''); h = statgen.load_genotype_cache('{cache}'); "
            "[msg, ~] = lastwarn(); assert(~isempty(strfind(msg, 'Missing or incorrectly sized'))); "
            "assert(h.num_snp == g.num_snp); failed = false; try; h.fetch_genotypes_int8(1); catch; failed = true; end; assert(failed); "
            f"failed = false; try; statgen.load_genotype_cache('{cache}', [], '{missing}'); catch; failed = true; end; assert(failed); "
            f"failed = false; try; statgen.load_genotype_cache('{cache}', {{'X'}}, '{bad}'); catch; failed = true; end; assert(failed); "
            f"fid = fopen('{beds[0]}', 'wb'); fwrite(fid, [1 2 3]); fclose(fid); lastwarn(''); "
            f"h = statgen.load_genotype_cache('{cache}'); [msg, ~] = lastwarn(); assert(~isempty(msg)); "
            f"header = statgen.load_genotype_cache('{cache}', {{'X'}}, '{bad_header}'); failed = false; "
            "try; header.fetch_genotypes_int8(1); catch; failed = true; end; assert(failed); "
            f"h = statgen.load_genotype_cache('{cache}', [], '{override}'); h.save_cache('{saved}'); "
            f"h = statgen.load_genotype_cache('{saved}'); assert(isequal(h.fetch_genotypes_int8([1 6]), int8(2*ones(4,2)))); "
            f"x = statgen.load_genotype_cache('{cache}', {{'X'}}, '{override}'); assert(isequal(x.fetch_genotypes_int8(1), int8(2*ones(4,1))));"
        )
        result = run_octave(script)
    assert result.returncode == 0, result.stderr
