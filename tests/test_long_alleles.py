"""Independent vectors, historical caches and bounded source processing."""
import csv
import warnings

import numpy as np
import pytest

from statgen import _utils
from statgen.reference import load_reference, load_reference_cache, save_reference_cache
from statgen.sumstats import load_sumstats
from statgen.genotype import load_genotype
from tests.conftest import (FIXTURES_DIR, R_PACKAGE_DIR, run_octave, run_rscript,
                            skipif_no_octave, skipif_no_rscript, matlab_output_text)

ROOT = FIXTURES_DIR / 'long_alleles'
CAP = 1_000_000
VECTORS = list(csv.DictReader((ROOT / 'vectors.tsv').open(), delimiter='\t'))


def sequences():
    return [(v['motif'] * ((int(v['length']) + len(v['motif']) - 1) // len(v['motif'])))[:int(v['length'])]
            for v in VECTORS]


def test_shared_vectors_and_historical_cache():
    expected = np.array([int(v['expected']) for v in VECTORS], dtype=np.uint64)
    with warnings.catch_warnings(record=True) as caught:
        np.testing.assert_array_equal(_utils.allele_hash64(sequences()), expected)
    assert not caught
    for v in VECTORS:
        assert (v['historical'] == v['expected']) == (int(v['length']) <= 150)
    old = load_reference_cache(ROOT / 'legacy_python.npz')
    new = load_reference(ROOT / 'legacy.bim')
    assert old.shards[0].checksum == new.shards[0].checksum
    assert old.a1_hash64[0] == int(next(v['historical'] for v in VECTORS if v['length'] == '193000'))
    assert old.a1_hash64[0] != new.a1_hash64[0]
    assert old.validate_checksums()
    np.testing.assert_array_equal(old.a1, new.a1)


@pytest.mark.parametrize('budget', [97, 8192, 262144])
def test_bounded_batches_and_long_outlier(monkeypatch, budget):
    monkeypatch.setattr(_utils, '_ALLELE_HASH_BATCH_BYTES', budget)
    # Ordinary SNPs/short indels plus rare long alleles; compare an independent oracle.
    chosen = [v for v in VECTORS if int(v['length']) in (1, 2, 150, 151, 20000)]
    alleles = [(v['motif'] * (int(v['length']) + 1))[:int(v['length'])] for v in chosen]
    flat_sizes = []
    original = np.frombuffer
    def track(data, *args, **kwargs):
        flat_sizes.append(len(data))
        return original(data, *args, **kwargs)
    monkeypatch.setattr(np, 'frombuffer', track)
    for copies in (1, 15):
        hashes = _utils.allele_hash64(alleles * copies)
        assert list(map(int, hashes)) == [int(v['expected']) for v in chosen] * copies
        for length in (1, 2):  # SNP-only and short-indel controls.
            control = [v for v in VECTORS if int(v['length']) == length]
            inputs = [(v['motif'] * 2)[:length] for v in control] * copies
            assert list(map(int, _utils.allele_hash64(inputs))) == [int(v['expected']) for v in control] * copies
    assert max(flat_sizes) <= budget


def write_sources(tmp_path, length=CAP + 1):
    a1, a2 = 'A' * (length - 1) + 'C', 'A' * (length - 1) + 'G'
    bim = tmp_path / 'source.bim'
    bim.write_text(f'1\tlong\t0\t100\t{a1}\t{a2}\n')
    # Ignored row and missing-p row exercise warning locations before filtering.
    tsv = tmp_path / 'source.tsv'
    tsv.write_text('chr\tbp\ta1\ta2\tp\tz\tn\nY\t10\tA\tG\t0.1\t1\t100\n'
                   f'1\t90\t{a1}\tG\tNA\t1\t100\n1\t100\t{a1}\t{a2}\t0.1\t2\t100\n')
    (tmp_path / 'source.fam').write_text('F\tI\t0\t0\t1\t-9\n')
    (tmp_path / 'source.bed').write_bytes(bytes([0x6c, 0x1b, 0x01, 0]))
    return bim, tsv, a1, a2


def test_source_warning_ownership_preservation_and_collision_roundtrip(tmp_path):
    bim, tsv, a1, a2 = write_sources(tmp_path)
    with pytest.warns(RuntimeWarning, match='hash cap') as caught:
        ref = load_reference(bim)
    assert len(caught) == 1
    assert '2 allele(s)' in str(caught[0].message)
    assert 'row 1, column a1, length 1000001' in str(caught[0].message)
    assert ref.a1[0] == a1 and ref.a2[0] == a2
    assert ref.a1_hash64[0] == ref.a2_hash64[0]
    cache = tmp_path / 'ref.npz'
    save_reference_cache(ref, cache)
    loaded = load_reference_cache(cache)
    assert loaded.validate_checksums()
    assert loaded.a1[0] == a1 and loaded.a2[0] == a2
    with pytest.warns(RuntimeWarning, match='hash cap') as caught:
        stats = load_sumstats(tsv, loaded)
    assert len(caught) == 1 and '3 allele(s)' in str(caught[0].message)
    assert 'row 3, column a1' in str(caught[0].message)
    assert stats.is_present[0]
    with pytest.warns(RuntimeWarning, match='hash cap') as caught:
        genotype = load_genotype(tmp_path / 'source', loaded)
    assert len(caught) == 1 and genotype.is_present[0]
    # Equal capped keys across reference rows still trigger duplicate validation.
    bim.write_text(bim.read_text() + f'1\tduplicate\t0\t100\t{a1}\t{a2}\n')
    with pytest.warns(RuntimeWarning, match='hash cap'), pytest.raises(ValueError, match='duplicate'):
        load_reference(bim)


def test_ten_megabyte_valid_and_invalid_suffix(tmp_path):
    bim, _, a1, _ = write_sources(tmp_path, 10_000_000)
    with pytest.warns(RuntimeWarning, match='hash cap') as caught:
        ref = load_reference(bim)
    assert len(caught) == 1 and ref.a1[0] == a1
    assert ref.a1_hash64[0] == _utils.allele_hash64(['A' * CAP])[0]
    bim.write_text(f'1\tinvalid\t0\t100\t{"A" * CAP}N\tG\n')
    with pytest.warns(RuntimeWarning, match='hash cap'), pytest.raises(ValueError, match='uppercase DNA'):
        load_reference(bim)


def vector_files(tmp_path):
    alleles = tmp_path / 'alleles.txt'
    alleles.write_text('\n'.join(sequences()) + '\n')
    return alleles


@pytest.mark.octave
@skipif_no_octave
def test_octave_shared_vectors_legacy_and_silent_helpers(tmp_path):
    alleles = vector_files(tmp_path)
    output = tmp_path / 'hashes.txt'
    result = run_octave(
        f"fid=fopen('{alleles}'); a=textscan(fid,'%s'); fclose(fid); lastwarn(''); "
        "h=statgen.internal.allele_hash64(a{1}); [msg, id]=lastwarn(); assert(isempty(msg)); "
        f"fid=fopen('{output}','w'); fprintf(fid,'%u\\n',h); fclose(fid); "
        f"old=statgen.load_reference_cache('{ROOT / 'legacy_matlab.mat'}'); "
        f"new=statgen.load_reference('{ROOT / 'legacy.bim'}'); "
        "assert(old.a1_hash64 ~= new.a1_hash64); assert(strcmp(old.shards{1}.checksum,new.shards{1}.checksum)); "
        "assert(isequal(old.a1,new.a1));", timeout=120)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [v['expected'] for v in VECTORS]


@pytest.mark.r
@skipif_no_rscript
def test_r_shared_vectors_legacy_and_silent_helpers(tmp_path):
    alleles = vector_files(tmp_path)
    output = tmp_path / 'hashes.txt'
    result = run_rscript(
        f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
        'options(warn=2); '
        f'h <- .allele_hash64(readLines("{alleles}")); writeLines(as.character(h),"{output}"); '
        f'old <- load_reference_cache("{ROOT / "legacy_r.rds"}"); new <- load_reference("{ROOT / "legacy.bim"}"); '
        'stopifnot(old$shards[[1]]$checksum == new$shards[[1]]$checksum, identical(a1(old),a1(new)), '
        'a1_hash64(old) != a1_hash64(new));', timeout=120)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [v['expected'] for v in VECTORS]


def _with_cap_warnings(script):
    # MATLAB tests share an engine: restore warning state on success and error.
    saved = "statgen_test_cap_warning"
    return (
        f"{saved}=warning('query','statgen:allele_hash'); "
        "warning('on','statgen:allele_hash'); try; " + script + "; "
        f"catch ME; warning({saved}.state,'statgen:allele_hash'); rethrow(ME); end; "
        f"warning({saved}.state,'statgen:allele_hash');"
    )


@pytest.mark.octave
@skipif_no_octave
def test_octave_loader_warnings_collision_roundtrip_and_full_validation(tmp_path):
    bim, tsv, _, _ = write_sources(tmp_path)
    cache = tmp_path / 'ref.mat'
    result = run_octave(_with_cap_warnings(
        f"ref=statgen.load_reference('{bim}'); "
        "assert(ref.a1_hash64 == ref.a2_hash64); "
        f"statgen.save_reference_cache(ref,'{cache}','format','v7'); lastwarn(''); "
        f"loaded=statgen.load_reference_cache('{cache}'); "
        "[msg,id]=lastwarn(); assert(isempty(msg)); "
        "assert(isequal(ref.a1,loaded.a1) && isequal(ref.a2,loaded.a2)); "
        "assert(loaded.validate_checksums()); "
        f"s=statgen.load_sumstats('{tsv}',loaded); assert(s.is_present); "
        f"g=statgen.load_genotype('{tmp_path / 'source'}',loaded); assert(g.is_present);"), timeout=120)
    assert result.returncode == 0, result.stderr
    console = matlab_output_text(result.stdout + result.stderr)
    assert console.count('allele(s) exceed hash cap') == 3
    assert 'row 3, column a1, length 1000001' in console
    # Validate oversized strings in blocks, including the suffix beyond the cap.
    result = run_octave(
        "a=repmat('A',1,10000000); lastwarn(''); "
        "[bad,lens]=statgen.internal.invalid_dna_allele({a}); assert(~bad && lens==10000000); "
        "h=statgen.internal.allele_hash64({a}); assert(h==statgen.internal.allele_hash64({a(1:1000000)})); "
        "a(end)='N'; assert(statgen.internal.invalid_dna_allele({a})); "
        "[msg,id]=lastwarn(); assert(isempty(msg));", timeout=120)
    assert result.returncode == 0, result.stderr
    # The full source path must also accept a valid 10 MB sequence and reject invalid suffixes.
    bim, _, _, _ = write_sources(tmp_path, 10_000_000)
    result = run_octave(_with_cap_warnings(
        f"ref=statgen.load_reference('{bim}'); assert(numel(ref.a1{{1}})==10000000);"), timeout=120)
    assert result.returncode == 0, result.stderr
    assert matlab_output_text(result.stdout + result.stderr).count('allele(s) exceed hash cap') == 1
    bim.write_text(f'1\tinvalid\t0\t100\t{"A" * CAP}N\tG\n')
    result = run_octave(
        f"failed=false; try; statgen.load_reference('{bim}'); catch ME; "
        "failed=~isempty(strfind(ME.message,'uppercase DNA')); end; assert(failed);", timeout=120)
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_loader_warnings_collision_roundtrip_and_full_validation(tmp_path):
    bim, tsv, _, _ = write_sources(tmp_path)
    cache = tmp_path / 'ref.rds'
    # Isolate post-reader behavior: unchanged read.table is too slow on these fields.
    setup = (f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
             'assignInNamespace("read.table", function(file, ...) '
             'data.table::fread(file,header=FALSE,quote="",colClasses="character",'
             'na.strings=NULL,data.table=FALSE,nThread=1L,showProgress=FALSE),ns="utils"); ')
    result = run_rscript(setup +
        'messages <- character(); withCallingHandlers({ '
        f'ref <- load_reference("{bim}"); stopifnot(a1_hash64(ref)==a2_hash64(ref)); '
        f'save_reference_cache(ref,"{cache}"); loaded <- load_reference_cache("{cache}"); '
        'stopifnot(identical(a1(ref),a1(loaded)),identical(a2(ref),a2(loaded)),validate_checksums(loaded)); '
        f's <- load_sumstats("{tsv}",loaded); stopifnot(is_present(s)); '
        f'g <- load_genotype("{tmp_path / "source"}",loaded); stopifnot(is_present(g)); '
        '}, warning=function(w) { messages <<- c(messages,conditionMessage(w)); invokeRestart("muffleWarning") }); '
        'stopifnot(length(messages)==3,grepl("row 3, column a1, length 1000001",messages[[2]],fixed=TRUE));', timeout=120)
    assert result.returncode == 0, result.stderr
    bim, _, _, _ = write_sources(tmp_path, 10_000_000)
    result = run_rscript(setup +
        'messages <- character(); withCallingHandlers({ '
        f'ref <- load_reference("{bim}"); '
        'stopifnot(nchar(a1(ref))==10000000,a1_hash64(ref)==.allele_hash64(strrep("A",1000000))); '
        '}, warning=function(w) { messages <<- c(messages,conditionMessage(w)); invokeRestart("muffleWarning") }); '
        'stopifnot(length(messages)==1);', timeout=120)
    assert result.returncode == 0, result.stderr
    bim.write_text(f'1\tinvalid\t0\t100\t{"A" * CAP}N\tG\n')
    result = run_rscript(setup +
        f'failed <- tryCatch({{suppressWarnings(load_reference("{bim}")); FALSE}}, '
        'error=function(e) grepl("uppercase DNA",conditionMessage(e),fixed=TRUE)); stopifnot(failed);', timeout=120)
    assert result.returncode == 0, result.stderr


def test_differences_after_former_prefix_and_loader_cap_boundaries(tmp_path):
    assert _utils.allele_hash64(['A' * 150 + 'C'])[0] != _utils.allele_hash64(['A' * 150 + 'G'])[0]
    for length in (CAP - 1, CAP, CAP + 1):
        bim = tmp_path / 'boundary.bim'
        bim.write_text(f'1\tboundary\t0\t100\t{"A" * length}\tG\n')
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            load_reference(bim)
        assert len(caught) == (length > CAP)


@pytest.mark.octave
@skipif_no_octave
def test_octave_budgets_and_former_prefix_differences(tmp_path):
    (tmp_path / 'char.m').write_text(
        "function out=char(varargin)\n"
        "out=builtin('char',varargin{:});\n"
        "global STATGEN_TEST_MAX_CHARS;\n"
        "STATGEN_TEST_MAX_CHARS=max(STATGEN_TEST_MAX_CHARS,numel(out));\nend\n")
    (tmp_path / 'cellfun.m').write_text(
        "function out=cellfun(varargin)\n"
        "out=builtin('cellfun',varargin{:});\n"
        "global STATGEN_TEST_MAX_PIECES;\n"
        "if iscell(out); sizes=builtin('cellfun','length',out); "
        "STATGEN_TEST_MAX_PIECES=max(STATGEN_TEST_MAX_PIECES,sum(sizes(:))); end\nend\n")
    expected_cap = int(next(v['expected'] for v in VECTORS if v['length'] == str(CAP)))
    result = run_octave(
        f"oldpath=path(); restore=onCleanup(@() path(oldpath)); try; addpath('{tmp_path}'); "
        "global STATGEN_TEST_MAX_CHARS STATGEN_TEST_MAX_PIECES; STATGEN_TEST_MAX_PIECES=0; "
        # SNV fast path, short alleles spanning batches, and a long outlier.
        "cases={repmat({'A';'C';'G';'T'},250,1), repmat({'A';'AC';'ACGT';'G'},250,1), "
        "{'A';'AC';repmat('ACGT',1,5000)}}; budgets=[65536 4096 4096]; "
        "for scenario=1:numel(cases); a=cases{scenario}; budget=budgets(scenario); "
        "lens=cellfun('length',a); STATGEN_TEST_MAX_CHARS=0; "
        "assert(~any(statgen.internal.invalid_dna_allele(a,lens,budget))); assert(STATGEN_TEST_MAX_CHARS<=budget/32); "
        "STATGEN_TEST_MAX_CHARS=0; h=statgen.internal.allele_hash64(a,budget); assert(STATGEN_TEST_MAX_CHARS<=budget/64); "
        "expected=statgen.internal.allele_hash64(a); assert(isequal(h,expected)); "
        "a{end}(end)='N'; bad=statgen.internal.invalid_dna_allele(a,lens,budget); assert(sum(bad)==1 && bad(end)); "
        "end; p=repmat('A',1,150); h=statgen.internal.allele_hash64({[p 'C'];[p 'G']}); assert(h(1)~=h(2)); "
        "seq=repmat('ACGTGCA',1,ceil(1000000/7)); seq=[seq(1:1000000) 'C']; "
        f"expected_cap=bitshift(uint64({expected_cap >> 32}),32)+uint64({expected_cap & 0xffffffff}); "
        "for budget=[65536 16*1024*1024]; count=1; if budget>65536; count=32; end; a=repmat({seq},count,1); "
        "STATGEN_TEST_MAX_CHARS=0; STATGEN_TEST_MAX_PIECES=0; h=statgen.internal.allele_hash64(a,budget); "
        "assert(all(h==expected_cap)); assert(STATGEN_TEST_MAX_CHARS<=budget/64); assert(STATGEN_TEST_MAX_PIECES<=budget/64); "
        "end; catch ME; clear restore; rethrow(ME); end; clear restore;", timeout=120)
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_budgets_and_former_prefix_differences():
    result = run_rscript(
        f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
        'seen <- 0; calls <- 0L; trace("charToRaw",where=baseenv(),print=FALSE, '
        'tracer=quote({ .GlobalEnv$seen <- max(.GlobalEnv$seen,nchar(x,type="bytes")); .GlobalEnv$calls <- .GlobalEnv$calls+1L })); '
        # SNV fast path, bounded short batches, and an isolated long block.
        'cases <- list(rep(c("A","C","G","T"),1000),rep(c("A","AC","ACGT","G"),1000), '
        'c("A","AC",strrep("ACGT",5000))); budgets <- c(524288L,97L,1024L); '
        'for(i in seq_along(cases)) { a <- cases[[i]]; expected <- .allele_hash64(a); '
        'budget <- budgets[[i]]; .allele_hash_batch_bytes <- budget; '
        'seen <- 0; stopifnot(identical(expected,.allele_hash64(a)),seen<=2*budget); }; '
        '.allele_hash_batch_bytes <- 8192L; calls <- 0L; h <- .allele_hash64(rep("AC",1000)); stopifnot(calls==1L); '
        'untrace("charToRaw",where=baseenv()); .warn_allele_cap(c(NA_character_,"A"),c("G","T"),"test",1:2); '
        'p <- strrep("A",150); h <- .allele_hash64(c(paste0(p,"C"),paste0(p,"G"))); stopifnot(h[[1]]!=h[[2]]);', timeout=120)
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_batch_offsets_beyond_integer_range():
    result = run_rscript(
        f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
        'options(warn=2); .allele_hash_batch_bytes <- 8192L; lengths <- c(rep.int(1000000L,2148L),rep.int(2L,20000L)); '
        'offsets <- cumsum(as.double(lengths)); stopifnot(tail(offsets,1)>.Machine$integer.max); '
        'b <- .allele_hash_batches(lengths); '
        'stopifnot(!anyNA(b$starts),!anyNA(b$ends),b$starts[[1]]==1L,tail(b$ends,1)==length(lengths), '
        'identical(b$starts[-1L],head(b$ends,-1L)+1L)); '
        'bytes <- offsets[b$ends]-c(0,offsets[head(b$ends,-1L)]); '
        'large <- lengths[b$starts]>.allele_hash_batch_bytes; '
        'stopifnot(sum(large)==2148L,all(b$starts[large]==b$ends[large]), '
        'all(bytes[!large]<=2*.allele_hash_batch_bytes),sum(!large)>1L);')
    assert result.returncode == 0, result.stderr


@pytest.mark.octave
@skipif_no_octave
def test_octave_cap_warning_settings_restored():
    success = _with_cap_warnings("warning('statgen:allele_hash','cap warning probe')")
    failure = _with_cap_warnings("error('statgen:warning_test','expected test failure')")
    result = run_octave(
        "original_cap_state=warning('query','statgen:allele_hash'); warning('off','statgen:allele_hash'); try; "
        + success + " state=warning('query','statgen:allele_hash'); assert(strcmp(state.state,'off')); "
        "failed=false; try; " + failure + " catch ME; failed=strcmp(ME.identifier,'statgen:warning_test'); end; assert(failed); "
        "state=warning('query','statgen:allele_hash'); assert(strcmp(state.state,'off')); "
        "catch ME; warning(original_cap_state.state,'statgen:allele_hash'); rethrow(ME); end; "
        "warning(original_cap_state.state,'statgen:allele_hash');")
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_batch_memory_and_exactness_limits_are_independent():
    result = run_rscript(
        f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
        'options(warn=2); .allele_hash_batch_bytes <- 524288L; '
        'b <- .allele_hash_batches(c(8192L,8193L,8192L,2L)); '
        'stopifnot(identical(b$starts,c(1L,2L,3L)),identical(b$ends,c(1L,2L,4L))); '
        'lengths <- rep.int(8192L,129L); b <- .allele_hash_batches(lengths); '
        'stopifnot(all((b$ends-b$starts+1L)*8192<=2*.allele_hash_batch_bytes)); '
        'writes <- 0L; trace("[<-.integer64",where=asNamespace("bit64"),print=FALSE, '
        'tracer=quote(.GlobalEnv$writes <- .GlobalEnv$writes+1L)); '
        # Multiple batches at small and large budgets, including both sides of
        # the exactness limit. Every lane write must use plain numeric vectors.
        'a <- rep(c("A","AC","G",strrep("ACGT",2048L),paste0(strrep("ACGT",2048L),"A")),20L); '
        'expected <- .allele_hash64(a); for(budget in c(1024L,524288L)) { '
        '.allele_hash_batch_bytes <- budget; stopifnot(identical(.allele_hash64(a),expected)); }; '
        'stopifnot(writes==0L); untrace("[<-.integer64",where=asNamespace("bit64"));', timeout=60)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('length', [0, 1, 2, 3, 4, 7, 8, 9, 31, 32, 33, 8191, 8192, 8193, 193000, CAP])
def test_python_power_doubling_matches_scalar_recurrence(length):
    for base in (257, 263):
        def scalar():
            value = 1
            for _ in range(length + 1):
                yield value
                value = value * base % 2147483647
        expected = np.fromiter(scalar(), dtype=np.int64, count=length + 1)
        actual = _utils._allele_hash_powers(base, length)
        np.testing.assert_array_equal(actual, expected)
        assert actual.dtype == np.int64


@pytest.mark.r
@skipif_no_rscript
def test_r_power_doubling_is_exact(tmp_path):
    points = sorted({0, 1, 2, CAP} | {2**k + delta for k in range(1, 20) for delta in (-1, 0, 1)})
    path = tmp_path / 'powers.tsv'
    path.write_text('exponent\th1\th2\n' + ''.join(
        f'{n}\t{pow(257,n,2147483647)}\t{pow(263,n,2147483647)}\n' for n in points))
    result = run_rscript(
        f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
        f'expected <- read.table("{path}",header=TRUE); options(warn=2); '
        'for(n in c(0L,1L,2L,3L,7L,8L,9L,1023L,1024L,1025L,1000000L)) { '
        'actual <- .allele_hash_powers(n); selected <- expected[expected$exponent<=n,]; '
        'stopifnot(identical(dim(actual),c(n+1L,2L)),typeof(actual)=="double", '
        'all(actual[selected$exponent+1L,1L]==selected$h1),all(actual[selected$exponent+1L,2L]==selected$h2)); }')
    assert result.returncode == 0, result.stderr


def write_filtered_sources(tmp_path):
    bim = tmp_path / 'reuse.bim'
    bim.write_text(f'1\tlong\t0\t100\t{"A" * 151}\tGT\n2\tshort\t0\t200\tACG\tT\n')
    (tmp_path / 'reuse.bed').write_bytes(bytes([0x6c,0x1b,0x01,0,0]))
    (tmp_path / 'reuse.fam').write_text('F\tI\t0\t0\t1\t-9\n')
    stats = tmp_path / 'reuse.tsv'
    stats.write_text('chr\tbp\ta1\ta2\tp\tz\tn\nY\t10\tA\tG\t0.1\t1\t100\n'
                     f'1\t90\t{"C" * 17}\tT\tNA\t1\t100\n2\t200\tACG\tT\t0.1\t3\t100\n'
                     f'1\t100\t{"A" * 151}\tGT\t0.1\t2\t100\n')
    return bim, stats


def test_python_source_loaders_hash_filtered_alleles(tmp_path):
    bim, stats = write_filtered_sources(tmp_path)
    ref = load_reference(bim)
    assert load_reference(bim, shards=['2']).num_snp == 1
    assert load_genotype(tmp_path / 'reuse', ref).is_present.all()
    np.testing.assert_array_equal(load_sumstats(stats, ref).zvec, [2,3])


@pytest.mark.r
@skipif_no_rscript
def test_r_source_loaders_hash_filtered_alleles(tmp_path):
    bim, stats = write_filtered_sources(tmp_path)
    result = run_rscript(
        f'for(f in list.files("{R_PACKAGE_DIR / "R"}",pattern="[.]R$",full.names=TRUE)) source(f); '
        f'ref <- load_reference("{bim}"); stopifnot(num_snp(load_reference("{bim}",shards="2"))==1L); '
        f'stopifnot(all(is_present(load_genotype("{tmp_path / "reuse"}",ref)))); '
        f'stopifnot(identical(zvec(load_sumstats("{stats}",ref)),c(2,3)));', timeout=60)
    assert result.returncode == 0, result.stderr
