"""Stage 3: R against the independent portable allele-swap oracle."""
import json
import shutil

import numpy as np
import pytest

from statgen.genotype import load_genotype
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, R_PACKAGE_DIR, run_rscript, skipif_no_rscript, write_plink_bed_calls

pytestmark = [pytest.mark.r, skipif_no_rscript]
ROOT = FIXTURES_DIR / 'allele_swaps'
EXPECTED = json.loads((ROOT / 'expected.json').read_text())
CACHE = FIXTURES_DIR / 'cache/genotype_r_genotype_cache_0_2.rds'


def _run(script):
    result = run_rscript(
        f'for (f in list.files({json.dumps(str(R_PACKAGE_DIR / "R"))}, pattern="[.]R$", full.names=TRUE)) source(f); '
        + script, timeout=60,
    )
    assert result.returncode == 0, f'{result.stdout}\n{result.stderr}'
    return result


def _fails(expr, pattern):
    return f'failed <- FALSE; tryCatch({{{expr}}}, error=function(e) failed <<- grepl({json.dumps(pattern)}, conditionMessage(e), fixed=TRUE)); stopifnot(failed); '


@pytest.mark.parametrize('enabled', [False, True])
def test_shared_oracle_roundtrips_and_python_parity(enabled, tmp_path):
    expected = EXPECTED['options'][str(enabled).lower()]
    requested = np.r_[np.flatnonzero(np.array(expected['source_row0']) >= 0)[::-1], [1, 13, 1]]
    indices = ','.join(str(i + 1) + 'L' for i in requested)
    output = tmp_path / 'aligned.json'
    cache, scache = tmp_path / 'genotype.rds', tmp_path / 'sumstats.rds'
    result = _run(
        f'ref <- load_reference("{ROOT / "reference.bim"}"); '
        f'g <- load_genotype("{ROOT / "genotype/@"}", ref, allow_indel_swaps={str(enabled).upper()}); '
        f's <- load_sumstats("{ROOT / "sumstats.tsv.gz"}", ref, allow_indel_swaps={str(enabled).upper()}); '
        f'indices <- c({indices}); raw <- fetch_genotypes_int8(g, indices); '
        'stopifnot(is.logical(is_allele_swapped(g)), is.null(dim(is_allele_swapped(g))), typeof(raw) == "integer", '
        'identical(raw, fetch_genotypes_int8(g, indices))); '
        f'save_cache(g, "{cache}"); h <- load_genotype_cache("{cache}"); '
        'stopifnot(identical(is_allele_swapped(h), is_allele_swapped(g)), identical(source_row0(h), source_row0(g)), '
        'identical(fetch_genotypes_int8(h, indices), raw), '
        'identical(fetch_genotypes(h, indices, haploid_mode="ploidy_scaled"), fetch_genotypes(g, indices, haploid_mode="ploidy_scaled"))); '
        f'gx <- select_shards(g, "X"); hx <- load_genotype_cache("{cache}", "X"); '
        'stopifnot(identical(is_allele_swapped(gx), is_allele_swapped(g)[14:16]), '
        'identical(is_allele_swapped(hx), is_allele_swapped(gx)), '
        'identical(fetch_genotypes_int8(hx, c(1L,2L,1L)), fetch_genotypes_int8(gx, c(1L,2L,1L)))); '
        + _fails('fetch_genotypes_int8(g, 9L)', 'not present') +
        f'save_cache(s, "{scache}"); sc <- load_sumstats_cache("{scache}"); '
        'stopifnot(identical(zvec(sc), zvec(s)), identical(eaf_vec(sc), eaf_vec(s))); '
        'out <- list(source_row0=source_row0(g), swapped=is_allele_swapped(g), present=is_present(g), raw=raw, '
        'floating=fetch_genotypes(g, indices), scaled=fetch_genotypes(g, indices, haploid_mode="ploidy_scaled"), '
        'logp=logpvec(s), z=zvec(s), n=nvec(s), beta=beta_vec(s), se=se_vec(s), eaf=eaf_vec(s), info=info_vec(s)); '
        f'jsonlite::write_json(out, "{output}", digits=NA, na="null");'
    )
    data = json.loads(output.read_text())
    np.testing.assert_array_equal(data['source_row0'], expected['source_row0'])
    np.testing.assert_array_equal(data['swapped'], expected['is_allele_swapped'])
    np.testing.assert_array_equal(data['present'], np.array(expected['source_row0']) >= 0)
    raw = np.array(expected['raw_counts'], dtype=float)[:, requested]
    np.testing.assert_array_equal(data['raw'], np.where(np.isnan(raw), -1, raw).astype(int))
    np.testing.assert_allclose(np.array(data['floating'], dtype=float), raw, equal_nan=True)
    np.testing.assert_allclose(np.array(data['scaled'], dtype=float), np.array(expected['scaled_counts'], dtype=float)[:, requested], equal_nan=True)
    for field in ('z', 'n', 'beta', 'se', 'eaf', 'info'):
        np.testing.assert_allclose(np.array(data[field], dtype=float), np.array(expected['aligned_sumstats'][field], dtype=float), equal_nan=True)
    np.testing.assert_allclose(np.array(data['logp'], dtype=float), -np.log10(np.array(expected['aligned_sumstats']['p'], dtype=float)), equal_nan=True)
    python = load_genotype(ROOT / 'genotype/@', load_reference(ROOT / 'reference.bim'), allow_indel_swaps=enabled)
    np.testing.assert_array_equal(data['raw'], python.fetch_genotypes_int8(requested))
    for label, counts in expected['reporting'].items():
        for path in (ROOT / 'sumstats.tsv.gz', ROOT / f'genotype/{label}.bim'):
            prefix = f'{path}: shard {label}: '
            for key, message in [('swaps', 'successful allele swap(s)'), ('strand_ambiguous', 'strand-ambiguous A/T or C/G swap(s)'), ('disabled_indels', 'indel row(s) unmatched because indel swaps are disabled')]:
                assert (prefix + f'{counts[key]} {message}' in result.stderr) == bool(counts[key])
    assert '0 strand-ambiguous' not in result.stderr
    assert '0 indel row(s)' not in result.stderr
    assert 'remain unmatched' not in result.stderr


def test_matcher_ambiguity_precedence_counts_and_hash_precision():
    _run(
        'h <- .allele_hash64(c("A","G","AG","AC","GT")); '
        'm <- function(bp,a,b,enabled) .match_variant_keys(100L,h[1],h[2],bp,a,b,"1","sumstats",enabled); '
        'for (enabled in c(FALSE,TRUE)) { for (swapped in c(FALSE,TRUE)) { '
        'a <- h[1L+swapped]; b <- h[2L-swapped]; '
        + _fails('m(c(100L,100L),c(a,a),c(b,b),enabled)', 'Ambiguous duplicate') +
        '}; x <- m(c(100L,100L,100L,999L,999L),h[c(2,1,2,1,1)],h[c(1,2,1,2,2)],enabled); '
        'stopifnot(identical(x$source_rows,2L), !x$is_allele_swapped, x$num_disabled_indels==0L); }; '
        'x <- .match_variant_keys(c(100L,200L,200L,300L,400L),h[c(1,1,3,1,4)],h[c(3,3,1,2,5)], '
        'c(100L,100L,200L,300L,400L,999L),h[c(3,3,1,2,5,3)],h[c(1,1,3,1,4,1)],"1","sumstats",FALSE); '
        'stopifnot(identical(x$source_rows,c(NA_integer_,3L,NA_integer_,4L,NA_integer_)), x$num_disabled_indels==2L); '
        + _fails('.match_variant_keys(100L,h[1],h[3],c(100L,100L),h[c(3,3)],h[c(1,1)],"1","sumstats",TRUE)', 'Ambiguous duplicate') +
        'for (value in c("1152921504606846977","9218868437227405313")) { '
        'a <- bit64::as.integer64(value); b <- a + bit64::as.integer64(1); '
        'stopifnot(suppressWarnings(as.double(a)==as.double(b))); '
        'x <- .match_variant_keys(c(100L,100L),c(a,b),c(b,a),c(100L,100L),c(b,a),c(a,b),"1","sumstats",TRUE); '
        'stopifnot(identical(x$source_rows,c(2L,1L)), !any(x$is_allele_swapped)); }; '
        'high <- .allele_hash64(c("GTAGAA","GTAGAC")); single <- h[c(1,1)]; '
        'stopifnot(all(high >= bit64::as.integer64("9218868437227405312"))); '
        'for (enabled in c(FALSE,TRUE)) { '
        'x <- .match_variant_keys(c(100L,100L),single,high,c(100L,100L),rev(high),single,"1","sumstats",enabled); '
        'if (enabled) { stopifnot(identical(x$source_rows,c(2L,1L)), all(x$is_allele_swapped), x$num_disabled_indels==0L); '
        '} else { stopifnot(all(is.na(x$source_rows)), !any(x$is_allele_swapped), x$num_disabled_indels==2L); }; }; '
        + _fails('.match_variant_keys(c(100L,100L),single,high,rep(100L,3),high[c(1,2,1)],h[c(1,1,1)],"1","sumstats",TRUE)', 'Ambiguous duplicate') +
        'x <- .match_variant_keys(integer(),h[integer()],h[integer()],100L,h[1],h[2],"1","sumstats"); '
        'stopifnot(length(x$source_rows)==0L, length(x$is_allele_swapped)==0L);'
    )


@pytest.mark.parametrize('alleles', [('A','G'),('A','T'),('C','G')])
def test_reference_rejects_reciprocal_snvs(alleles, tmp_path):
    a, b = alleles
    path = tmp_path / 'reference.bim'
    path.write_text('1\tlong\t0\t50\tAC\tGT\n1\tsnv\t0\t75\tA\tC\n' + f'1\ta\t0\t100\t{a}\t{b}\n1\tb\t0\t100\t{b}\t{a}\n')
    _run(_fails(f'load_reference("{path}")', 'reference shard 1: reciprocal SNV keys are not allowed at bp=100'))


def test_cached_reference_is_trusted_without_string_decoding(tmp_path):
    path, source, cache = tmp_path / 'ref.bim', tmp_path / 'source.tsv', tmp_path / 'ref.rds'
    path.write_text('1\ta\t0\t100\tA\tG\n1\tb\t0\t100\tC\tT\n')
    source.write_text('chr\tbp\ta1\ta2\tp\tz\tn\n1\t100\tA\tG\t0.1\t1\t100\n1\t100\tG\tA\t0.2\t2\t100\n')
    _run(f'ref <- load_reference("{path}"); save_cache(ref,"{cache}"); p <- readRDS("{cache}"); '
         'p$a1_hash64[2] <- p$a2_hash64[1]; p$a2_hash64[2] <- p$a1_hash64[1]; '
         f'p$a1_text_by_shard <- p$a2_text_by_shard <- "invalid"; saveRDS(p,"{cache}"); '
         f'ref <- load_reference_cache("{cache}"); s <- load_sumstats("{source}",ref); stopifnot(identical(zvec(s), c(1,2)));')


def test_quiet_and_optional_values(tmp_path):
    source = tmp_path / 'source.tsv'
    source.write_text('chr\tbp\ta1\ta2\tp\n1\t200\tG\tA\t0\n1\t100\tA\tG\t0.1\n1\t100\tG\tA\t0.2\n1\t100\tG\tA\t0.3\n1\t9999\tA\tC\t0.1\n1\t9999\tA\tC\t0.1\n')
    result = _run(f'set_verbosity("quiet"); ref <- load_reference("{ROOT / "reference.bim"}"); '
                  f's <- load_sumstats("{source}",ref); stopifnot(is_present(s)[2], is.infinite(logpvec(s)[2]), is.null(zvec(s)), is.null(eaf_vec(s))); '
                  f'g <- load_genotype("{ROOT / "genotype/@"}",ref);')
    assert 'successful allele swap' not in result.stderr
    assert 'indel row(s)' not in result.stderr
    assert 'strand-ambiguous' not in result.stderr
    assert 'zvec' in result.stderr and 'nvec' in result.stderr
    source.write_text('chr\tbp\ta1\ta2\tp\tz\tn\tbeta\tse\teaf\tinfo\n1\t200\tG\tA\t0.1\tNA\tNA\tNA\tNA\tNA\tNA\n')
    _run(f'ref <- load_reference("{ROOT / "reference.bim"}", "1"); s <- load_sumstats("{source}",ref); '
         'stopifnot(is_present(s)[2], logpvec(s)[2]==1, all(is.nan(c(zvec(s)[2],nvec(s)[2],beta_vec(s)[2],se_vec(s)[2],eaf_vec(s)[2],info_vec(s)[2]))));')


def test_cache_fixture_legacy_migration_and_swapped_overrides(tmp_path):
    source = tmp_path / 'source'
    shutil.copytree(ROOT / 'genotype', source)
    cache, migrated = tmp_path / 'cache.rds', tmp_path / 'migrated.rds'
    _run(f'ref <- load_reference("{ROOT / "reference.bim"}"); g <- load_genotype("{source / "@"}",ref,TRUE); '
         f'fixture <- load_genotype_cache("{CACHE}"); '
         'stopifnot(identical(is_allele_swapped(fixture),is_allele_swapped(g)), identical(source_row0(fixture),source_row0(g))); '
         'indices <- c(which(is_present(g)),2L,5L,14L,16L,2L); expected <- fetch_genotypes_int8(g,indices); '
         'stopifnot(identical(fetch_genotypes_int8(fixture,indices),expected)); '
         f'save_cache(g,"{cache}"); unlink(list.files("{source}",full.names=TRUE)); h <- load_genotype_cache("{cache}"); '
         f'stopifnot(identical(fetch_genotypes_int8(h,indices,"{ROOT / "genotype/@.bed"}"),expected), '
         f'identical(fetch_genotypes(h,indices,"{ROOT / "genotype/@.bed"}","ploidy_scaled"),fetch_genotypes(fixture,indices,haploid_mode="ploidy_scaled"))); '
         f'old <- load_genotype_cache("{FIXTURES_DIR / "cache/genotype_r_genotype_cache_0_1.rds"}"); '
         'stopifnot(!any(is_allele_swapped(old)), identical(fetch_genotypes_int8(old,c(1L,6L,1L)),matrix(2L,4,3))); '
         f'save_cache(old,"{migrated}"); m <- load_genotype_cache("{migrated}"); '
         'stopifnot(!any(is_allele_swapped(m)), identical(source_row0(m),source_row0(old)));')


@pytest.mark.parametrize('mutation', [
    'p$is_allele_swapped <- NULL', 'p$is_allele_swapped <- p$is_allele_swapped[-1L]',
    'dim(p$is_allele_swapped) <- c(4L,4L)', 'p$is_allele_swapped <- as.numeric(p$is_allele_swapped)',
    'p$is_allele_swapped[9L] <- TRUE', 'p$is_allele_swapped[1L] <- NA',
])
def test_cache_rejects_invalid_swap_flags(mutation, tmp_path):
    path = tmp_path / 'bad.rds'
    _run(f'p <- readRDS("{CACHE}"); {mutation}; saveRDS(p,"{path}"); ' + _fails(f'load_genotype_cache("{path}")', 'is_allele_swapped'))


def test_swap_state_independent_of_missing_hardcalls(tmp_path):
    for suffix in ('.bim','.fam','.bed'):
        shutil.copyfile((ROOT / 'genotype/1').with_suffix(suffix), (tmp_path / 'source').with_suffix(suffix))
    bed = tmp_path / 'source.bed'
    payload = bytearray(bed.read_bytes())
    payload[5] = 0x55
    bed.write_bytes(payload)
    _run(f'ref <- load_reference("{ROOT / "reference.bim"}","1"); g <- load_genotype("{tmp_path / "source"}",ref); '
         'stopifnot(is_present(g)[2],is_allele_swapped(g)[2],identical(fetch_genotypes_int8(g,c(2L,2L)),matrix(-1L,4,2)));')


@pytest.mark.parametrize('value', ['1', 'c(TRUE,FALSE)', 'NA', 'NULL', '"true"'])
def test_option_requires_logical_scalar(value):
    _run(f'ref <- load_reference("{ROOT / "reference.bim"}"); '
         + _fails(f'load_genotype("{ROOT / "genotype/@"}",ref,{value})', 'allow_indel_swaps')
         + _fails(f'load_sumstats("{ROOT / "sumstats.tsv.gz"}",ref,{value})', 'allow_indel_swaps'))


def test_physical_rows_and_ploidy_survive_reference_projection(tmp_path):
    prefix = tmp_path / 'combined'
    prefix.with_suffix('.bim').write_text('1\tignored\t0\t1\tA\tC\nX\tsnv\t0\t100\tG\tA\nX\texact\t0\t200\tC\tT\nX\tindel\t0\t300\tAT\tA\n')
    shutil.copyfile(ROOT / 'genotype/1.fam', prefix.with_suffix('.fam'))
    prefix.with_suffix('.ploidy').write_text('2\t2\n1\t2\n2\t2\n0\t1\n')
    write_plink_bed_calls(prefix.with_suffix('.bed'), [[2,2,2,2],[0,1,2,-1],[2,1,0,-1],[0,1,2,-1]])
    _run(f'ref <- load_reference("{ROOT / "reference.bim"}","X"); g <- load_genotype("{prefix}",ref,TRUE); '
         'stopifnot(identical(source_row0(g),1:3),identical(is_allele_swapped(g),c(TRUE,FALSE,TRUE)), '
         'identical(ploidy_male(g),c(1,2,0)),identical(ploidy_female(g),c(2,2,1)), '
         'identical(fetch_genotypes_int8(g,c(1L,3L,1L)),matrix(rep(c(2L,1L,0L,-1L),3),4,3)), '
         'identical(fetch_genotypes(g,c(1L,3L),haploid_mode="ploidy_scaled"),matrix(c(1,1,0,NaN,0,0.5,0,NaN),4,2)));')


@pytest.mark.parametrize('verbosity', ['info','quiet'])
@pytest.mark.parametrize('row,warn', [('X\t100\tA\tG',True),('1\t9999\tA\tC',True),('1\t500\tAG\tA',True),('1\t200\tG\tA',False)])
def test_no_matching_sumstats_warning(row, warn, verbosity, tmp_path):
    path = tmp_path / 'source.tsv'
    path.write_text('chr\tbp\ta1\ta2\tp\tz\tn\n' + row + '\t0.1\t1\t100\n')
    result = _run(f'set_verbosity("{verbosity}"); ref <- load_reference("{ROOT / "reference.bim"}","1"); '
                  f's <- load_sumstats("{path}",ref); stopifnot(any(is_present(s))=={str(not warn).upper()});')
    assert (f"{path}: no matching source rows for reference shard '1'" in result.stderr) == warn


def test_shard_swap_validation_names_shard():
    _run(f'ref <- load_reference("{ROOT / "reference.bim"}"); '
         f'g <- load_genotype("{ROOT / "genotype/@"}",ref); '
         's <- shards(g)[[1]]; args <- s[names(formals(.new_genotype_shard))]; '
         'args$is_allele_swapped <- as.integer(s$is_allele_swapped); '
         + _fails('do.call(.new_genotype_shard,args)', 'GenotypeShard 1: is_allele_swapped must be logical and requires is_present') +
         'args$is_allele_swapped <- s$is_allele_swapped; args$is_allele_swapped[9] <- TRUE; '
         + _fails('do.call(.new_genotype_shard,args)', 'GenotypeShard 1: is_allele_swapped must be logical and requires is_present'))
