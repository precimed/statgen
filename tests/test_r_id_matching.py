"""R source-loader contract for exact ID and allele-group matching."""
import gzip
import json

import pytest

from tests.conftest import skipif_no_rscript
from tests.test_id_matching import bfile, reference, stats_file
from tests.test_r_allele_swaps import _fails as fails, _run as run

pytestmark = [pytest.mark.r, skipif_no_rscript]


def load_ref(tmp_path, rows, cached=True):
    reference(tmp_path, rows)
    script = f'ref <- load_reference("{tmp_path / "reference.bim"}"); '
    if cached:
        script += (f'save_cache(ref, "{tmp_path / "reference.rds"}"); '
                   f'ref <- load_reference_cache("{tmp_path / "reference.rds"}"); ')
    return script


@pytest.mark.parametrize('cached', [False, True])
@pytest.mark.parametrize('pair,pairs', [
    (('A', 'G'), [('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')]),
    (('A', 'C'), [('A', 'C'), ('C', 'A'), ('T', 'G'), ('G', 'T')]),
])
def test_orientations_and_cache_roundtrips(tmp_path, cached, pair, pairs):
    script = load_ref(tmp_path, [('1', f'id{i}', i + 1, *pair) for i in range(4)], cached)
    path = stats_file(tmp_path, [(f'id{i}', *p) for i, p in enumerate(pairs)],
                      beta=[2.] * 4, eaf=[.2] * 4, se=[.5] * 4, info=[.9] * 4)
    prefix = bfile(tmp_path, [('1', f'id{i}', 100 + i, *p) for i, p in enumerate(pairs)], calls=[0, 0, 2, -1])
    result = run(script +
        f's <- load_sumstats("{path}", ref, match_by="id"); '
        f'g <- load_genotype("{prefix}", ref, match_by="id"); '
        'stopifnot(identical(zvec(s), c(3,-3,3,-3)), identical(beta_vec(s), c(2,-2,2,-2)), '
        'max(abs(eaf_vec(s) - c(.2,.8,.2,.8))) < 1e-12, all(logpvec(s)==2), all(nvec(s)==100), '
        'all(se_vec(s)==.5), all(info_vec(s)==.9), identical(is_allele_swapped(g), c(FALSE,TRUE,FALSE,TRUE)), '
        'identical(fetch_genotypes_int8(g, c(1,2,3,4,2)), matrix(c(0L,2L,2L,-1L,2L),1))); '
        f'save_cache(g,"{tmp_path / "g.rds"}"); h <- load_genotype_cache("{tmp_path / "g.rds"}"); '
        'stopifnot(identical(source_row0(h), source_row0(g)), identical(is_allele_swapped(h), is_allele_swapped(g)), '
        'identical(fetch_genotypes_int8(h,c(4,2,1)), matrix(c(-1L,2L,0L),1))); '
        f'save_cache(s,"{tmp_path / "s.rds"}"); t <- load_sumstats_cache("{tmp_path / "s.rds"}"); '
        'stopifnot(identical(zvec(t),zvec(s)), identical(eaf_vec(t),eaf_vec(s)), '
        'is_object_compatible(ref,g), is_object_compatible(ref,s)); ')
    assert result.stderr.count('2 matched row(s) requiring inversion') == 2
    assert result.stderr.count('2 matched row(s) requiring strand complementation') == 2


def test_same_id_distinct_groups_and_source_order(tmp_path):
    script = load_ref(tmp_path, [('1', 'same', 1, 'A', 'G'), ('2', 'same', 2, 'A', 'C')])
    path = stats_file(tmp_path, [('same', 'G', 'T'), ('same', 'T', 'C')], z=[5, 7])
    prefix = bfile(tmp_path, [('1', 'same', 10, 'T', 'C'), ('2', 'same', 20, 'G', 'T')])
    run(script + f's <- load_sumstats("{path}",ref,match_by="id"); '
        f'g <- load_genotype("{prefix}",ref,match_by="id"); '
        'stopifnot(identical(zvec(s),c(7,-5)), identical(source_row0(g),c(0L,1L)));')


@pytest.mark.parametrize('sharded', [False, True])
def test_panel_wide_reference_duplicates_and_subset(tmp_path, sharded):
    script = load_ref(tmp_path, [('1', 'dup', 1, 'A', 'G'), ('2', 'dup', 2, 'T', 'C')])
    path = stats_file(tmp_path, [('dup', 'A', 'G'), ('dup', 'G', 'A')])
    if sharded:
        bfile(tmp_path, [('1', 'dup', 10, 'A', 'G')], name='source1')
        bfile(tmp_path, [('2', 'dup', 20, 'T', 'C')], name='source2')
        prefix = tmp_path / 'source@'
    else:
        prefix = bfile(tmp_path, [('1', 'dup', 10, 'A', 'G'), ('2', 'dup', 20, 'T', 'C')])
    result = run(script + f's <- load_sumstats("{path}",ref,match_by="id"); '
        f'g <- load_genotype("{prefix}",ref,match_by="id"); '
        'stopifnot(!any(is_present(s)),!any(is_present(g))); sub <- select_shards(ref,"1"); '
        f'g <- load_genotype("{prefix}",sub,match_by="id"); stopifnot(all(is_present(g))); '
        + fails(f'load_sumstats("{path}",sub,match_by="id")', 'Ambiguous'))
    assert result.stderr.count('2 slots in 1 duplicate ID key groups') == 2
    assert result.stderr.count('no matching source rows') == 2


def test_source_ambiguity_and_chromosome_scope(tmp_path):
    script = load_ref(tmp_path, [('1', 'conflict', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('conflict', 'A', 'G'), ('conflict', 'T', 'C')], CHR=['1', '2'])
    prefix = bfile(tmp_path, [('1', 'conflict', 10, 'A', 'G'), ('2', 'conflict', 20, 'T', 'C')])
    run(script + fails(f'load_sumstats("{path}",ref,match_by="id")', 'conflict') +
        f'g <- load_genotype("{prefix}",ref,match_by="id"); stopifnot(source_row0(g)==0L);')
    prefix = bfile(tmp_path, [('1', 'conflict', 10, 'A', 'G'), ('1', 'conflict', 20, 'T', 'C')], calls=[-1, -1])
    run(script + fails(f'load_genotype("{prefix}",ref,match_by="id")', 'conflict'))


@pytest.mark.parametrize('coordinates', [{}, {'CHR': ['chrBAD'], 'POS': ['nonsense']}])
def test_ignored_coordinates_and_missing_p(tmp_path, coordinates):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'A', 'G')], **coordinates)
    zipped = tmp_path / 'stats.tsv.gz'
    with gzip.open(zipped, 'wb') as out:
        out.write(path.read_bytes())
    run(script + f's <- load_sumstats("{zipped}",ref,match_by="id"); stopifnot(is_present(s));')
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'T', 'C')], P=[.01, 'NA'])
    run(script + f's <- load_sumstats("{path}",ref,match_by="id"); stopifnot(is_present(s));')
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'a', 'C')], P=[.01, 'NA'])
    run(script + fails(f'load_sumstats("{path}",ref,match_by="id")', 'uppercase DNA'))


@pytest.mark.parametrize('pair', [('A', 'T'), ('C', 'G'), ('A', 'AG'), ('AC', 'GT')])
def test_allele_exclusions_warn_empty_shard(tmp_path, pair):
    script = load_ref(tmp_path, [('1', 'id', 1, *pair)])
    path = stats_file(tmp_path, [('id', *pair)])
    prefix = bfile(tmp_path, [('1', 'id', 10, *pair)])
    result = run(script + f's <- load_sumstats("{path}",ref,match_by="id"); '
        f'g <- load_genotype("{prefix}",ref,match_by="id"); stopifnot(!any(is_present(s)), !any(is_present(g)));')
    assert 'no matching source rows' in result.stderr
    assert 'reference panel: excluded 1' in result.stderr
    assert f'{path}: excluded 1' in result.stderr


@pytest.mark.parametrize('marker,eligible', [
    ('!', True), ('~', True), ('x' * 32, True), ('x' * 33, False),
    ('.', False), ('é', False), ('\x7f', False), ('NA', False), ('NaN', False),
    ('na', True), ('nan', True), ('0', True),
])
def test_literal_id_boundaries(tmp_path, marker, eligible):
    script = load_ref(tmp_path, [('1', marker, 1, 'A', 'G')])
    path = stats_file(tmp_path, [(marker, 'A', 'G')])
    prefix = bfile(tmp_path, [('1', marker, 10, 'A', 'G')])
    result = run(script + f's <- load_sumstats("{path}",ref,match_by="id"); '
        f'g <- load_genotype("{prefix}",ref,match_by="id"); '
        f'stopifnot(is_present(s)=={str(eligible).upper()}, identical(is_present(s),is_present(g)), '
        f'identical(snp(ref), {json.dumps(marker)}));')
    assert ('unusable' in result.stderr) == (not eligible)


# fread strips surrounding whitespace by default; that divergence is accepted.
@pytest.mark.parametrize('marker', ['', 'i d'])
def test_source_ids_not_trimmed(tmp_path, marker):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [(marker, 'A', 'G')])
    result = run(script + f's <- load_sumstats("{path}",ref,match_by="id"); stopifnot(!is_present(s));')
    assert 'unusable' in result.stderr


def test_case_sensitive_no_fallback_and_quiet_warnings(tmp_path):
    script = load_ref(tmp_path, [('1', 'rs123', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('RS123', 'A', 'G')], CHR=['1'], POS=[1])
    run(script + f's <- load_sumstats("{path}",ref,match_by="id"); stopifnot(!is_present(s)); '
        f's <- load_sumstats("{path}",ref); stopifnot(is_present(s));')
    path = stats_file(tmp_path, [('.', 'A', 'G'), ('pal', 'A', 'T')])
    result = run(script + 'set_verbosity("quiet"); ' +
        f's <- load_sumstats("{path}",ref,match_by="id"); stopifnot(!is_present(s));')
    assert 'unusable' in result.stderr and 'excluded 1 palindromic' not in result.stderr


@pytest.mark.parametrize('mode', ['genomic_position', 'id'])
def test_identical_alleles_rejected_before_eligibility(tmp_path, mode):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('.', 'A', 'A')], CHR=['1'], POS=[1], P=['NA'])
    prefix = bfile(tmp_path, [('1', '.', 1, 'A', 'A')])
    for loader, source in [('load_sumstats', path), ('load_genotype', prefix)]:
        script += fails(f'{loader}("{source}",ref,match_by="{mode}")', 'a1 and a2 must differ')
    run(script)


def test_invalid_options_before_io():
    script = ''
    for loader in ('load_sumstats', 'load_genotype'):
        for options in ('match_by="ID"', 'match_by=NULL', 'match_by=NA', 'match_by=c("id","id")',
                        'match_by="id",allow_indel_swaps=TRUE'):
            script += fails(f'{loader}("does-not-exist",NULL,{options})', 'match_by')
    run(script)


def test_unused_duplicates_and_empty_source(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'T', 'C'), ('id', 'A', 'T'), ('id', 'AC', 'GT'),
                                ('unused', 'A', 'G'), ('unused', 'G', 'A'), ('.', 'A', 'G')],
                      P=[.01] * 5 + ['NA'], beta=[float('nan')] * 6)
    prefix = bfile(tmp_path, [('1', 'id', 10, 'T', 'C'), ('1', 'unused', 20, 'A', 'G'),
                             ('1', 'unused', 30, 'G', 'A')])
    result = run(script + f's <- load_sumstats("{path}",ref,match_by="id"); '
        f'g <- load_genotype("{prefix}",ref,match_by="id"); '
        'stopifnot(is_present(s), is_present(g), is.nan(beta_vec(s)), is.null(eaf_vec(s)));')
    assert 'unusable' not in result.stderr
    path = stats_file(tmp_path, [('id', 'A', 'G')], P=['NA'])
    result = run(script + f's <- load_sumstats("{path}",ref,match_by="id"); stopifnot(!is_present(s));')
    assert 'no matching source rows' in result.stderr


def test_sharded_offsets_chr_x_samples_and_ploidy(tmp_path):
    script = load_ref(tmp_path, [('1', 'first', 1, 'A', 'G'), ('X', 'second', 1, 'A', 'C')])
    bfile(tmp_path, [('1', 'unused', 10, 'A', 'G'), ('1', 'first', 20, 'A', 'G')], name='source1')
    x = bfile(tmp_path, [('X', 'second', 10, 'G', 'T')], name='sourceX')
    x.with_suffix('.ploidy').write_text('1\t2\n')
    (tmp_path / 'source1.fam').write_text('f i 0 0 1 -9\nf j 0 0 2 -9\n')
    run(script + f'g <- load_genotype("{tmp_path / "source@"}",ref,match_by="id"); '
        'stopifnot(identical(source_row0(g),c(1L,0L)), identical(is_allele_swapped(g),c(FALSE,TRUE)), '
        'identical(fetch_genotypes(g,c(2,1),haploid_mode="ploidy_scaled"), matrix(c(1,NaN,0,2),2)), '
        'is_allele_swapped(select_shards(g,"X")));')


def test_genotype_exclusion_diagnostics_chromosome_scope(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    prefix = bfile(tmp_path, [('1', 'id', 10, 'A', 'G'), ('1', '.', 20, 'A', 'G'),
                             ('1', 'pal', 30, 'A', 'T'), ('2', '.', 10, 'A', 'G'),
                             ('2', 'pal', 20, 'A', 'T'), ('2', 'multi', 30, 'A', 'AG')])
    result = run(script + f'g <- load_genotype("{prefix}",ref,match_by="id"); stopifnot(is_present(g));')
    assert f'{prefix}.bim: excluded 1 unusable missing, empty, or placeholder ID(s)' in result.stderr
    assert f'{prefix}.bim: excluded 1 palindromic SNV(s)' in result.stderr
    assert 'excluded 2' not in result.stderr and f'{prefix}.bim: excluded 1 multi-base' not in result.stderr


def test_only_id_mode_decodes_ids_and_neither_mode_decodes_alleles(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'C', 'T')], CHR=['1'], POS=[1])
    prefix = bfile(tmp_path, [('1', 'id', 1, 'C', 'T')])
    run(script + 'ref$shards[[1]]$a1_text <- NULL; ref$shards[[1]]$a2_text <- NULL; '
        'original <- .decode_text_payload; calls <- 0L; .decode_text_payload <- function(payload, expected_len, field) { '
        'stopifnot(field=="snp"); calls <<- calls + 1L; original(payload,expected_len,field) }; '
        f's <- suppressWarnings(load_sumstats("{path}",ref)); g <- load_genotype("{prefix}",ref); stopifnot(calls==0L); '
        f's <- load_sumstats("{path}",ref,match_by="id"); stopifnot(calls==1L,zvec(s)== -3); '
        f'g <- load_genotype("{prefix}",ref,match_by="id"); stopifnot(calls==2L,is_allele_swapped(g));')


def test_missing_snp_column(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = tmp_path / 'stats.tsv'
    path.write_text('CHR\tPOS\ta1\ta2\tP\n1\t1\tA\tG\t.01\n')
    run(script + fails(f'load_sumstats("{path}",ref,match_by="id")', 'missing required columns: snp'))
