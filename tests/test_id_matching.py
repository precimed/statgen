"""Python source-loader contract for exact ID and allele-group matching."""
import warnings

import numpy as np
import pandas as pd
import pytest

from statgen import get_verbosity, set_verbosity
from statgen.genotype import load_genotype, load_genotype_cache
from statgen.reference import load_reference, load_reference_cache
from statgen.sumstats import load_sumstats, load_sumstats_cache


def reference(tmp_path, rows):
    path = tmp_path / 'reference.bim'
    path.write_text(''.join(f'{chrom}\t{name}\t0\t{bp}\t{a1}\t{a2}\n'
                            for chrom, name, bp, a1, a2 in rows))
    return load_reference(path)


def stats_file(tmp_path, rows, **columns):
    path = tmp_path / 'stats.tsv'
    frame = pd.DataFrame(rows, columns=['SNP', 'a1', 'a2'])
    frame['P'] = 0.01
    frame['z'] = 3.0
    frame['n'] = 100.0
    for name, values in columns.items():
        frame[name] = values
    frame.to_csv(path, sep='\t', index=False)
    return path


def bfile(tmp_path, rows, calls=None, name='source'):
    prefix = tmp_path / name
    prefix.with_suffix('.bim').write_text(''.join(
        f'{chrom}\t{marker}\t0\t{bp}\t{a1}\t{a2}\n'
        for chrom, marker, bp, a1, a2 in rows))
    prefix.with_suffix('.fam').write_text('f i 0 0 1 -9\n')
    calls = [0] * len(rows) if calls is None else calls
    encoding = {0: 3, 1: 2, 2: 0, -1: 1}
    prefix.with_suffix('.bed').write_bytes(b'\x6c\x1b\x01' + bytes(encoding[x] for x in calls))
    return prefix


@pytest.fixture(autouse=True)
def restore_verbosity():
    previous = get_verbosity()
    set_verbosity('info')
    yield
    set_verbosity(previous)


@pytest.mark.parametrize('ref_pair,pairs', [
    (('A', 'G'), [('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')]),
    (('A', 'C'), [('A', 'C'), ('C', 'A'), ('T', 'G'), ('G', 'T')]),
])
def test_orientations_and_cache_roundtrips(tmp_path, capsys, ref_pair, pairs):
    ref = reference(tmp_path, [('1', f'id{i}', i + 1, *ref_pair) for i in range(4)])
    ref.save_cache(tmp_path / 'ref.npz')
    ref = load_reference_cache(tmp_path / 'ref.npz')
    checksum = ref.shards[0].checksum
    path = stats_file(tmp_path, [(f'id{i}', *pair) for i, pair in enumerate(pairs)],
                      beta=[2.0] * 4, eaf=[0.2] * 4, se=[0.5] * 4, info=[0.9] * 4)
    stats = load_sumstats(path, ref, match_by='id')
    np.testing.assert_allclose(stats.zvec, [3, -3, 3, -3])
    np.testing.assert_allclose(stats.beta_vec, [2, -2, 2, -2])
    np.testing.assert_allclose(stats.eaf_vec, [.2, .8, .2, .8])
    for accessor, value in [('logpvec', 2), ('nvec', 100), ('se_vec', .5), ('info_vec', .9)]:
        np.testing.assert_allclose(getattr(stats, accessor), [value] * 4)
    prefix = bfile(tmp_path, [('1', f'id{i}', 101 + i, *pair) for i, pair in enumerate(pairs)],
                   calls=[0, 0, 2, -1])
    geno = load_genotype(prefix, ref, match_by='id')
    np.testing.assert_array_equal(geno.is_allele_swapped, [False, True, False, True])
    np.testing.assert_array_equal(geno.fetch_genotypes_int8([0, 1, 2, 3]), [[0, 2, 2, -1]])
    stats.save_cache(tmp_path / 'stats.npz')
    geno.save_cache(tmp_path / 'geno.npz')
    np.testing.assert_allclose(load_sumstats_cache(tmp_path / 'stats.npz').zvec, stats.zvec)
    restored = load_genotype_cache(tmp_path / 'geno.npz')
    np.testing.assert_array_equal(restored.source_row0, geno.source_row0)
    np.testing.assert_array_equal(restored.is_allele_swapped, geno.is_allele_swapped)
    np.testing.assert_array_equal(restored.fetch_genotypes_int8([0, 1, 2, 3]), [[0, 2, 2, -1]])
    assert ref.shards[0].checksum == checksum
    output = capsys.readouterr().err
    assert output.count('2 matched row(s) requiring inversion') == 2
    assert output.count('2 matched row(s) requiring strand complementation') == 2


def test_same_id_distinct_allele_groups_and_reference_order(tmp_path):
    ref = reference(tmp_path, [('1', 'same', 1, 'A', 'G'), ('1', 'same', 2, 'A', 'C')])
    path = stats_file(tmp_path, [('same', 'G', 'T'), ('same', 'T', 'C')], z=[5, 7])
    np.testing.assert_allclose(load_sumstats(path, ref, match_by='id').zvec, [7, -5])
    prefix = bfile(tmp_path, [('1', 'same', 100, 'G', 'T'), ('1', 'same', 200, 'T', 'C')])
    np.testing.assert_array_equal(load_genotype(prefix, ref, match_by='id').source_row0, [1, 0])


@pytest.mark.parametrize('sharded', [False, True])
def test_reference_duplicates_are_panel_wide_and_subset_dependent(tmp_path, sharded):
    ref = reference(tmp_path, [('1', 'dup', 1, 'A', 'G'), ('2', 'dup', 2, 'T', 'C')])
    path = stats_file(tmp_path, [('dup', 'A', 'G'), ('dup', 'G', 'A')])
    with pytest.warns(RuntimeWarning) as caught:
        stats = load_sumstats(path, ref, match_by='id')
    assert any('2 slots in 1 duplicate ID key groups' in str(w.message) for w in caught)
    assert sum('no matching source rows' in str(w.message) for w in caught) == 2
    assert not stats.is_present.any()  # Source duplicates for excluded slots do not fail.
    if sharded:
        bfile(tmp_path, [('1', 'dup', 10, 'A', 'G')], name='source1')
        bfile(tmp_path, [('2', 'dup', 20, 'T', 'C')], name='source2')
        prefix = tmp_path / 'source@'
    else:
        prefix = bfile(tmp_path, [('1', 'dup', 10, 'A', 'G'), ('2', 'dup', 20, 'T', 'C')])
    with pytest.warns(RuntimeWarning, match='duplicate ID key groups'):
        assert not load_genotype(prefix, ref, match_by='id').is_present.any()
    subset = ref.select_shards(['1'])
    assert load_genotype(prefix, subset, match_by='id').is_present.all()
    unique = stats_file(tmp_path, [('dup', 'T', 'C')])
    assert load_sumstats(unique, subset, match_by='id').is_present.all()


def test_source_ambiguity_and_chromosome_scope(tmp_path):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'T', 'C')], CHR=['1', '2'])
    with pytest.raises(ValueError, match='Ambiguous'):
        load_sumstats(path, ref, match_by='id')
    prefix = bfile(tmp_path, [('1', 'id', 100, 'A', 'G'), ('2', 'id', 200, 'T', 'C')])
    np.testing.assert_array_equal(load_genotype(prefix, ref, match_by='id').source_row0, [0])
    prefix = bfile(tmp_path, [('1', 'id', 100, 'A', 'G'), ('1', 'id', 200, 'T', 'C')], calls=[-1, -1])
    with pytest.raises(ValueError, match='Ambiguous'):
        load_genotype(prefix, ref, match_by='id')


@pytest.mark.parametrize('coordinates', [{}, {'CHR': ['chrBAD'], 'POS': ['nonsense']}])
def test_ignored_coordinates_and_missing_p_filtering(tmp_path, coordinates):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'A', 'G')], **coordinates)
    assert load_sumstats(path, ref, match_by='id').is_present.all()
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'T', 'C')], P=[.01, 'NA'])
    assert load_sumstats(path, ref, match_by='id').is_present.all()
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'a', 'C')], P=[.01, 'NA'])
    with pytest.raises(ValueError, match='uppercase DNA'):
        load_sumstats(path, ref, match_by='id')


@pytest.mark.parametrize('pair', [('A', 'T'), ('C', 'G'), ('A', 'AG'), ('AC', 'GT')])
def test_allele_exclusions_warn_empty_shard(tmp_path, capsys, pair):
    ref = reference(tmp_path, [('1', 'id', 1, *pair)])
    path = stats_file(tmp_path, [('id', *pair)])
    with pytest.warns(RuntimeWarning, match='no matching source rows'):
        assert not load_sumstats(path, ref, match_by='id').is_present.any()
    prefix = bfile(tmp_path, [('1', 'id', 1, *pair)])
    assert not load_genotype(prefix, ref, match_by='id').is_present.any()
    output = capsys.readouterr().err
    assert 'reference panel: excluded 1' in output
    assert f'{path}: excluded 1' in output


@pytest.mark.parametrize('marker,eligible', [
    ('!', True), ('~', True), ('x' * 32, True), ('x' * 33, False),
    ('.', False), ('é', False), ('\x7f', False), ('NA', False), ('NaN', False),
    ('na', True), ('nan', True), ('0', True),
])
def test_id_boundaries_in_both_loaders(tmp_path, marker, eligible):
    ref = reference(tmp_path, [('1', marker, 1, 'A', 'G')])
    assert ref.snp[0] == marker
    path = stats_file(tmp_path, [(marker, 'A', 'G')])
    prefix = bfile(tmp_path, [('1', marker, 2, 'A', 'G')])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        stats = load_sumstats(path, ref, match_by='id')
        geno = load_genotype(prefix, ref, match_by='id')
    assert stats.is_present[0] == eligible
    assert geno.is_present[0] == eligible
    assert any('unusable' in str(w.message) for w in caught) == (not eligible)


@pytest.mark.parametrize('marker', ['', ' id', 'id ', 'i d'])
def test_source_ids_are_not_trimmed(tmp_path, marker):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [(marker, 'A', 'G')])
    with pytest.warns(RuntimeWarning) as caught:
        assert not load_sumstats(path, ref, match_by='id').is_present.any()
    assert any('unusable' in str(w.message) for w in caught)


def test_case_sensitive_no_coordinate_fallback_and_quiet_warnings(tmp_path, capsys):
    ref = reference(tmp_path, [('1', 'rs123', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('RS123', 'A', 'G')], CHR=['1'], POS=[1])
    with pytest.warns(RuntimeWarning, match='no matching source rows'):
        assert not load_sumstats(path, ref, match_by='id').is_present.any()
    assert load_sumstats(path, ref).is_present.all()
    capsys.readouterr()
    set_verbosity('quiet')
    path = stats_file(tmp_path, [('.', 'A', 'G')])
    with pytest.warns(RuntimeWarning) as caught:
        load_sumstats(path, ref, match_by='id')
    assert any('unusable' in str(w.message) for w in caught)
    assert capsys.readouterr().err == ''


@pytest.mark.parametrize('mode', ['genomic_position', 'id'])
def test_identical_alleles_rejected_before_eligibility(tmp_path, mode):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('.', 'A', 'A')], CHR=['1'], POS=[1], P=['NA'])
    prefix = bfile(tmp_path, [('1', '.', 1, 'A', 'A')])
    for loader, source in [(load_sumstats, path), (load_genotype, prefix)]:
        with pytest.raises(ValueError, match='a1 and a2 must differ'):
            loader(source, ref, match_by=mode)


@pytest.mark.parametrize('loader', [load_sumstats, load_genotype])
@pytest.mark.parametrize('mode,swaps', [('ID', False), ('position', False), (None, False), ([], False), ('id', True)])
def test_invalid_options_before_io(loader, mode, swaps):
    with pytest.raises(ValueError, match='match_by'):
        loader('does-not-exist', None, match_by=mode, allow_indel_swaps=swaps)


def test_position_mode_does_not_access_reference_ids_or_alleles(tmp_path, monkeypatch):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('ignored', 'A', 'G')], CHR=['1'], POS=[1])
    prefix = bfile(tmp_path, [('1', 'ignored', 1, 'A', 'G')])
    def unexpected(_):
        raise AssertionError('lazy strings accessed')
    for cls in (type(ref), type(ref.shards[0])):
        for field in ('snp', 'a1', 'a2'):
            monkeypatch.setattr(cls, field, property(unexpected))
    assert load_sumstats(path, ref).is_present.all()
    assert load_genotype(prefix, ref).is_present.all()
    assert load_sumstats(path, ref, match_by='genomic_position').is_present.all()


def test_ineligible_and_unused_source_duplicates_do_not_compete(tmp_path, capsys):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'T', 'C'), ('id', 'A', 'T'), ('id', 'A', 'AG'),
                                ('unused', 'A', 'G'), ('unused', 'G', 'A'), ('.', 'A', 'G')],
                      P=[.01, .01, .01, .01, .01, 'NA'], beta=[np.nan] * 6, eaf=[np.nan] * 6)
    with warnings.catch_warnings(record=True) as caught:
        stats = load_sumstats(path, ref, match_by='id')
    assert not caught  # Missing-P exclusion happens before unusable-ID reporting.
    assert stats.is_present.all()
    assert np.isnan(stats.beta_vec[0]) and np.isnan(stats.eaf_vec[0])
    output = capsys.readouterr().err
    assert 'excluded 1 palindromic SNV' in output
    assert 'excluded 1 multi-base variant' in output
    prefix = bfile(tmp_path, [('1', 'id', 10, 'T', 'C'), ('1', 'id', 20, 'A', 'T'),
                              ('1', 'unused', 30, 'A', 'G'), ('1', 'unused', 40, 'G', 'A')])
    assert load_genotype(prefix, ref, match_by='id').is_present.all()


def test_empty_source_after_missing_p_filter(tmp_path):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'A', 'G')], P=['NA'])
    with pytest.warns(RuntimeWarning, match='no matching source rows'):
        assert not load_sumstats(path, ref, match_by='id').is_present.any()


def test_id_mode_uses_only_reference_allele_hashes(tmp_path, monkeypatch):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'C', 'T')])
    prefix = bfile(tmp_path, [('1', 'id', 10, 'C', 'T')])
    def unexpected(_):
        raise AssertionError('allele strings accessed')
    for cls in (type(ref), type(ref.shards[0])):
        for field in ('a1', 'a2'):
            monkeypatch.setattr(cls, field, property(unexpected))
    np.testing.assert_allclose(load_sumstats(path, ref, match_by='id').zvec, [-3])
    assert load_genotype(prefix, ref, match_by='id').is_allele_swapped.all()


def test_genotype_sharded_offsets_and_ploidy(tmp_path):
    ref = reference(tmp_path, [('1', 'first', 1, 'A', 'G'), ('X', 'second', 1, 'A', 'C')])
    bfile(tmp_path, [('1', 'unused', 10, 'A', 'G'), ('1', 'first', 20, 'A', 'G')], name='source1')
    prefix_x = bfile(tmp_path, [('X', 'second', 10, 'G', 'T')], name='sourceX')
    prefix_x.with_suffix('.ploidy').write_text('1\t2\n')
    geno = load_genotype(tmp_path / 'source@', ref, match_by='id')
    np.testing.assert_array_equal(geno.source_row0, [1, 0])
    np.testing.assert_array_equal(geno.is_allele_swapped, [False, True])
    np.testing.assert_array_equal(geno.fetch_genotypes_int8([1, 0]), [[2, 0]])
    np.testing.assert_allclose(geno.fetch_genotypes([1, 0], haploid_mode='ploidy_scaled'), [[1, 0]])


def test_id_mode_requires_snp_column(tmp_path):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'A', 'G')])
    frame = pd.read_csv(path, sep='\t').drop(columns='SNP')
    frame.to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match='missing required columns: snp'):
        load_sumstats(path, ref, match_by='id')


def test_genotype_exclusion_diagnostics_use_reference_chromosome_scope(tmp_path, capsys):
    ref = reference(tmp_path, [('1', 'id', 1, 'A', 'G')])
    prefix = bfile(tmp_path, [
        ('1', 'id', 10, 'A', 'G'), ('1', '.', 20, 'A', 'G'),
        ('1', 'pal', 30, 'A', 'T'), ('1', 'multi', 40, 'A', 'AG'),
        ('2', '.', 10, 'A', 'G'), ('2', 'x' * 33, 20, 'A', 'G'),
        ('2', 'pal', 30, 'A', 'T'), ('2', 'multi', 40, 'A', 'AG'),
    ])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        geno = load_genotype(prefix, ref, match_by='id')
    assert geno.is_present.all()
    assert len(caught) == 1
    assert 'excluded 1 unusable missing, empty, or placeholder ID(s)' in str(caught[0].message)
    output = capsys.readouterr().err
    assert f'{prefix}.bim: excluded 1 palindromic SNV(s)' in output
    assert f'{prefix}.bim: excluded 1 multi-base variant(s)' in output
    assert 'excluded 2' not in output

    # Ordinary input validation still covers chromosomes outside the projection.
    bim_path = prefix.with_suffix('.bim')
    bim_path.write_text(bim_path.read_text().replace('2\tmulti\t0\t40\tA\tAG',
                                                   '2\tmulti\t0\t40\tA\tag'))
    with pytest.raises(ValueError, match='uppercase DNA'):
        load_genotype(prefix, ref, match_by='id')
