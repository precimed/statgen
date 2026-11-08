"""MATLAB/Octave ID matching, including packed reference-name conversion."""
import gzip

import pytest

from tests.conftest import matlab_output_text, run_octave, skipif_no_octave
from tests.test_id_matching import bfile, reference, stats_file

pytestmark = [pytest.mark.octave, skipif_no_octave]


def run(script):
    result = run_octave(script, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    return matlab_output_text(result.stdout + result.stderr)


def fails(call, message):
    message = message.replace("'", "''")
    return (f"failed = false; try; {call}; catch ME; "
            f"failed = ~isempty(strfind(ME.message, '{message}')); end; assert(failed); ")


def load_ref(tmp_path, rows, cached=True):
    reference(tmp_path, rows)
    script = f"ref = statgen.load_reference('{tmp_path / 'reference.bim'}'); "
    if cached:
        script += (f"ref.save_cache('{tmp_path / 'reference.mat'}'); "
                   f"ref = statgen.load_reference_cache('{tmp_path / 'reference.mat'}'); ")
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
    prefix = bfile(tmp_path, [('1', f'id{i}', 100 + i, *p) for i, p in enumerate(pairs)],
                   calls=[0, 0, 2, -1])
    output = run(script +
        f"s = statgen_load_sumstats('{path}', ref, 'match_by', 'id'); "
        f"g = statgen_load_genotype('{prefix}', ref, 'match_by', 'id'); "
        "assert(isequal(s.zvec, [3;-3;3;-3])); assert(isequal(s.beta_vec, [2;-2;2;-2])); "
        "assert(max(abs(s.eaf_vec - [.2;.8;.2;.8])) < 1e-12); "
        "assert(isequal(s.logpvec, 2 * ones(4,1))); assert(all(s.nvec == 100)); "
        "assert(all(s.se_vec == .5) && all(s.info_vec == .9)); "
        "assert(isequal(g.is_allele_swapped, logical([0;1;0;1]))); "
        "assert(isequal(g.fetch_genotypes_int8([1 2 3 4 2]), int8([0 2 2 -1 2]))); "
        f"g.save_cache('{tmp_path / 'g.mat'}'); h = statgen.load_genotype_cache('{tmp_path / 'g.mat'}'); "
        "assert(isequal(h.is_allele_swapped, g.is_allele_swapped)); "
        "assert(isequal(h.source_row0, g.source_row0)); "
        "assert(isequal(h.fetch_genotypes_int8([4 2 1]), int8([-1 2 0]))); "
        f"s.save_cache('{tmp_path / 's.mat'}'); t = statgen.load_sumstats_cache('{tmp_path / 's.mat'}'); "
        "assert(isequaln(t.zvec, s.zvec) && isequaln(t.eaf_vec, s.eaf_vec)); "
        "assert(ref.is_object_compatible(g)); assert(ref.is_object_compatible(s));")
    assert output.count('2 matched row(s) requiring inversion') == 2
    assert output.count('2 matched row(s) requiring strand complementation') == 2


def test_distinct_groups_source_order_and_shared_width(tmp_path):
    script = load_ref(tmp_path, [('1', 'same', 1, 'A', 'G'), ('2', 'same', 2, 'A', 'C'),
                                 ('2', 'x' * 32, 3, 'A', 'G')])
    path = stats_file(tmp_path, [('same', 'G', 'T'), ('x' * 32, 'A', 'G'), ('same', 'T', 'C')], z=[5, 9, 7])
    prefix = bfile(tmp_path, [('1', 'same', 10, 'T', 'C'), ('2', 'same', 20, 'G', 'T'),
                             ('2', 'x' * 32, 30, 'A', 'G')])
    run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); "
        f"g = statgen.load_genotype('{prefix}', ref, 'match_by', 'id'); "
        "assert(isequal(s.zvec, [7;-5;9])); assert(isequal(g.source_row0, [0;1;2]));")


@pytest.mark.parametrize('sharded', [False, True])
def test_panel_wide_duplicates_and_subset(tmp_path, sharded):
    script = load_ref(tmp_path, [('1', 'dup', 1, 'A', 'G'), ('2', 'dup', 2, 'T', 'C')])
    path = stats_file(tmp_path, [('dup', 'A', 'G'), ('dup', 'G', 'A')])
    if sharded:
        bfile(tmp_path, [('1', 'dup', 10, 'A', 'G')], name='source1')
        bfile(tmp_path, [('2', 'dup', 20, 'T', 'C')], name='source2')
        prefix = tmp_path / 'source@'
    else:
        prefix = bfile(tmp_path, [('1', 'dup', 10, 'A', 'G'), ('2', 'dup', 20, 'T', 'C')])
    output = run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); "
        f"g = statgen.load_genotype('{prefix}', ref, 'match_by', 'id'); "
        "assert(~any(s.is_present) && ~any(g.is_present)); sub = ref.select_shards({'1'}); "
        f"g = statgen.load_genotype('{prefix}', sub, 'match_by', 'id'); assert(all(g.is_present)); "
        + fails(f"statgen.load_sumstats('{path}', sub, 'match_by', 'id')", 'Ambiguous'))
    assert output.count('2 slots in 1 duplicate ID key groups') == 2
    assert output.count('no matching source rows') == 2


def test_source_ambiguity_chromosome_scope_and_missing_p(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'T', 'C')], CHR=['1', '2'])
    prefix = bfile(tmp_path, [('1', 'id', 10, 'A', 'G'), ('2', 'id', 20, 'T', 'C')])
    run(script + fails(f"statgen.load_sumstats('{path}', ref, 'match_by', 'id')", "for ID 'id'") +
        f"g = statgen.load_genotype('{prefix}', ref, 'match_by', 'id'); assert(g.source_row0 == 0);")
    prefix = bfile(tmp_path, [('1', 'id', 10, 'A', 'G'), ('1', 'id', 20, 'T', 'C')], calls=[-1, -1])
    path = stats_file(tmp_path, [('id', 'A', 'G'), ('id', 'T', 'C')], P=[.01, 'NA'])
    run(script + fails(f"statgen.load_genotype('{prefix}', ref, 'match_by', 'id')", "for ID 'id'") +
        f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); assert(s.is_present);")


def test_sharded_ambiguity_names_id_and_source(tmp_path):
    script = load_ref(tmp_path, [('1', 'first', 1, 'A', 'G'), ('2', 'conflict', 1, 'A', 'C')])
    bfile(tmp_path, [('1', 'first', 10, 'A', 'G')], name='source1')
    bfile(tmp_path, [('2', 'conflict', 10, 'A', 'C'), ('2', 'conflict', 20, 'G', 'T')], name='source2')
    call = f"statgen.load_genotype('{tmp_path / 'source@'}', ref, 'match_by', 'id')"
    run(script + fails(call, "for ID 'conflict'") + fails(call, str(tmp_path / 'source2.bim')))


def test_allele_exclusions_unused_duplicates_and_empty_source(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G'), ('1', 'pal', 2, 'C', 'G'),
                                 ('2', 'multi', 1, 'A', 'AG')])
    path = stats_file(tmp_path, [('id', 'T', 'C'), ('id', 'A', 'T'), ('id', 'AC', 'GT'),
                                ('unused', 'A', 'G'), ('unused', 'G', 'A'), ('.', 'A', 'G')],
                      P=[.01] * 5 + ['NA'], beta=[float('nan')] * 6)
    prefix = bfile(tmp_path, [('1', 'id', 10, 'T', 'C'), ('1', 'pal', 20, 'C', 'G'),
                             ('2', 'multi', 10, 'A', 'AG')])
    output = run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); "
        f"g = statgen.load_genotype('{prefix}', ref, 'match_by', 'id'); "
        "assert(isequal(s.is_present, logical([1;0;0]))); assert(isempty(s.eaf_vec)); "
        "assert(all(isnan(s.beta_vec))); assert(isequal(g.is_present, s.is_present));")
    assert 'unusable' not in output
    assert 'reference panel: excluded 1 palindromic SNV(s)' in output
    assert f'{path}: excluded 1 multi-base variant(s)' in output
    assert 'no matching source rows' in output
    path = stats_file(tmp_path, [('id', 'A', 'G')], P=['NA'])
    run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); assert(~any(s.is_present));")


@pytest.mark.parametrize('marker,eligible', [
    ('!', True), ('~', True), ('x' * 32, True), ('x' * 33, False),
    ('.', False), ('é', False), ('\x7f', False), ('NA', False), ('NaN', False),
    ('na', True), ('nan', True), ('0', True),
])
def test_literal_id_boundaries(tmp_path, marker, eligible):
    script = load_ref(tmp_path, [('1', marker, 1, 'A', 'G')])
    path = stats_file(tmp_path, [(marker, 'A', 'G')])
    prefix = bfile(tmp_path, [('1', marker, 10, 'A', 'G')])
    output = run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); "
        f"g = statgen.load_genotype('{prefix}', ref, 'match_by', 'id'); "
        f"assert(s.is_present == {str(eligible).lower()}); assert(g.is_present == s.is_present);")
    assert ('unusable' in output) == (not eligible)


@pytest.mark.parametrize('marker', ['', ' id', 'id ', 'i d'])
def test_source_ids_not_trimmed(tmp_path, marker):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = stats_file(tmp_path, [(marker, 'A', 'G')])
    output = run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); assert(~s.is_present);")
    assert 'unusable' in output


def test_ignored_coordinates_gzip_case_and_quiet(tmp_path):
    script = load_ref(tmp_path, [('1', 'rs123', 1, 'A', 'G')])
    path = stats_file(tmp_path, [('rs123', 'A', 'G')], CHR=['chrBAD'], POS=['nonsense'])
    zipped = tmp_path / 'stats.tsv.gz'
    with gzip.open(zipped, 'wb') as out:
        out.write(path.read_bytes())
    run(script + f"s = statgen.load_sumstats('{zipped}', ref, 'match_by', 'id'); assert(s.is_present);")
    path = stats_file(tmp_path, [('RS123', 'A', 'G')], CHR=['1'], POS=[1])
    run(script + f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); assert(~s.is_present); "
        f"s = statgen.load_sumstats('{path}', ref); assert(s.is_present);")
    path = stats_file(tmp_path, [('.', 'A', 'G'), ('pal', 'A', 'T')])
    output = run(script + "statgen.set_verbosity('quiet'); " +
        f"s = statgen.load_sumstats('{path}', ref, 'match_by', 'id'); assert(~s.is_present);")
    assert 'unusable' in output and 'excluded 1 palindromic' not in output


def test_invalid_options_and_validation_before_eligibility(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    for loader in ('load_sumstats', 'load_genotype'):
        for option in ("'match_by', 'ID'", "'match_by', []", "'match_by', 1",
                       "'match_by', 'id', 'allow_indel_swaps', true",
                       "'allow_indel_swaps', true, 'match_by', 'id'"):
            script += fails(f"statgen.{loader}('does-not-exist', ref, {option})", 'match_by')
    path = stats_file(tmp_path, [('.', 'A', 'A')], CHR=['1'], POS=[1], P=['NA'])
    prefix = bfile(tmp_path, [('1', '.', 1, 'A', 'A')])
    for mode in ('id', 'genomic_position'):
        for loader, source in [('load_sumstats', path), ('load_genotype', prefix)]:
            script += fails(f"statgen.{loader}('{source}', ref, 'match_by', '{mode}')", 'a1 and a2 must differ')
    run(script)
    path = stats_file(tmp_path, [('id', 'a', 'C')], P=['NA'])
    run(load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')]) +
        fails(f"statgen.load_sumstats('{path}', ref, 'match_by', 'id')", 'uppercase DNA'))


def test_sharded_offsets_chr_x_samples_and_ploidy(tmp_path):
    script = load_ref(tmp_path, [('1', 'first', 1, 'A', 'G'), ('X', 'second', 1, 'A', 'C')])
    bfile(tmp_path, [('1', 'unused', 10, 'A', 'G'), ('1', 'first', 20, 'A', 'G')], name='source1')
    x = bfile(tmp_path, [('X', 'second', 10, 'G', 'T')], name='sourceX')
    x.with_suffix('.ploidy').write_text('1\t2\n')
    # Autosomal sample axis has a second subject absent from chrX.
    (tmp_path / 'source1.fam').write_text('f i 0 0 1 -9\nf j 0 0 2 -9\n')
    run(script + f"g = statgen.load_genotype('{tmp_path / 'source@'}', ref, 'match_by', 'id'); "
        "assert(isequal(g.source_row0, [1;0])); assert(isequal(g.is_allele_swapped, logical([0;1]))); "
        "assert(isequaln(g.fetch_genotypes([2 1], 'ploidy_scaled'), [1 0;NaN 2])); "
        "h = g.select_shards({'X'}); assert(h.is_allele_swapped);")


def test_genotype_diagnostic_chromosome_scope(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    prefix = bfile(tmp_path, [('1', 'id', 10, 'A', 'G'), ('1', '.', 20, 'A', 'G'),
                             ('1', 'pal', 30, 'A', 'T'), ('2', '.', 10, 'A', 'G'),
                             ('2', 'pal', 20, 'A', 'T'), ('2', 'multi', 30, 'A', 'AG')])
    output = run(script + f"g = statgen.load_genotype('{prefix}', ref, 'match_by', 'id'); assert(g.is_present);")
    assert f'{prefix}.bim: excluded 1 unusable missing, empty, or placeholder ID(s)' in output
    assert f'{prefix}.bim: excluded 1 palindromic SNV(s)' in output
    assert 'excluded 2' not in output and f'{prefix}.bim: excluded 1 multi-base' not in output


def test_missing_snp_and_position_mode_skips_ids(tmp_path):
    script = load_ref(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path = tmp_path / 'stats.tsv'
    path.write_text('CHR\tPOS\ta1\ta2\tP\n1\t1\tA\tG\t.01\n')
    run(script + fails(f"statgen.load_sumstats('{path}', ref, 'match_by', 'id')", 'missing required column: snp'))
    # Malformed packed strings must remain untouched by position matching;
    # allele strings are never decoded by either matcher.
    prefix = bfile(tmp_path, [('1', 'id', 1, 'A', 'G')])
    path.write_text('SNP\tCHR\tPOS\ta1\ta2\tP\nignored\t1\t1\tA\tG\t.01\n')
    run(script +
        "r = ref.shards{1}; bad = statgen.ReferenceShard.from_cache_text('1', 1, r.bp, "
        "r.a1_hash64, r.a2_hash64, r.checksum, sprintf('bad\\nids'), '', ''); "
        "broken = statgen.ReferencePanel({bad}); "
        f"s = statgen.load_sumstats('{path}', broken); assert(s.is_present); "
        f"g = statgen.load_genotype('{prefix}', broken); assert(g.is_present); "
        "good = statgen.ReferenceShard.from_cache_text('1', 1, r.bp, r.a1_hash64, r.a2_hash64, r.checksum, 'id', '', ''); "
        "hashed = statgen.ReferencePanel({good}); "
        f"g = statgen.load_genotype('{prefix}', hashed, 'match_by', 'id'); assert(g.is_present);")


def test_packed_conversion_batches_and_invalid_names():
    # A second batch, uneven lengths, empty and invalid IDs, and no truncation.
    run("values = repmat({'id'}, 65540, 1); values{1} = ''; values{2} = '.'; "
        "values{3} = repmat('x', 1, 33); values{4} = 'id '; values{65537} = 'last'; "
        "payload = strjoin(values', sprintf('\\n')); "
        "[rows, names, missing, invalid] = statgen.internal.id_name_matrix(payload, numel(values), true(numel(values),1), true); "
        "assert(missing == 2 && invalid == 2); assert(isequal(rows, (5:65540)')); "
        "assert(size(names, 2) == 32); assert(strcmp(deblank(names(end-3,:)), 'last')); "
        "[r, m, nm, ni] = statgen.internal.id_name_matrix('', 1, true, true); "
        "assert(isempty(r) && isempty(m) && nm == 1 && ni == 0); "
        "[r, m] = statgen.internal.id_name_matrix('', 0, false(0,1), true); assert(isempty(r) && isempty(m));")
