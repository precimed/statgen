import pytest

from tests.conftest import FIXTURES_DIR, copy_sharded_genotype, matlab_data_lines, run_octave, skipif_no_octave


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_cache_roundtrip_preserves_metadata_and_fetch(tmp_path):
    cache = tmp_path / "genotype_cache.mat"
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        f"statgen.save_genotype_cache(g, '{cache}', 'format', 'v5'); "
        f"h = statgen.load_genotype_cache('{cache}'); "
        "same_meta = strcmp(h.source_layout, g.source_layout) && h.num_snp == g.num_snp && h.num_sample == g.num_sample; "
        "same_meta = same_meta && isequal(h.is_present, g.is_present) && isequaln(h.ploidy_male, g.ploidy_male); "
        "same_meta = same_meta && isequal(h.source_row0, g.source_row0) && isequal(h.fid, g.fid); "
        "same_meta = same_meta && isequal(h.is_subject_present('X'), g.is_subject_present('X')); "
        "same_fetch = isequal(h.fetch_genotypes_int8([6 1 6 3]), g.fetch_genotypes_int8([6 1 6 3])); "
        "gf = h.fetch_genotypes([6 1 6 3]); "
        "fprintf('%d %d %s %.0f\\n', same_meta, same_fetch, class(gf), sum(isnan(gf(:))));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1 double 0"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_cache_top_level_variables_and_format_option(tmp_path):
    cache = tmp_path / "genotype_cache_v5.mat"
    expected_vars = [
        "father_id",
        "fid",
        "iid",
        "is_allele_swapped",
        "is_female",
        "is_male",
        "is_present",
        "metadata",
        "mother_id",
        "ploidy_female",
        "ploidy_male",
        "sex",
        "source_row0",
        "source_subject_row0",
        "subject_present",
    ]
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        f"statgen.save_genotype_cache(g, '{cache}', 'format', 'v5'); "
        f"payload = load('{cache}'); "
        "names = sort(fieldnames(payload)); "
        "fprintf('%s\\n', strjoin(names, ',')); "
        "fprintf('%s %d %d %d %d\\n', payload.metadata.schema, payload.metadata.n_shards, "
        "numel(payload.is_present), size(payload.subject_present, 1), size(payload.subject_present, 2));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        ",".join(expected_vars),
        "genotype_cache/0.2 2 8 4 2",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_cache_loads_without_source_sidecars_and_fetches_with_override(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    copy_sharded_genotype(source_dir)
    cache = tmp_path / "genotype_cache.mat"

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{source_dir / '@'}', ref); "
        f"statgen.save_genotype_cache(g, '{cache}'); "
        f"delete('{source_dir / '1.bed'}'); delete('{source_dir / '1.bim'}'); delete('{source_dir / '1.fam'}'); "
        f"delete('{source_dir / 'X.bed'}'); delete('{source_dir / 'X.bim'}'); delete('{source_dir / 'X.fam'}'); "
        f"delete('{source_dir / 'X.ploidy'}'); "
        f"h = statgen.load_genotype_cache('{cache}'); "
        "ok_load = strcmp(h.source_layout, 'sharded') && h.num_snp == 8 && h.num_sample == 4; "
        f"gi = h.fetch_genotypes_int8([1 6], '{FIXTURES_DIR / 'genotype/sharded/@.bed'}'); "
        "fprintf('%d %d %d\\n', ok_load, size(gi, 1), size(gi, 2)); "
        "fprintf('%d ', gi); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == [
        "1 4 2",
        "2 2 2 2 2 2 2 2",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_cache_subset_preserves_sample_axis_source_layout_and_fetch(tmp_path):
    cache = tmp_path / "genotype_cache.mat"
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        f"statgen.save_genotype_cache(g, '{cache}'); "
        f"h = statgen.load_genotype_cache('{cache}', {{'X'}}); "
        f"gi = h.fetch_genotypes_int8([1 3 1], '{FIXTURES_DIR / 'genotype/sharded/@.bed'}'); "
        "fprintf('%s %s %d %d\\n', h.source_layout, h.shards{1}.label, h.num_snp, h.num_sample); "
        "fprintf('%d ', h.is_subject_present('X')); fprintf('\\n'); "
        "fprintf('%d ', gi); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == [
        "sharded X 3 4",
        "1 1 1 1",
        "2 2 2 2 2 2 2 2 2 2 2 2",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_cache_validation_failures(tmp_path):
    bad_schema = tmp_path / "bad_schema.mat"
    bad_offsets = tmp_path / "bad_offsets.mat"
    bad_zero_shards = tmp_path / "bad_zero_shards.mat"
    unknown = tmp_path / "valid.mat"
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        f"statgen.save_genotype_cache(g, '{bad_schema}'); "
        f"statgen.save_genotype_cache(g, '{bad_offsets}'); "
        f"statgen.save_genotype_cache(g, '{bad_zero_shards}'); "
        f"statgen.save_genotype_cache(g, '{unknown}'); "
        f"payload = load('{bad_schema}'); payload.metadata.schema = 'bad'; "
        f"save('{bad_schema}', '-struct', 'payload'); "
        f"payload = load('{bad_offsets}'); payload.metadata.shard_start0(2) = 99; "
        f"save('{bad_offsets}', '-struct', 'payload'); "
        f"payload = load('{bad_zero_shards}'); payload.metadata.n_shards = 0; "
        f"payload.metadata.shard_labels = {{}}; payload.metadata.shard_checksums = {{}}; "
        f"payload.metadata.shard_start0 = []; payload.metadata.shard_stop0 = []; "
        f"save('{bad_zero_shards}', '-struct', 'payload'); "
        f"e1 = 0; try; statgen.load_genotype_cache('{bad_schema}'); catch; e1 = 1; end; "
        f"e2 = 0; try; statgen.load_genotype_cache('{bad_offsets}'); catch; e2 = 1; end; "
        f"e3 = 0; try; statgen.load_genotype_cache('{unknown}', {{'2'}}); catch; e3 = 1; end; "
        f"e4 = 0; try; statgen.load_genotype_cache('{bad_zero_shards}'); catch ME; e4 = ~isempty(strfind(ME.message, 'n_shards must be positive')); end; "
        "fprintf('%d %d %d %d\\n', e1, e2, e3, e4);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1 1 1"]
