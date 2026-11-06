function panel = load_genotype_cache(path, shards)
%LOAD_GENOTYPE_CACHE Load cached genotype metadata.
%
%   genotype = statgen.load_genotype_cache(path)
%   genotype = statgen.load_genotype_cache(path, shards)
%
% Loads one MATLAB .mat GenotypePanel cache file produced by
% statgen.save_genotype_cache. The cache stores metadata and source BED paths,
% not dense genotype calls. Optional shards return a logical subset of the
% cached shard labels. Fetching genotypes still requires the original BED files,
% or an explicit BED path override in fetch_genotypes.
%
% See also statgen.GenotypePanel, statgen.load_genotype,
% statgen.save_genotype_cache.
    if nargin < 2
        shards = [];
    end
    path = char(path);
    loaded = load(path);
    require_fields_(loaded, {'metadata', 'is_present', 'ploidy_male', 'ploidy_female', ...
        'source_row0', 'subject_present', 'source_subject_row0', 'fid', 'iid', ...
        'father_id', 'mother_id', 'sex', 'is_male', 'is_female'});

    meta = loaded.metadata;
    is_present = logical(loaded.is_present(:));
    ploidy_male = double(loaded.ploidy_male(:));
    ploidy_female = double(loaded.ploidy_female(:));
    source_row0 = double(loaded.source_row0(:));
    num_snp = numel(is_present);
    if numel(ploidy_male) ~= num_snp || numel(ploidy_female) ~= num_snp || numel(source_row0) ~= num_snp
        error('statgen:cache', 'Invalid genotype cache: SNP-axis vector lengths mismatch');
    end

    [labels, checksums, start0, stop0, source_layout, bed_paths, bed_file_sizes, source_num_snp, source_num_sample, num_sample] = validate_metadata_(meta, num_snp);

    if strcmp(meta.schema, 'genotype_cache/0.2')
        require_fields_(loaded, {'is_allele_swapped'});
        flags = loaded.is_allele_swapped;
        if ~islogical(flags) || ~isvector(flags) || numel(flags) ~= num_snp
            error('statgen:cache', 'Invalid genotype cache: is_allele_swapped must be a logical SNP-axis vector');
        end
        is_allele_swapped = flags(:);
        if any(is_allele_swapped & ~is_present)
            error('statgen:cache', 'Invalid genotype cache: is_allele_swapped requires is_present');
        end
    else
        is_allele_swapped = false(num_snp, 1);
    end

    fid = statgen.internal.ensure_cell_col(loaded.fid);
    iid = statgen.internal.ensure_cell_col(loaded.iid);
    father_id = statgen.internal.ensure_cell_col(loaded.father_id);
    mother_id = statgen.internal.ensure_cell_col(loaded.mother_id);
    sex = double(loaded.sex(:));
    is_male = logical(loaded.is_male(:));
    is_female = logical(loaded.is_female(:));
    if numel(fid) ~= num_sample || numel(iid) ~= num_sample || numel(father_id) ~= num_sample || numel(mother_id) ~= num_sample || numel(sex) ~= num_sample || numel(is_male) ~= num_sample || numel(is_female) ~= num_sample
        error('statgen:cache', 'Invalid genotype cache: FAM vector length mismatch');
    end
    if any(sex ~= floor(sex) | ~ismember(sex, [0; 1; 2]))
        error('statgen:cache', 'Invalid genotype cache: sex values must be 0, 1, or 2');
    end
    if ~isequal(is_male, sex == 1)
        error('statgen:cache', 'Invalid genotype cache: is_male does not match sex');
    end
    if ~isequal(is_female, sex == 2)
        error('statgen:cache', 'Invalid genotype cache: is_female does not match sex');
    end
    if ~isequal(size(loaded.subject_present), [num_sample, numel(labels)])
        error('statgen:cache', 'Invalid genotype cache: subject_present shape mismatch');
    end
    if ~isequal(size(loaded.source_subject_row0), [num_sample, numel(labels)])
        error('statgen:cache', 'Invalid genotype cache: source_subject_row0 shape mismatch');
    end
    subject_present = logical(loaded.subject_present);
    source_subject_row0 = double(loaded.source_subject_row0);

    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_genotype_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);
        out_shards{i} = statgen.GenotypeShard( ...
            labels{idx}, bed_paths{idx}, bed_file_sizes(idx), ...
            source_num_snp(idx), source_num_sample(idx), ...
            source_row0(ix), subject_present(:, idx), source_subject_row0(:, idx), ...
            is_present(ix), is_allele_swapped(ix), ploidy_male(ix), ploidy_female(ix), checksums{idx});
    end

    panel = statgen.GenotypePanel(out_shards, fid, iid, father_id, mother_id, sex, source_layout);
end

function require_fields_(loaded, fields)
    for i = 1:numel(fields)
        if ~isfield(loaded, fields{i})
            error('statgen:cache', 'Invalid genotype cache: missing field %s', fields{i});
        end
    end
end

function [labels, checksums, start0, stop0, source_layout, bed_paths, bed_file_sizes, source_num_snp, source_num_sample, num_sample] = validate_metadata_(meta, num_snp)
    if ~isfield(meta, 'schema') || ~any(strcmp(meta.schema, {'genotype_cache/0.1', 'genotype_cache/0.2'}))
        if isfield(meta, 'schema')
            schema = meta.schema;
        else
            schema = '<missing>';
        end
        error('statgen:cache', 'Unsupported genotype cache schema: %s', schema);
    end
    need = {'n_shards', 'shard_labels', 'shard_checksums', 'shard_start0', ...
        'shard_stop0', 'source_layout', 'bed_paths', 'bed_file_sizes', ...
        'source_num_snp', 'source_num_sample', 'num_sample'};
    for i = 1:numel(need)
        if ~isfield(meta, need{i})
            error('statgen:cache', 'Invalid genotype cache: missing metadata field %s', need{i});
        end
    end

    [labels, checksums, start0, stop0, n_shards] = statgen.internal.validate_cache_shard_metadata( ...
        meta, meta.schema, 'genotype', num_snp, true);
    bed_paths = statgen.internal.ensure_cell_col(meta.bed_paths);
    bed_file_sizes = double(meta.bed_file_sizes(:));
    source_num_snp = double(meta.source_num_snp(:));
    source_num_sample = double(meta.source_num_sample(:));
    source_layout = char(meta.source_layout);
    num_sample = double(meta.num_sample);

    if numel(bed_paths) ~= n_shards || numel(bed_file_sizes) ~= n_shards || ...
            numel(source_num_snp) ~= n_shards || numel(source_num_sample) ~= n_shards
        error('statgen:cache', 'Invalid genotype cache: shard metadata length mismatch');
    end
    statgen.internal.validate_requested_shards(labels, labels, 'load_genotype_cache');
    if ~strcmp(source_layout, 'non_sharded') && ~strcmp(source_layout, 'sharded')
        error('statgen:cache', 'Invalid genotype cache: source_layout must be ''non_sharded'' or ''sharded''');
    end
    if ~isscalar(num_sample) || num_sample < 0 || num_sample ~= floor(num_sample)
        error('statgen:cache', 'Invalid genotype cache: num_sample must be a non-negative integer');
    end

    for i = 1:n_shards
        expected = statgen.internal.bfile_expected_bed_size(source_num_sample(i), source_num_snp(i));
        if bed_file_sizes(i) ~= expected
            error('statgen:cache', 'Invalid genotype cache: bed_file_size mismatch for shard %s', labels{i});
        end
    end
end
