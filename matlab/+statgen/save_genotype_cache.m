function save_genotype_cache(panel, path, varargin)
%SAVE_GENOTYPE_CACHE Save genotype metadata to a MATLAB .mat cache.
%
%   statgen.save_genotype_cache(genotype, path)
%   statgen.save_genotype_cache(genotype, path, 'format', format)
%
% Saves a GenotypePanel metadata cache, including source BED paths. Dense
% genotype calls are not stored in the cache.
%
% Default format is 'v7.3' in MATLAB and 'v7' in Octave.
%
% See also statgen.GenotypePanel, statgen.load_genotype_cache.
    path = char(path);
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_genotype_cache', statgen.internal.default_cache_mat_format(), {'v7', 'v7.3', 'v5'}, varargin{:});

    n_shards = numel(panel.shards);
    metadata.schema = 'genotype_cache/0.2';
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);
    metadata.source_layout = panel.source_layout;
    metadata.bed_paths = cell(n_shards, 1);
    metadata.bed_file_sizes = zeros(n_shards, 1);
    metadata.source_num_snp = zeros(n_shards, 1);
    metadata.source_num_sample = zeros(n_shards, 1);
    metadata.num_sample = panel.num_sample;

    subject_present = false(panel.num_sample, n_shards);
    source_subject_row0 = -ones(panel.num_sample, n_shards);
    for i = 1:n_shards
        s = panel.shards{i};
        off = panel.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.reference_checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
        metadata.bed_paths{i} = s.bed_path;
        metadata.bed_file_sizes(i) = s.bed_file_size;
        metadata.source_num_snp(i) = s.source_num_snp;
        metadata.source_num_sample(i) = s.source_num_sample;
        subject_present(:, i) = s.subject_present;
        source_subject_row0(:, i) = s.source_subject_row0;
    end

    is_present = panel.is_present;
    is_allele_swapped = panel.is_allele_swapped;
    ploidy_male = panel.ploidy_male;
    ploidy_female = panel.ploidy_female;
    source_row0 = panel.source_row0;
    fid = panel.fid;
    iid = panel.iid;
    father_id = panel.father_id;
    mother_id = panel.mother_id;
    sex = panel.sex;
    is_male = panel.is_male;
    is_female = panel.is_female;

    payload = struct();
    payload.metadata = metadata;
    payload.is_present = is_present;
    payload.is_allele_swapped = is_allele_swapped;
    payload.ploidy_male = ploidy_male;
    payload.ploidy_female = ploidy_female;
    payload.source_row0 = source_row0;
    payload.subject_present = subject_present;
    payload.source_subject_row0 = source_subject_row0;
    payload.fid = fid;
    payload.iid = iid;
    payload.father_id = father_id;
    payload.mother_id = mother_id;
    payload.sex = sex;
    payload.is_male = is_male;
    payload.is_female = is_female;
    statgen.internal.save_cache_mat(path, payload, save_arg);
end
