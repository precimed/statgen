function save_sumstats_cache(sumstats, path, varargin)
%SAVE_SUMSTATS_CACHE Save summary statistics to a MATLAB .mat cache.
%
%   statgen.save_sumstats_cache(sumstats, path)
%   statgen.save_sumstats_cache(sumstats, path, 'format', format)
%
% Saves a Sumstats object for faster reload in MATLAB or Octave.
%
% Default format is 'v7.3' in MATLAB and 'v7' in Octave.
%
% See also statgen.Sumstats, statgen.load_sumstats_cache.
    path = char(path);
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_sumstats_cache', statgen.internal.default_cache_mat_format(), {'v7', 'v7.3', 'v5'}, varargin{:});

    n_shards = numel(sumstats.shards);
    metadata.schema = 'sumstats_cache/0.1';
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);
    metadata.has_z = ~isempty(sumstats.zvec);
    metadata.has_n = ~isempty(sumstats.nvec);
    metadata.has_beta = ~isempty(sumstats.beta_vec);
    metadata.has_se = ~isempty(sumstats.se_vec);
    metadata.has_eaf = ~isempty(sumstats.eaf_vec);
    metadata.has_info = ~isempty(sumstats.info_vec);

    for i = 1:n_shards
        s = sumstats.shards{i};
        off = sumstats.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.reference_checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
    end

    logpvec = sumstats.logpvec;
    zvec = sumstats.zvec;
    nvec = sumstats.nvec;
    beta_vec = sumstats.beta_vec;
    se_vec = sumstats.se_vec;
    eaf_vec = sumstats.eaf_vec;
    info_vec = sumstats.info_vec;

    payload = struct();
    payload.metadata = metadata;
    payload.zvec = zvec;
    payload.nvec = nvec;
    payload.logpvec = logpvec;
    payload.beta_vec = beta_vec;
    payload.se_vec = se_vec;
    payload.eaf_vec = eaf_vec;
    payload.info_vec = info_vec;
    statgen.internal.save_cache_mat(path, payload, save_arg);
end
