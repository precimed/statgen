function save_sumstats_cache(sumstats, path)
% Save Sumstats to MATLAB binary .mat cache file.
    path = char(path);
    n_shards = numel(sumstats.shards);

    cache_meta.schema = 'sumstats_cache/0.1';
    cache_meta.n_shards = n_shards;
    cache_meta.shard_labels = cell(n_shards, 1);
    cache_meta.shard_checksums = cell(n_shards, 1);
    cache_meta.has_beta = ~isempty(sumstats.beta_vec);
    cache_meta.has_se = ~isempty(sumstats.se_vec);
    cache_meta.has_eaf = ~isempty(sumstats.eaf_vec);
    cache_meta.has_info = ~isempty(sumstats.info_vec);

    cache_shards = struct('zvec', {}, 'nvec', {}, 'logpvec', {}, ...
        'beta_vec', {}, 'se_vec', {}, 'eaf_vec', {}, 'info_vec', {});

    for i = 1:n_shards
        s = sumstats.shards{i};
        cache_meta.shard_labels{i} = s.label;
        cache_meta.shard_checksums{i} = s.checksum;
        cache_shards(i).zvec = s.zvec;
        cache_shards(i).nvec = s.nvec;
        cache_shards(i).logpvec = s.logpvec;
        cache_shards(i).beta_vec = s.beta_vec;
        cache_shards(i).se_vec = s.se_vec;
        cache_shards(i).eaf_vec = s.eaf_vec;
        cache_shards(i).info_vec = s.info_vec;
    end

    save(path, 'cache_meta', 'cache_shards', '-v7');
end
