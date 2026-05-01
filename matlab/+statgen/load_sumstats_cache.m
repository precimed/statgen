function sumstats = load_sumstats_cache(path, shards)
% Load Sumstats from MATLAB .mat cache file with optional shard subsetting.
    if nargin < 2
        shards = [];
    end
    path = char(path);
    loaded = load(path, 'cache_meta', 'cache_shards');
    meta = loaded.cache_meta;
    shards_data = loaded.cache_shards;

    if ~strcmp(meta.schema, 'sumstats_cache/0.1')
        error('statgen:cache', 'Unsupported sumstats cache schema: %s', meta.schema);
    end

    labels = statgen.ensure_cell_col(meta.shard_labels);
    checksums = statgen.ensure_cell_col(meta.shard_checksums);
    if numel(labels) ~= numel(checksums)
        error('statgen:cache', 'Invalid sumstats cache: shard_labels and shard_checksums length mismatch');
    end
    if numel(shards_data) ~= numel(labels)
        error('statgen:cache', 'Invalid sumstats cache: cache_shards length does not match shard labels');
    end

    selected = statgen.validate_requested_shards(shards, labels, 'load_sumstats_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        sd = shards_data(idx);
        out_shards{i} = statgen.SumstatsShard( ...
            labels{idx}, checksums{idx}, ...
            sd.zvec(:), sd.nvec(:), sd.logpvec(:), ...
            pick_opt_(sd, 'beta_vec', meta.has_beta), ...
            pick_opt_(sd, 'se_vec', meta.has_se), ...
            pick_opt_(sd, 'eaf_vec', meta.has_eaf), ...
            pick_opt_(sd, 'info_vec', meta.has_info));
    end
    sumstats = statgen.Sumstats(out_shards);
end

function out = pick_opt_(sd, name, has_field)
    if has_field
        out = sd.(name)(:);
    else
        out = [];
    end
end
