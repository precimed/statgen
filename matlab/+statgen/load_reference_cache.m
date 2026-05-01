function panel = load_reference_cache(path, shards)
% Load a ReferencePanel from a MATLAB binary .mat cache file.
%
%   panel = statgen.load_reference_cache(path)
%   panel = statgen.load_reference_cache(path, shards)
    if nargin < 2
        shards = [];
    end

    path = char(path);
    loaded = load(path, 'cache_meta', 'cache_shards');

    meta   = loaded.cache_meta;
    shards_data = loaded.cache_shards;

    if ~strcmp(meta.schema, 'reference_cache/0.1')
        error('statgen:cache', 'Unsupported reference cache schema: %s', meta.schema);
    end

    labels = statgen.ensure_cell_col(meta.shard_labels);
    checksums = statgen.ensure_cell_col(meta.shard_checksums);
    if numel(labels) ~= numel(checksums)
        error('statgen:cache', ...
            'Invalid reference cache: shard_labels and shard_checksums length mismatch');
    end
    if numel(shards_data) ~= numel(labels)
        error('statgen:cache', ...
            'Invalid reference cache: cache_shards length does not match shard_labels');
    end

    selected = statgen.validate_requested_shards(shards, labels, 'load_reference_cache');
    shard_objs = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');

        sd   = shards_data(idx);
        chr_v = statgen.ensure_cell_col(sd.chr);
        snp_v = statgen.ensure_cell_col(sd.snp);
        a1_v  = statgen.ensure_cell_col(sd.a1);
        a2_v  = statgen.ensure_cell_col(sd.a2);
        shard_obj = statgen.ReferenceShard( ...
            labels{idx}, chr_v, snp_v, sd.bp(:), a1_v, a2_v);
        if ~strcmp(shard_obj.checksum, checksums{idx})
            error('statgen:cache', ...
                'Invalid reference cache: shard %s checksum mismatch', labels{idx});
        end
        shard_objs{i} = shard_obj;
    end

    panel = statgen.ReferencePanel(shard_objs);
end
