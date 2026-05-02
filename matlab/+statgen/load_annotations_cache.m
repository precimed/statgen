function panel = load_annotations_cache(path, shards)
% Load AnnotationPanel from MATLAB .mat cache file with optional shard subset.
    if nargin < 2
        shards = [];
    end

    path = char(path);
    loaded = load(path, 'cache_meta', 'cache_shards');
    meta = loaded.cache_meta;
    shards_data = loaded.cache_shards;

    if ~strcmp(meta.schema, 'annotations_cache/0.1')
        error('statgen:cache', 'Unsupported annotations cache schema: %s', meta.schema);
    end

    labels = statgen.ensure_cell_col(meta.shard_labels);
    checksums = statgen.ensure_cell_col(meta.shard_checksums);
    if ~isfield(meta, 'n_shards') || ~isscalar(meta.n_shards) || meta.n_shards ~= numel(labels)
        error('statgen:cache', 'Invalid annotations cache: n_shards mismatch');
    end
    if numel(labels) ~= numel(checksums)
        error('statgen:cache', 'Invalid annotations cache: shard_labels and shard_checksums length mismatch');
    end
    if numel(shards_data) ~= numel(labels)
        error('statgen:cache', ...
            'Invalid annotations cache: cache_shards length does not match shard labels');
    end

    annonames = statgen.ensure_cell_col(meta.annonames);
    if isempty(annonames)
        error('statgen:cache', 'Invalid annotations cache: annonames must be non-empty');
    end
    if any(cellfun('isempty', annonames))
        error('statgen:cache', 'Invalid annotations cache: annonames must not contain empty names');
    end
    if numel(unique(annonames)) ~= numel(annonames)
        error('statgen:cache', 'Invalid annotations cache: annonames must be unique');
    end

    selected = statgen.validate_requested_shards(shards, labels, 'load_annotations_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        sd = shards_data(idx);
        out_shards{i} = statgen.AnnotationShard(labels{idx}, checksums{idx}, sd.annomat);
    end

    panel = statgen.AnnotationPanel(out_shards, annonames);
end
