function panel = load_reference_cache(path)
% Load a ReferencePanel from a MATLAB binary .mat cache file.
    path = char(path);
    loaded = load(path, 'cache_meta', 'cache_shards');

    meta   = loaded.cache_meta;
    shards = loaded.cache_shards;

    if ~strcmp(meta.schema, 'reference_cache/0.1')
        error('statgen:cache', 'Unsupported reference cache schema: %s', meta.schema);
    end

    shard_objs = cell(meta.n_shards, 1);
    for i = 1:meta.n_shards
        sd   = shards(i);
        chr_v = ensure_cell_col_(sd.chr);
        snp_v = ensure_cell_col_(sd.snp);
        a1_v  = ensure_cell_col_(sd.a1);
        a2_v  = ensure_cell_col_(sd.a2);
        shard_objs{i} = statgen.ReferenceShard( ...
            meta.shard_labels{i}, chr_v, snp_v, sd.bp(:), a1_v, a2_v);
    end

    panel = statgen.ReferencePanel(shard_objs);
end

function out = ensure_cell_col_(x)
    if iscell(x)
        out = x(:);
    elseif ischar(x)
        out = {x};
    else
        n = size(x, 1);
        out = cell(n, 1);
        for i = 1:n
            out{i} = strtrim(x(i, :));
        end
    end
end
