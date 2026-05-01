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

    labels = ensure_cell_col_(meta.shard_labels);
    checksums = ensure_cell_col_(meta.shard_checksums);
    if numel(labels) ~= numel(checksums)
        error('statgen:cache', ...
            'Invalid reference cache: shard_labels and shard_checksums length mismatch');
    end
    if numel(shards_data) ~= numel(labels)
        error('statgen:cache', ...
            'Invalid reference cache: cache_shards length does not match shard_labels');
    end

    selected = validate_requested_shards_(shards, labels, 'load_reference_cache');
    shard_objs = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');

        sd   = shards_data(idx);
        chr_v = ensure_cell_col_(sd.chr);
        snp_v = ensure_cell_col_(sd.snp);
        a1_v  = ensure_cell_col_(sd.a1);
        a2_v  = ensure_cell_col_(sd.a2);
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

function out = ensure_cell_col_(x)
    if iscell(x)
        out = x(:);
    elseif ischar(x)
        out = strtrim(cellstr(x));
        out = out(:);
    else
        out = strtrim(cellstr(x));
        out = out(:);
    end
end

function selected = validate_requested_shards_(requested, available, where)
    if nargin < 1 || isempty(requested)
        selected = available;
        return
    end

    if ischar(requested)
        requested = {requested};
    elseif isnumeric(requested)
        error('statgen:shards', ...
            '%s: shards must be a non-empty list of unique canonical contig labels', ...
            where);
    end

    requested = requested(:)';
    if isempty(requested)
        error('statgen:shards', ...
            '%s: shards must be a non-empty list of unique canonical contig labels', ...
            where);
    end

    canonical = [arrayfun(@num2str, 1:22, 'UniformOutput', false), {'X'}];
    canonical_idx = zeros(1, numel(requested));
    seen = {};
    for i = 1:numel(requested)
        label = char(requested{i});
        idx = find(strcmp(canonical, label), 1, 'first');
        if isempty(idx)
            error('statgen:shards', ...
                '%s: unsupported shard label %s; expected canonical labels 1-22 or X', ...
                where, label);
        end
        if any(strcmp(seen, label))
            error('statgen:shards', '%s: duplicate shard label %s in shards', where, label);
        end
        if ~any(strcmp(available, label))
            error('statgen:shards', '%s: requested shard %s is not present', where, label);
        end
        seen{end+1} = label; %#ok<AGROW>
        canonical_idx(i) = idx;
    end

    if any(diff(canonical_idx) <= 0)
        error('statgen:shards', '%s: shards must be in canonical subsequence order', where);
    end

    selected = requested;
end
