function selected = validate_requested_shards(requested, available, where)
% Validate shard subset request against canonical ordering and availability.
    if nargin < 1 || isempty(requested)
        selected = available;
        return
    end

    if ischar(requested)
        requested = {requested};
    elseif isnumeric(requested)
        error('statgen:shards', '%s: shards must be a non-empty list of unique canonical contig labels', where);
    end

    requested = requested(:)';
    if isempty(requested)
        error('statgen:shards', '%s: shards must be a non-empty list of unique canonical contig labels', where);
    end

    canonical = statgen.canonical_labels();
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
