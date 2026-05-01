function selected = validate_requested_shards(requested, available, where)
% Validate shard subset request against canonical ordering and availability.
    selected = statgen.validate_requested_shards(requested, available, where);
end
