function save_reference_cache(panel, path)
% Save a ReferencePanel to a MATLAB binary .mat cache file.
% SNP-axis numeric vectors (bp) are stored as native binary arrays.
% String fields (chr, snp, a1, a2) are stored as cell arrays of chars.
    path = char(path);
    n_shards = numel(panel.shards);

    % Small metadata: schema, labels, checksums
    cache_meta.schema   = 'reference_cache/0.1';
    cache_meta.n_shards = n_shards;
    shard_labels    = cell(n_shards, 1);
    shard_checksums = cell(n_shards, 1);
    for i = 1:n_shards
        shard_labels{i}    = panel.shards{i}.label;
        shard_checksums{i} = panel.shards{i}.checksum;
    end
    cache_meta.shard_labels    = shard_labels;
    cache_meta.shard_checksums = shard_checksums;

    % Per-shard arrays: numeric vectors stored as native binary in .mat
    cache_shards = struct('chr', {}, 'snp', {}, 'bp', {}, 'a1', {}, 'a2', {});
    for i = 1:n_shards
        s = panel.shards{i};
        cache_shards(i).chr = s.chr;  % cell array of strings
        cache_shards(i).snp = s.snp;  % cell array of strings
        cache_shards(i).bp  = s.bp;   % double column vector (binary)
        cache_shards(i).a1  = s.a1;   % cell array of strings
        cache_shards(i).a2  = s.a2;   % cell array of strings
    end

    save(path, 'cache_meta', 'cache_shards', '-v7');
end
