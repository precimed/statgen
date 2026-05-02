function save_annotations_cache(panel, path)
% Save AnnotationPanel to MATLAB binary .mat cache file.
    path = char(path);
    n_shards = numel(panel.shards);

    cache_meta.schema = 'annotations_cache/0.1';
    cache_meta.n_shards = n_shards;
    cache_meta.shard_labels = cell(n_shards, 1);
    cache_meta.shard_checksums = cell(n_shards, 1);
    cache_meta.annonames = panel.annonames;

    cache_shards = struct('annomat', {});
    for i = 1:n_shards
        s = panel.shards{i};
        cache_meta.shard_labels{i} = s.label;
        cache_meta.shard_checksums{i} = s.checksum;
        cache_shards(i).annomat = s.annomat;
    end

    save(path, 'cache_meta', 'cache_shards', '-v7');
end
