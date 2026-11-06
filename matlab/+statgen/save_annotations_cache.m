function save_annotations_cache(panel, path, varargin)
%SAVE_ANNOTATIONS_CACHE Save annotations to a MATLAB .mat cache.
%
%   statgen.save_annotations_cache(annotations, path)
%   statgen.save_annotations_cache(annotations, path, 'format', format)
%
% Saves an AnnotationPanel for faster reload in MATLAB or Octave.
%
% Default format is 'v7.3' in MATLAB and 'v7' in Octave.
%
% See also statgen.AnnotationPanel, statgen.load_annotations_cache.
    path = char(path);
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_annotations_cache', statgen.internal.default_cache_mat_format(), {'v7', 'v7.3', 'v5'}, varargin{:});

    n_shards = numel(panel.shards);
    metadata.schema = 'annotations_cache/0.2';
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);

    for i = 1:n_shards
        s = panel.shards{i};
        off = panel.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.reference_checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
    end

    annomat = panel.annomat;
    annonames = panel.annonames;
    is_binary = logical(panel.is_binary(:));
    annotation_metadata = statgen.internal.ensure_cell_col(panel.annotation_metadata);

    payload = struct();
    payload.metadata = metadata;
    payload.annomat = annomat;
    payload.annonames = annonames;
    payload.is_binary = is_binary;
    payload.annotation_metadata = annotation_metadata;
    statgen.internal.save_cache_mat(path, payload, save_arg);
end
