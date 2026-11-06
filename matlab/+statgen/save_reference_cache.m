function save_reference_cache(panel, path, varargin)
%SAVE_REFERENCE_CACHE Save a ReferencePanel to a MATLAB .mat cache.
%
%   statgen.save_reference_cache(panel, path)
%   statgen.save_reference_cache(panel, path, 'format', 'v7')
%
% String fields are stored as per-shard newline-delimited character payloads
% and are decoded lazily by cache-loaded ReferenceShard objects.
%
% Default format is 'v7.3' in MATLAB and 'v7' in Octave.
%
% See also statgen.ReferencePanel, statgen.load_reference_cache.
    path = char(path);
    save_arg = parse_args_(varargin{:});

    n_shards = numel(panel.shards);
    metadata.schema = 'reference_cache/0.1';
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);

    for i = 1:n_shards
        s = panel.shards{i};
        off = panel.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
    end

    bp = panel.bp;
    a1_hash64 = panel.a1_hash64;
    a2_hash64 = panel.a2_hash64;

    snp_text_by_shard = cell(n_shards, 1);
    a1_text_by_shard = cell(n_shards, 1);
    a2_text_by_shard = cell(n_shards, 1);
    for i = 1:n_shards
        [snp_text_by_shard{i}, a1_text_by_shard{i}, a2_text_by_shard{i}] = ...
            panel.shards{i}.cache_text_payloads();
    end

    payload = struct();
    payload.metadata = metadata;
    payload.bp = bp;
    payload.snp_text_by_shard = snp_text_by_shard;
    payload.a1_text_by_shard = a1_text_by_shard;
    payload.a2_text_by_shard = a2_text_by_shard;
    payload.a1_hash64 = a1_hash64;
    payload.a2_hash64 = a2_hash64;
    statgen.internal.save_cache_mat(path, payload, save_arg);
end

function save_arg = parse_args_(varargin)
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_reference_cache', statgen.internal.default_cache_mat_format(), {'v7', 'v7.3', 'v5'});

    if mod(numel(varargin), 2) ~= 0
        error('statgen:io', 'save_reference_cache options must be name-value pairs');
    end
    seen_format = false;
    for i = 1:2:numel(varargin)
        if ~(ischar(varargin{i}) || isstring(varargin{i}))
            error('statgen:io', 'save_reference_cache option names must be strings');
        end
        name = char(varargin{i});
        if strcmpi(name, 'format')
            if seen_format
                error('statgen:io', 'save_reference_cache format option specified more than once');
            end
            seen_format = true;
            [~, save_arg] = statgen.internal.parse_mat_format( ...
                'save_reference_cache', varargin{i + 1}, {'v7', 'v7.3', 'v5'});
        else
            error('statgen:io', 'save_reference_cache unknown option: %s', name);
        end
    end
end
