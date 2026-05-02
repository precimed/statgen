function panel = create_annotations(reference, annomat, annonames)
% Create AnnotationPanel from a reference-aligned binary matrix.
    names = ensure_names_(annonames);

    if issparse(annomat)
        A = sparse(annomat);
    else
        A = double(annomat);
    end
    if ndims(A) ~= 2
        error('statgen:annotations', 'annomat must be a 2D matrix');
    end

    n = reference.num_snp;
    k = numel(names);
    if size(A, 1) ~= n || size(A, 2) ~= k
        error('statgen:annotations', ...
            'annomat shape mismatch: expected (%d, %d), got (%d, %d)', ...
            n, k, size(A, 1), size(A, 2));
    end

    vals = nonzeros(sparse(A));
    bad = (vals ~= 0) & (vals ~= 1);
    if any(bad)
        v = vals(find(bad, 1, 'first'));
        error('statgen:annotations', 'annomat contains non-binary value: %g', v);
    end

    A = spones(sparse(A));

    shards = cell(numel(reference.shards), 1);
    for i = 1:numel(reference.shards)
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        shards{i} = statgen.AnnotationShard(s_ref.label, s_ref.checksum, A(ix, :));
    end

    panel = statgen.AnnotationPanel(shards, names);
end

function out = ensure_names_(names)
    if ischar(names) || isstring(names)
        names = cellstr(names(:));
    end
    if ~iscell(names)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end

    out = statgen.ensure_cell_col(names);
    if isempty(out)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end
    if any(cellfun('isempty', out))
        error('statgen:annotations', 'annotation names must not contain empty strings');
    end
    if numel(unique(out)) ~= numel(out)
        error('statgen:annotations', 'annotation names must be unique');
    end
end
