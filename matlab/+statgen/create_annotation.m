function panel = create_annotation(reference, annovec, annoname)
% Create single-column AnnotationPanel from a binary aligned vector.
    if nargin < 3
        error('statgen:arg', 'create_annotation requires reference, annovec, annoname');
    end

    name = char(annoname);
    if isempty(name)
        error('statgen:annotations', 'annoname must be non-empty');
    end

    n = reference.num_snp;
    vec = annovec;
    if ndims(vec) > 2 || (ismatrix(vec) && ~any(size(vec) == 1))
        error('statgen:annotations', 'annovec must be a vector with length %d', n);
    end
    vec = double(vec(:));
    if numel(vec) ~= n
        error('statgen:annotations', 'annovec length mismatch: expected %d, got %d', n, numel(vec));
    end

    bad = (vec ~= 0) & (vec ~= 1);
    if any(bad)
        i = find(bad, 1, 'first');
        error('statgen:annotations', 'annovec(%d) must be binary (0/1)', i);
    end

    panel = statgen.create_annotations(reference, sparse(vec), {name});
end
