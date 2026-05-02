classdef AnnotationShard
% Immutable annotation matrix for one reference shard.
    properties (SetAccess = private)
        label
        checksum
        num_snp
        num_annot
        annomat
    end

    methods
        function obj = AnnotationShard(label, checksum, annomat)
            if nargin == 0, return; end
            obj.label = char(label);
            obj.checksum = char(checksum);
            obj.annomat = ensure_sparse_binary_(annomat);
            obj.num_snp = size(obj.annomat, 1);
            obj.num_annot = size(obj.annomat, 2);
        end
    end
end

function out = ensure_sparse_binary_(x)
    if issparse(x)
        out = sparse(x);
    else
        out = sparse(double(x));
    end
    vals = nonzeros(out);
    if any((vals ~= 0) & (vals ~= 1))
        bad = vals(find((vals ~= 0) & (vals ~= 1), 1, 'first'));
        error('statgen:annotations', 'annomat contains non-binary value: %g', bad);
    end
    out = spones(out);
end
