classdef ReferenceShard
% Immutable in-memory representation of one .bim shard.
    properties (SetAccess = private)
        label      % chromosome label string
        num_snp    % scalar integer
        chr        % n×1 cell array of strings
        snp        % n×1 cell array of strings
        bp         % n×1 double vector
        a1         % n×1 cell array of strings
        a2         % n×1 cell array of strings
        checksum   % lowercase MD5 hex string
    end

    methods
        function obj = ReferenceShard(label, chr_vec, snp_vec, bp_vec, a1_vec, a2_vec)
            if nargin == 0, return; end
            obj.label   = char(label);
            obj.num_snp = numel(chr_vec);
            obj.chr     = chr_vec(:);
            obj.snp     = snp_vec(:);
            obj.bp      = double(bp_vec(:));
            obj.a1      = a1_vec(:);
            obj.a2      = a2_vec(:);
            % MD5 over 'chr:bp:a1:a2\n' lines in row order
            n = obj.num_snp;
            parts = cell(n, 1);
            for i = 1:n
                parts{i} = sprintf('%s:%d:%s:%s\n', ...
                    char(chr_vec{i}), round(bp_vec(i)), char(a1_vec{i}), char(a2_vec{i}));
            end
            obj.checksum = hash('md5', strjoin(parts, ''));
        end
    end
end
