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
            bp_str = cellstr(num2str(round(obj.bp), '%d'));
            parts = strcat(obj.chr, {':'}, bp_str, {':'}, obj.a1, {':'}, obj.a2, {sprintf('\n')});
            text_payload = [parts{:}];
            computed = md5_hex_(text_payload);
            obj.checksum = computed;
        end
    end
end

function out = md5_hex_(text_payload)
    try
        out = lower(hash('md5', text_payload));
        return
    catch
        % MATLAB path: use Java MessageDigest when hash(...) is unavailable.
    end

    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(text_payload));
    d = typecast(md.digest(), 'uint8');
    out = lower(reshape(dec2hex(d)', 1, []));
end
