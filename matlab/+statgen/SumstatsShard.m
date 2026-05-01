classdef SumstatsShard
% Immutable in-memory aligned sumstats vectors for one reference shard.
    properties (SetAccess = private)
        label
        checksum
        num_snp
        zvec
        nvec
        logpvec
        beta_vec
        se_vec
        eaf_vec
        info_vec
    end

    methods
        function obj = SumstatsShard(label, checksum, zvec, nvec, logpvec, beta_vec, se_vec, eaf_vec, info_vec)
            if nargin == 0, return; end
            obj.label = char(label);
            obj.checksum = char(checksum);
            obj.zvec = double(zvec(:));
            obj.nvec = double(nvec(:));
            obj.logpvec = double(logpvec(:));
            if nargin < 6 || isempty(beta_vec), obj.beta_vec = []; else, obj.beta_vec = double(beta_vec(:)); end
            if nargin < 7 || isempty(se_vec),   obj.se_vec = [];   else, obj.se_vec = double(se_vec(:));   end
            if nargin < 8 || isempty(eaf_vec),  obj.eaf_vec = [];  else, obj.eaf_vec = double(eaf_vec(:)); end
            if nargin < 9 || isempty(info_vec), obj.info_vec = []; else, obj.info_vec = double(info_vec(:)); end
            obj.num_snp = numel(obj.zvec);
        end
    end
end
