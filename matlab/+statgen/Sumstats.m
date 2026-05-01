classdef Sumstats
% Ordered collection of SumstatsShard objects with genome-wide accessors.
    properties (SetAccess = private)
        num_snp
        zvec
        nvec
        logpvec
        beta_vec
        se_vec
        eaf_vec
        info_vec
        shard_offsets
        shards
    end

    methods
        function obj = Sumstats(shards_cell)
            if nargin == 0, return; end
            obj.shards = shards_cell;
            n_shards = numel(shards_cell);

            total = 0;
            z_all = [];
            n_all = [];
            logp_all = [];
            beta_all = [];
            se_all = [];
            eaf_all = [];
            info_all = [];
            has_beta = true;
            has_se = true;
            has_eaf = true;
            has_info = true;

            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});

            for i = 1:n_shards
                s = shards_cell{i};
                n_i = s.num_snp;
                offsets(i).shard_label = s.label;
                offsets(i).start0 = total;
                offsets(i).stop0 = total + n_i;

                z_all = [z_all; s.zvec];
                n_all = [n_all; s.nvec];
                logp_all = [logp_all; s.logpvec];

                if isempty(s.beta_vec), has_beta = false; else, beta_all = [beta_all; s.beta_vec]; end
                if isempty(s.se_vec), has_se = false; else, se_all = [se_all; s.se_vec]; end
                if isempty(s.eaf_vec), has_eaf = false; else, eaf_all = [eaf_all; s.eaf_vec]; end
                if isempty(s.info_vec), has_info = false; else, info_all = [info_all; s.info_vec]; end

                total = total + n_i;
            end

            obj.num_snp = total;
            obj.zvec = z_all;
            obj.nvec = n_all;
            obj.logpvec = logp_all;
            if has_beta, obj.beta_vec = beta_all; else, obj.beta_vec = []; end
            if has_se, obj.se_vec = se_all; else, obj.se_vec = []; end
            if has_eaf, obj.eaf_vec = eaf_all; else, obj.eaf_vec = []; end
            if has_info, obj.info_vec = info_all; else, obj.info_vec = []; end
            obj.shard_offsets = offsets;
        end

        function out = select_shards(obj, shards)
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.validate_requested_shards(shards, available, 'Sumstats.select_shards');
            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.Sumstats(out_shards);
        end
    end
end
