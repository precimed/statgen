classdef Sumstats
% Ordered collection of SumstatsShard objects with genome-wide accessors.
    properties (SetAccess = private)
        num_snp
        shard_offsets
        shards
    end
    properties (Dependent)
        zvec
        nvec
        logpvec
        beta_vec
        se_vec
        eaf_vec
        info_vec
    end

    methods
        function obj = Sumstats(shards_cell)
            if nargin == 0, return; end
            obj.shards = shards_cell;
            n_shards = numel(shards_cell);

            total = 0;
            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});

            for i = 1:n_shards
                s = shards_cell{i};
                n_i = s.num_snp;
                offsets(i).shard_label = s.label;
                offsets(i).start0 = total;
                offsets(i).stop0 = total + n_i;
                total = total + n_i;
            end

            obj.num_snp = total;
            obj.shard_offsets = offsets;
        end

        function out = get.zvec(obj)
            out = concat_required_(obj.shards, 'zvec');
        end

        function out = get.nvec(obj)
            out = concat_required_(obj.shards, 'nvec');
        end

        function out = get.logpvec(obj)
            out = concat_required_(obj.shards, 'logpvec');
        end

        function out = get.beta_vec(obj)
            out = concat_optional_(obj.shards, 'beta_vec');
        end

        function out = get.se_vec(obj)
            out = concat_optional_(obj.shards, 'se_vec');
        end

        function out = get.eaf_vec(obj)
            out = concat_optional_(obj.shards, 'eaf_vec');
        end

        function out = get.info_vec(obj)
            out = concat_optional_(obj.shards, 'info_vec');
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

function out = concat_required_(shards, field_name)
    if isempty(shards)
        out = [];
        return
    end
    vals = cell(numel(shards), 1);
    for i = 1:numel(shards)
        vals{i} = shards{i}.(field_name);
    end
    out = vertcat(vals{:});
end

function out = concat_optional_(shards, field_name)
    if isempty(shards)
        out = [];
        return
    end
    vals = cell(numel(shards), 1);
    for i = 1:numel(shards)
        v = shards{i}.(field_name);
        if isempty(v)
            out = [];
            return
        end
        vals{i} = v;
    end
    out = vertcat(vals{:});
end
