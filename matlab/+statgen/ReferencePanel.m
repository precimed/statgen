classdef ReferencePanel
% Ordered collection of ReferenceShard objects with genome-wide accessors.
    properties (SetAccess = private)
        num_snp        % total SNP count (scalar)
        shard_offsets  % struct array: shard_label, start0, stop0 (zero-based half-open)
        shards         % cell array of ReferenceShard objects
    end
    properties (Dependent)
        chr            % num_snp×1 cell array of chromosome labels
        snp            % num_snp×1 cell array of SNP identifiers
        bp             % num_snp×1 double vector of base-pair positions
        a1             % num_snp×1 cell array
        a2             % num_snp×1 cell array
    end

    methods
        function obj = ReferencePanel(shards_cell)
            if nargin == 0, return; end
            obj.shards = shards_cell;
            n_shards = numel(shards_cell);

            total   = 0;
            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});

            for i = 1:n_shards
                s = shards_cell{i};
                n_i = s.num_snp;
                offsets(i).shard_label = s.label;
                offsets(i).start0      = total;
                offsets(i).stop0       = total + n_i;
                total   = total + n_i;
            end

            obj.num_snp       = total;
            obj.shard_offsets = offsets;
        end

        function out = get.chr(obj)
            if isempty(obj.shards), out = {}; return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.chr;
            end
            out = vertcat(vals{:});
        end

        function out = get.snp(obj)
            if isempty(obj.shards), out = {}; return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.snp;
            end
            out = vertcat(vals{:});
        end

        function out = get.bp(obj)
            if isempty(obj.shards), out = []; return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.bp;
            end
            out = vertcat(vals{:});
        end

        function out = get.a1(obj)
            if isempty(obj.shards), out = {}; return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.a1;
            end
            out = vertcat(vals{:});
        end

        function out = get.a2(obj)
            if isempty(obj.shards), out = {}; return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.a2;
            end
            out = vertcat(vals{:});
        end

        function out = select_shards(obj, shards)
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.validate_requested_shards( ...
                shards, available, 'ReferencePanel.select_shards');

            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.ReferencePanel(out_shards);
        end

        function ok = is_object_compatible(obj, other)
            ok = true;

            try
                other_shards = other.shards;
            catch
                statgen.ReferencePanel.compat_warn_( ...
                    'statgen: is_object_compatible: object has no shards property');
                ok = false; return;
            end

            if numel(other_shards) ~= numel(obj.shards)
                statgen.ReferencePanel.compat_warn_( ...
                    'statgen: is_object_compatible: shard count mismatch: ref=%d, obj=%d', ...
                    numel(obj.shards), numel(other_shards));
                ok = false; return;
            end

            for i = 1:numel(obj.shards)
                rs = obj.shards{i};
                try
                    os = other_shards{i};
                catch
                    os = other_shards(i);
                end

                try
                    os_label = os.label;
                catch
                    statgen.ReferencePanel.compat_warn_( ...
                        'statgen: is_object_compatible: shard %s: object shard has no label', ...
                        rs.label);
                    ok = false; continue;
                end

                if ~strcmp(os_label, rs.label)
                    statgen.ReferencePanel.compat_warn_( ...
                        'statgen: is_object_compatible: shard label mismatch: ref=%s, obj=%s', ...
                        rs.label, os_label);
                    ok = false; continue;
                end

                try
                    os_num = os.num_snp;
                catch
                    statgen.ReferencePanel.compat_warn_( ...
                        'statgen: is_object_compatible: shard %s: object shard has no num_snp', ...
                        rs.label);
                    ok = false; continue;
                end

                if os_num ~= rs.num_snp
                    statgen.ReferencePanel.compat_warn_( ...
                        'statgen: is_object_compatible: shard %s: row count mismatch: ref=%d, obj=%d', ...
                        rs.label, rs.num_snp, os_num);
                    ok = false; continue;
                end

                try
                    os_chk = os.checksum;
                    if ~strcmp(os_chk, rs.checksum)
                        statgen.ReferencePanel.compat_warn_( ...
                            'statgen: is_object_compatible: shard %s: checksum mismatch', ...
                            rs.label);
                        ok = false;
                    end
                catch
                    % no checksum field — skip check
                end
            end
        end
    end

    methods (Static, Access = private)
        function compat_warn_(fmt, varargin)
            msg = sprintf(fmt, varargin{:});
            warning('statgen:compat', '%s', msg);
        end
    end
end
