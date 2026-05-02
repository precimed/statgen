classdef AnnotationPanel
% Ordered collection of AnnotationShard objects with genome-wide accessors.
    properties (SetAccess = private)
        num_snp
        num_annot
        annonames
        shard_offsets
        shards
    end
    properties (Dependent)
        annomat
    end

    methods
        function obj = AnnotationPanel(shards_cell, annonames)
            if nargin == 0, return; end
            obj.shards = shards_cell;
            obj.annonames = ensure_names_(annonames);
            obj.num_annot = numel(obj.annonames);

            n_shards = numel(shards_cell);
            total = 0;
            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});

            for i = 1:n_shards
                s = shards_cell{i};
                if s.num_annot ~= obj.num_annot
                    error('statgen:annotations', 'all shards must share identical annotation columns');
                end
                n_i = s.num_snp;
                offsets(i).shard_label = s.label;
                offsets(i).start0 = total;
                offsets(i).stop0 = total + n_i;
                total = total + n_i;
            end

            obj.num_snp = total;
            obj.shard_offsets = offsets;
        end

        function out = get.annomat(obj)
            if isempty(obj.shards)
                out = sparse([], [], [], 0, obj.num_annot);
                return
            end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.annomat;
            end
            out = vertcat(vals{:});
        end

        function out = select_shards(obj, shards)
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.validate_requested_shards(shards, available, 'AnnotationPanel.select_shards');
            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.AnnotationPanel(out_shards, obj.annonames);
        end

        function out = select_annotations(obj, names)
            req = ensure_names_(names);
            [tf, idx] = ismember(req, obj.annonames);
            if ~all(tf)
                missing = req(~tf);
                error('statgen:annotations', 'unknown annotation name(s): %s', strjoin(missing, ', '));
            end

            out_shards = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                s = obj.shards{i};
                out_shards{i} = statgen.AnnotationShard(s.label, s.checksum, s.annomat(:, idx));
            end
            out = statgen.AnnotationPanel(out_shards, req);
        end

        function out = union_annotations(obj, other, mode)
            if nargin < 3 || isempty(mode)
                mode = 'by_name';
            end
            mode = char(mode);
            if ~strcmp(mode, 'by_name')
                error('statgen:annotations', 'union_annotations supports only mode=''by_name''');
            end

            try
                rhs_names = ensure_names_(other.annonames);
                rhs_shards = other.shards;
            catch
                error('statgen:annotations', 'union_annotations requires another AnnotationPanel-like object');
            end

            overlap = intersect(obj.annonames, rhs_names, 'stable');
            if ~isempty(overlap)
                error('statgen:annotations', 'annotation name collision(s): %s', strjoin(overlap, ', '));
            end

            if numel(rhs_shards) ~= numel(obj.shards)
                error('statgen:annotations', 'union_annotations requires matching shard structure');
            end

            out_shards = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                a = obj.shards{i};
                b = rhs_shards{i};
                if ~strcmp(a.label, b.label)
                    error('statgen:annotations', 'union_annotations requires matching shard labels');
                end
                if a.num_snp ~= b.num_snp
                    error('statgen:annotations', 'union_annotations requires matching shard row counts');
                end
                b_has_checksum = false;
                if isobject(b)
                    b_has_checksum = isprop(b, 'checksum');
                elseif isstruct(b)
                    b_has_checksum = isfield(b, 'checksum');
                end
                if b_has_checksum
                    b_checksum = b.checksum;
                    if ~strcmp(a.checksum, b_checksum)
                        error('statgen:annotations', ...
                            'union_annotations requires checksum-compatible reference alignment');
                    end
                end

                out_shards{i} = statgen.AnnotationShard( ...
                    a.label, a.checksum, [a.annomat, sparse(double(b.annomat))]);
            end

            out = statgen.AnnotationPanel(out_shards, [obj.annonames; rhs_names]);
        end
    end
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
