classdef AnnotationPanel
%STATGEN.ANNOTATIONPANEL Reference-aligned SNP annotation matrix.
%
%   annotations = statgen.load_annotations(bed_paths, reference)
%   annotations = statgen.load_annotation(path, reference, ...)
%   annotations = statgen.load_annotations_cache(path)
%   annotations = statgen.create_annotations(reference, annotation_matrix, annotation_names)
%
% An AnnotationPanel stores one or more SNP annotation columns aligned to a
% ReferencePanel. annomat is a num_snp-by-num_annot sparse numeric matrix in panel
% order.
%
% Common properties:
%   num_snp              Number of SNPs in the aligned reference.
%   num_annot            Number of annotation columns.
%   annonames            Annotation names.
%   is_binary            Logical vector marking binary columns.
%   annotation_metadata  Opaque metadata strings.
%   annomat              Sparse numeric annotation matrix.
%
% Common methods:
%   select_shards       Restrict the panel to selected shards.
%   select_annotations  Restrict the panel to selected annotation columns.
%   union_annotations   Combine non-overlapping annotation columns.
%   save_cache          Save the panel to a MATLAB .mat cache.
%
% See also statgen.load_annotations, statgen.load_annotation, statgen.create_annotations,
% statgen.create_annotation, statgen.load_annotations_cache.
    properties (SetAccess = private)
        num_snp
        num_annot
        annonames
        is_binary
        annotation_metadata
        shard_offsets
        shards
    end
    properties (Dependent)
        annomat
    end

    methods
        function obj = AnnotationPanel(shards_cell, annonames, is_binary, annotation_metadata)
            if nargin == 0, return; end
            if nargin < 3, is_binary = []; end
            if nargin < 4, annotation_metadata = []; end
            if isempty(shards_cell)
                error('statgen:annotations', 'AnnotationPanel requires at least one shard');
            end
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
            if nargin < 3 || isempty(is_binary)
                error('statgen:annotations', 'is_binary must be supplied');
            end
            obj.is_binary = coerce_is_binary_(is_binary, obj.num_annot);
            if isempty(annotation_metadata)
                obj.annotation_metadata = repmat({''}, obj.num_annot, 1);
            else
                obj.annotation_metadata = statgen.internal.annotation_coerce_metadata_vector( ...
                    annotation_metadata, obj.num_annot, 'annotation_metadata');
            end
        end

        function display(obj)
            name = inputname(1);
            if ~isempty(name)
                fprintf('%s =\n\n', name);
            end
            disp(obj);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  statgen.AnnotationPanel array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            nnz_total = count_nnz_(obj.shards);
            denom = max(1, obj.num_snp * obj.num_annot);
            fprintf('  statgen.AnnotationPanel object\n\n');
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    num_annot: %d\n', obj.num_annot);
            fprintf('    shards: %d\n', numel(obj.shards));
            fprintf('    shard_labels: %s\n', shard_labels_string_(obj.shards));
            fprintf('    annonames: %s\n', statgen.internal.display_join_strings(obj.annonames));
            fprintf('    annomat: %d-by-%d sparse numeric, nnz=%d, density=%.4g\n', ...
                obj.num_snp, obj.num_annot, nnz_total, nnz_total / denom);
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
        %SELECT_SHARDS Return annotations restricted to selected shards.
        %
        %   out = annotations.select_shards(shards)
        %
        % shards is a cell array or string array of canonical shard labels, for
        % example {'21', '22'}. The returned AnnotationPanel preserves the
        % requested shard order.
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.internal.validate_requested_shards(shards, available, 'AnnotationPanel.select_shards');
            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.AnnotationPanel(out_shards, obj.annonames, obj.is_binary, obj.annotation_metadata);
        end

        function out = select_annotations(obj, names)
        %SELECT_ANNOTATIONS Return selected annotation columns by name.
        %
        %   out = annotations.select_annotations(names)
        %
        % names is a character array, string array, or cell array naming
        % annotation columns in annotations.annonames. The output preserves the
        % requested annotation order.
        %
        % See also statgen.AnnotationPanel.select_shards,
        % statgen.AnnotationPanel.union_annotations.
            req = ensure_names_(names);
            [tf, idx] = ismember(req, obj.annonames);
            if ~all(tf)
                missing = req(~tf);
                error('statgen:annotations', 'unknown annotation name(s): %s', strjoin(missing, ', '));
            end

            out_shards = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                s = obj.shards{i};
                out_shards{i} = statgen.AnnotationShard(s.label, s.reference_checksum, s.annomat(:, idx));
            end
            out = statgen.AnnotationPanel(out_shards, req, obj.is_binary(idx), obj.annotation_metadata(idx));
        end

        function out = union_annotations(obj, other, mode)
        %UNION_ANNOTATIONS Combine non-overlapping annotation columns.
        %
        %   out = annotations.union_annotations(other)
        %
        % Returns an AnnotationPanel with columns from both inputs. The inputs
        % must use the same reference alignment and distinct annotation names.
        %
        % See also statgen.AnnotationPanel.select_annotations.
            if nargin < 3 || isempty(mode)
                mode = 'by_name';
            end
            mode = char(mode);
            if ~strcmp(mode, 'by_name')
                error('statgen:annotations', 'union_annotations supports only mode=''by_name''');
            end
            if ~isa(other, 'statgen.AnnotationPanel')
                error('statgen:annotations', 'union_annotations requires another AnnotationPanel');
            end
            rhs_names = ensure_names_(other.annonames);
            rhs_shards = other.shards;
            rhs_is_binary = coerce_is_binary_(other.is_binary, numel(rhs_names));
            rhs_metadata = statgen.internal.annotation_coerce_metadata_vector( ...
                other.annotation_metadata, numel(rhs_names), 'annotation_metadata');

            overlap = intersect(obj.annonames, rhs_names, 'stable');
            if ~isempty(overlap)
                error('statgen:annotations', 'annotation name collision(s): %s', strjoin(overlap, ', '));
            end

            if numel(rhs_shards) ~= numel(obj.shards)
                error('statgen:annotations', 'union_annotations shard count mismatch');
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
                b_checksum = b.reference_checksum;
                if ~strcmp(a.reference_checksum, b_checksum)
                    error('statgen:annotations', 'union_annotations reference_checksum mismatch');
                end

                out_shards{i} = statgen.AnnotationShard( ...
                    a.label, a.reference_checksum, [a.annomat, b.annomat]);
            end

            out = statgen.AnnotationPanel(out_shards, [obj.annonames; rhs_names], ...
                [obj.is_binary; rhs_is_binary], [obj.annotation_metadata; rhs_metadata]);
        end

        function save_cache(obj, path, varargin)
        %SAVE_CACHE Save annotations to a MATLAB .mat cache.
        %
        %   annotations.save_cache(path)
        %   annotations.save_cache(path, 'format', format)
        %
        % Saves the panel in the same format as statgen.save_annotations_cache.
        %
        % Default format is 'v7.3' in MATLAB and 'v7' in Octave.
        %
        % See also statgen.save_annotations_cache,
        % statgen.load_annotations_cache.
            statgen.save_annotations_cache(obj, path, varargin{:});
        end
    end
end

function out = coerce_is_binary_(value, expected_len)
    if ischar(value) || isstring(value)
        error('statgen:annotations', 'is_binary must be a logical vector with length %d', expected_len);
    end
    arr = value(:);
    if numel(arr) ~= expected_len
        error('statgen:annotations', 'is_binary length mismatch: expected %d, got %d', expected_len, numel(arr));
    end
    if ~islogical(arr)
        bad = ~(arr == 0 | arr == 1);
        if any(bad)
            i = find(bad, 1, 'first');
            error('statgen:annotations', 'is_binary(%d) must be logical', i);
        end
    end
    out = logical(arr);
end

function out = ensure_names_(names)
    if ischar(names) || isstring(names)
        names = cellstr(names(:));
    end
    if ~iscell(names)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end

    out = statgen.internal.ensure_cell_col(names);
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

function out = shard_labels_string_(shards)
    labels = cell(numel(shards), 1);
    for i = 1:numel(shards)
        labels{i} = shards{i}.label;
    end
    out = statgen.internal.display_join_strings(labels);
end

function n = count_nnz_(shards)
    n = 0;
    for i = 1:numel(shards)
        n = n + nnz(shards{i}.annomat);
    end
end
