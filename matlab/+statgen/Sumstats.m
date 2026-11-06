classdef Sumstats
%STATGEN.SUMSTATS Reference-aligned summary statistics for one trait.
%
%   sumstats = statgen.load_sumstats(path, reference)
%   sumstats = statgen.load_sumstats_cache(path)
%   sumstats = statgen.create_sumstats(reference, p, ...)
%
% A Sumstats object stores GWAS summary statistics aligned to a ReferencePanel.
% SNP-axis properties are returned as num_snp-by-1 vectors in panel order.
% Variants without valid p after alignment have NaN in logpvec and present
% optional vectors, with is_present false. A source p-value of 0 is represented
% as Inf in logpvec.
% Optional statistic vectors are [] when the source field was absent.
% zvec and beta_vec are oriented to reference a1 (ALT); eaf_vec is its frequency.
%
% Common properties:
%   num_snp     Number of SNPs in the aligned reference.
%   logpvec     -log10(p) values.
%   zvec        Optional z-score vector.
%   nvec        Optional sample-size vector.
%   beta_vec, se_vec, eaf_vec, info_vec
%               Further optional statistic vectors.
%   is_present  Logical vector marking variants with non-missing p.
%
% Common methods:
%   select_shards Restrict the object to selected shards.
%   save_cache    Save the object to a MATLAB .mat cache.
%
% See also statgen.load_sumstats, statgen.load_sumstats_cache,
% statgen.create_sumstats, statgen.save_sumstats_cache.
    properties (SetAccess = private)
        num_snp
        shard_offsets
        shards
    end
    properties (Dependent)
        zvec
        nvec
        logpvec
        is_present
        beta_vec
        se_vec
        eaf_vec
        info_vec
    end

    methods
        function obj = Sumstats(shards_cell)
            if nargin == 0, return; end
            if isempty(shards_cell)
                error('statgen:sumstats', 'Sumstats requires at least one shard');
            end
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

        function display(obj)
            name = inputname(1);
            if ~isempty(name)
                fprintf('%s =\n\n', name);
            end
            disp(obj);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  statgen.Sumstats array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            fprintf('  statgen.Sumstats object\n\n');
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    shards: %d\n', numel(obj.shards));
            fprintf('    shard_labels: %s\n', shard_labels_string_(obj.shards));
            fprintf('    present_snps: %d\n', count_present_(obj.shards));
            fprintf('    logpvec: %d-by-1 double\n', obj.num_snp);
            [present, absent] = optional_field_groups_(obj.shards, ...
                {'zvec', 'nvec', 'beta_vec', 'se_vec', 'eaf_vec', 'info_vec'});
            fprintf('    optional_present: %s\n', statgen.internal.display_join_strings(present));
            fprintf('    optional_absent: %s\n', statgen.internal.display_join_strings(absent));
        end

        function out = get.zvec(obj)
            out = concat_optional_(obj.shards, 'zvec');
        end

        function out = get.nvec(obj)
            out = concat_optional_(obj.shards, 'nvec');
        end

        function out = get.logpvec(obj)
            out = concat_required_(obj.shards, 'logpvec');
        end

        function out = get.is_present(obj)
            out = ~isnan(obj.logpvec);
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
        %SELECT_SHARDS Return summary statistics restricted to selected shards.
        %
        %   out = sumstats.select_shards(shards)
        %
        % shards is a cell array or string array of canonical shard labels, for
        % example {'21', '22'}. The returned Sumstats object preserves the
        % requested shard order.
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.internal.validate_requested_shards(shards, available, 'Sumstats.select_shards');
            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.Sumstats(out_shards);
        end

        function save_cache(obj, path, varargin)
        %SAVE_CACHE Save summary statistics to a MATLAB .mat cache.
        %
        %   sumstats.save_cache(path)
        %   sumstats.save_cache(path, 'format', format)
        %
        % Saves the object in the same format as statgen.save_sumstats_cache.
        %
        % Default format is 'v7.3' in MATLAB and 'v7' in Octave.
        %
        % See also statgen.save_sumstats_cache, statgen.load_sumstats_cache.
            statgen.save_sumstats_cache(obj, path, varargin{:});
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

function out = shard_labels_string_(shards)
    labels = cell(numel(shards), 1);
    for i = 1:numel(shards)
        labels{i} = shards{i}.label;
    end
    out = statgen.internal.display_join_strings(labels);
end

function n = count_present_(shards)
    n = 0;
    for i = 1:numel(shards)
        n = n + sum(shards{i}.is_present);
    end
end

function [present, absent] = optional_field_groups_(shards, fields)
    present = {};
    absent = {};
    for j = 1:numel(fields)
        field = fields{j};
        has_field = true;
        for i = 1:numel(shards)
            if isempty(shards{i}.(field))
                has_field = false;
                break;
            end
        end
        if has_field
            present{end + 1, 1} = field;
        else
            absent{end + 1, 1} = field;
        end
    end
end
