classdef ReferencePanel
%STATGEN.REFERENCEPANEL Reference SNP coordinate system for aligned objects.
%
%   reference = statgen.load_reference(path)
%   reference = statgen.load_reference(path, shards)
%   reference = statgen.load_reference_cache(path)
%
% A ReferencePanel defines the SNP order, chromosome labels, base-pair
% positions, and alleles used by other statgen objects. SNP-axis properties are
% returned as num_snp-by-1 vectors in panel order.
%
% Cache-loaded ReferencePanel objects may keep SNP identifiers and alleles in
% encoded shard-local text payloads and decode them lazily when accessed.
%
% Common properties:
%   num_snp    Number of SNPs in the panel.
%   chr        Chromosome labels.
%   snp        SNP identifiers.
%   bp         Base-pair positions.
%   a1, a2     Alleles from the BIM input.
%   is_single_nucleotide_variant
%              True for variants where both alleles are single nucleotides.
%   is_strand_ambiguous
%              True for A/T, T/A, C/G, and G/C allele pairs.
%
% Common methods:
%   select_shards           Restrict the reference to selected shards.
%   is_object_compatible    Check whether another object uses this reference.
%   save_cache              Save the reference to a MATLAB .mat cache.
%
% See also statgen.load_reference, statgen.load_reference_cache,
% statgen.save_reference_cache.
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
        a1_hash64      % num_snp×1 uint64 vector
        a2_hash64      % num_snp×1 uint64 vector
        is_single_nucleotide_variant  % num_snp×1 logical vector
        is_strand_ambiguous     % num_snp×1 logical vector
    end

    methods
        function obj = ReferencePanel(shards_cell)
            if nargin == 0, return; end
            if isempty(shards_cell)
                error('statgen:reference', 'ReferencePanel requires at least one shard');
            end
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

        function display(obj)
            name = inputname(1);
            if ~isempty(name)
                fprintf('%s =\n\n', name);
            end
            disp(obj);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  statgen.ReferencePanel array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            fprintf('  statgen.ReferencePanel object\n\n');
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    shards: %d\n', numel(obj.shards));
            if ~isempty(obj.shards)
                labels = cell(numel(obj.shards), 1);
                for i = 1:numel(obj.shards)
                    labels{i} = obj.shards{i}.label;
                end
                fprintf('    shard_labels: %s\n', statgen.internal.display_join_strings(labels));
            end
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

        function out = get.a1_hash64(obj)
            if isempty(obj.shards), out = zeros(0, 1, 'uint64'); return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.a1_hash64;
            end
            out = vertcat(vals{:});
        end

        function out = get.a2_hash64(obj)
            if isempty(obj.shards), out = zeros(0, 1, 'uint64'); return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.a2_hash64;
            end
            out = vertcat(vals{:});
        end

        function out = get.is_single_nucleotide_variant(obj)
            if isempty(obj.shards), out = false(0, 1); return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.is_single_nucleotide_variant;
            end
            out = vertcat(vals{:});
        end

        function out = get.is_strand_ambiguous(obj)
            if isempty(obj.shards), out = false(0, 1); return; end
            vals = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                vals{i} = obj.shards{i}.is_strand_ambiguous;
            end
            out = vertcat(vals{:});
        end

        function out = select_shards(obj, shards)
        %SELECT_SHARDS Return a reference restricted to selected shards.
        %
        %   out = reference.select_shards(shards)
        %
        % shards is a cell array or string array of canonical shard labels, for
        % example {'21', '22'}. The returned ReferencePanel preserves the
        % requested shard order.
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.internal.validate_requested_shards( ...
                shards, available, 'ReferencePanel.select_shards');

            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.ReferencePanel(out_shards);
        end

        function ok = validate_checksums(obj)
        %VALIDATE_CHECKSUMS Validate stored reference checksums.
        %
        %   ok = reference.validate_checksums()
        %
        % Returns true when the loaded reference fields match their stored
        % per-shard checksums. Cache-loaded references decode allele fields on
        % demand for this explicit validation path.
            for i = 1:numel(obj.shards)
                s = obj.shards{i};
                validate_reference_checksum_(s);
            end
            ok = true;
        end

        function ok = is_object_compatible(obj, other)
        %IS_OBJECT_COMPATIBLE Check whether another object uses this reference.
        %
        %   ok = reference.is_object_compatible(other)
        %
        % Returns true when other has the same shard labels, SNP counts, and
        % reference checksums as this ReferencePanel. Mismatches are reported as
        % compatibility warnings and return false, not an exception.
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
                    os_chk = os.reference_checksum;
                catch
                    try
                        os_chk = os.checksum;
                    catch
                        os_chk = [];
                    end
                end
                if ~isempty(os_chk) && ~strcmp(os_chk, rs.checksum)
                    statgen.ReferencePanel.compat_warn_( ...
                        'statgen: is_object_compatible: shard %s: reference_checksum mismatch', ...
                        rs.label);
                    ok = false;
                end
            end
        end

        function save_cache(obj, path, varargin)
        %SAVE_CACHE Save the reference to a MATLAB .mat cache.
        %
        %   reference.save_cache(path)
        %   reference.save_cache(path, 'format', format)
        %
        % Saves the reference in the same format as statgen.save_reference_cache.
        %
        % Default format is 'v7.3' in MATLAB and 'v7' in Octave.
        %
        % See also statgen.save_reference_cache, statgen.load_reference_cache.
            statgen.save_reference_cache(obj, path, varargin{:});
        end
    end

    methods (Static, Access = private)
        function compat_warn_(fmt, varargin)
            msg = sprintf(fmt, varargin{:});
            warning('statgen:compat', '%s', msg);
        end
    end
end

function validate_reference_checksum_(shard)
    try
        bp_str = cellstr(num2str(round(shard.bp), '%d'));
        parts = strcat(shard.chr, {':'}, bp_str, {':'}, shard.a1, {':'}, shard.a2, {sprintf('\n')});
        text_payload = [parts{:}];
        computed = statgen.internal.md5_hex(text_payload);
    catch ME
        error('statgen:cache', ...
            'Reference checksum validation failed while materializing reference fields (%s)', ...
            ME.message);
    end
    if ~strcmp(computed, shard.checksum)
        error('statgen:cache', 'Reference checksum mismatch for shard %s', shard.label);
    end
end
