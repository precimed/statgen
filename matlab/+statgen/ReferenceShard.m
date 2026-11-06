classdef ReferenceShard
% Immutable in-memory representation of one .bim shard.
    properties (SetAccess = private)
        label      % chromosome label string
        num_snp    % scalar integer
        bp         % n×1 double vector
        a1_hash64  % n×1 uint64 vector
        a2_hash64  % n×1 uint64 vector
        checksum   % lowercase MD5 hex string
    end
    properties (Dependent)
        chr        % n×1 cell array of strings
        snp        % n×1 cell array of strings
        a1         % n×1 cell array of strings
        a2         % n×1 cell array of strings
        is_single_nucleotide_variant  % n×1 logical vector
        is_strand_ambiguous     % n×1 logical vector
    end
    properties (Access = private)
        chr_data
        snp_data
        a1_data
        a2_data
        snp_text
        a1_text
        a2_text
        is_text_cache = false
    end

    methods
        function obj = ReferenceShard(label, chr_vec, snp_vec, bp_vec, a1_vec, a2_vec, checksum)
            if nargin == 0, return; end
            if nargin < 7, checksum = ''; end
            obj.label   = char(label);
            obj.num_snp = numel(chr_vec);
            obj.chr_data = statgen.internal.ensure_cell_col(chr_vec);
            obj.snp_data = statgen.internal.ensure_cell_col(snp_vec);
            obj.bp      = double(bp_vec(:));
            obj.a1_data = statgen.internal.ensure_cell_col(a1_vec);
            obj.a2_data = statgen.internal.ensure_cell_col(a2_vec);
            if numel(obj.chr_data) ~= obj.num_snp || numel(obj.snp_data) ~= obj.num_snp || numel(obj.bp) ~= obj.num_snp || numel(obj.a1_data) ~= obj.num_snp || numel(obj.a2_data) ~= obj.num_snp
                error('statgen:reference', 'ReferenceShard vector lengths must match');
            end
            bad = ~strcmp(obj.chr_data, obj.label);
            if any(bad)
                error('statgen:reference', 'Reference shard %s contains multiple chr labels', obj.label);
            end
            obj.a1_hash64 = statgen.internal.allele_hash64(obj.a1_data);
            obj.a2_hash64 = statgen.internal.allele_hash64(obj.a2_data);
            if ~isempty(checksum)
                obj.checksum = char(checksum);
            else
                obj.checksum = reference_checksum_(obj.chr_data, obj.bp, obj.a1_data, obj.a2_data);
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
                fprintf('  statgen.ReferenceShard array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            fprintf('  statgen.ReferenceShard object\n\n');
            fprintf('    label: %s\n', obj.label);
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    checksum: %s\n', obj.checksum);
            if obj.is_text_cache
                fprintf('    chr: synthesized, snp/a1/a2: lazy text cache payloads\n');
            else
                fprintf('    chr, snp, a1, a2: %d-by-1 cell\n', obj.num_snp);
            end
        end

        function out = get.chr(obj)
            if obj.is_text_cache || isempty(obj.chr_data)
                out = repmat({obj.label}, obj.num_snp, 1);
            else
                out = obj.chr_data;
            end
        end

        function out = get.snp(obj)
            if obj.is_text_cache
                out = decode_text_vector_(obj.snp_text, obj.num_snp, 'snp', obj.label);
            else
                out = obj.snp_data;
            end
        end

        function out = get.a1(obj)
            if obj.is_text_cache
                out = decode_text_vector_(obj.a1_text, obj.num_snp, 'a1', obj.label);
            else
                out = obj.a1_data;
            end
        end

        function out = get.a2(obj)
            if obj.is_text_cache
                out = decode_text_vector_(obj.a2_text, obj.num_snp, 'a2', obj.label);
            else
                out = obj.a2_data;
            end
        end

        function out = get.is_single_nucleotide_variant(obj)
            [out, ~] = statgen.internal.variant_hash_masks(obj.a1_hash64, obj.a2_hash64);
        end

        function out = get.is_strand_ambiguous(obj)
            [~, out] = statgen.internal.variant_hash_masks(obj.a1_hash64, obj.a2_hash64);
        end

        function [snp_text, a1_text, a2_text] = cache_text_payloads(obj)
        %CACHE_TEXT_PAYLOADS Return per-shard newline-delimited cache strings.
            if obj.is_text_cache
                snp_text = obj.snp_text;
                a1_text = obj.a1_text;
                a2_text = obj.a2_text;
                return;
            end
            snp_text = encode_text_vector_(obj.snp, obj.num_snp, 'snp', obj.label);
            a1_text = encode_text_vector_(obj.a1, obj.num_snp, 'a1', obj.label);
            a2_text = encode_text_vector_(obj.a2, obj.num_snp, 'a2', obj.label);
        end
    end

    methods (Static)
        function obj = from_cache_text(label, num_snp, bp_vec, a1_hash64, a2_hash64, checksum, snp_text, a1_text, a2_text)
            obj = statgen.ReferenceShard();
            obj.label = char(label);
            obj.num_snp = double(num_snp);
            obj.bp = double(bp_vec(:));
            obj.a1_hash64 = uint64(a1_hash64(:));
            obj.a2_hash64 = uint64(a2_hash64(:));
            obj.checksum = char(checksum);
            obj.snp_text = char(snp_text);
            obj.a1_text = char(a1_text);
            obj.a2_text = char(a2_text);
            obj.is_text_cache = true;
            if numel(obj.bp) ~= obj.num_snp || numel(obj.a1_hash64) ~= obj.num_snp || numel(obj.a2_hash64) ~= obj.num_snp
                error('statgen:cache', 'Invalid reference cache: payload vector lengths mismatch');
            end
        end
    end
end

function out = encode_text_vector_(values, expected_n, field_name, label)
    values = statgen.internal.ensure_cell_col(values);
    if numel(values) ~= expected_n
        error('statgen:cache', 'Reference shard %s cache field %s length mismatch', label, field_name);
    end
    if expected_n == 0
        out = '';
        return;
    end
    newline = sprintf('\n');
    has_newline = cellfun(@(x) ~isempty(strfind(x, newline)), values);
    if any(has_newline)
        idx = find(has_newline, 1, 'first');
        error('statgen:cache', 'Reference shard %s field %s row %d contains a newline', label, field_name, idx);
    end
    out = strjoin(values(:)', newline);
end

function out = decode_text_vector_(payload, expected_n, field_name, label)
    if expected_n == 0
        out = cell(0, 1);
        return;
    end
    out = strsplit(char(payload), sprintf('\n'))';
    if numel(out) ~= expected_n
        error('statgen:cache', 'Invalid reference cache: shard %s field %s decoded length mismatch', label, field_name);
    end
    empty = cellfun('isempty', out);
    if any(empty)
        idx = find(empty, 1, 'first');
        error('statgen:cache', 'Invalid reference cache: shard %s field %s row %d is empty', label, field_name, idx);
    end
end

function out = reference_checksum_(chr, bp, a1, a2)
    % MD5 over 'chr:bp:a1:a2\n' lines in row order
    bp_str = strtrim(cellstr(num2str(round(bp(:)), '%d')));
    parts = strcat(chr(:), {':'}, bp_str, {':'}, a1(:), {':'}, a2(:), {sprintf('\n')});
    text_payload = [parts{:}];
    out = statgen.internal.md5_hex(text_payload);
end
