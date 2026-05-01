function panel = load_reference(path, shards)
% Load a reference panel from a .bim file or sharded .bim template.
%
%   panel = statgen.load_reference(path)
%   panel = statgen.load_reference(path, shards)
%
% If path contains '@', it is a sharded template (e.g. 'chr@.bim') and
% shard discovery uses canonical-label substitution only (1-22, X).
% Otherwise path is a single file split by chr column into canonical shards.
    if nargin < 2
        shards = [];
    end

    path = char(path);
    if ~isempty(strfind(path, '@'))
        panel = load_sharded_(path, shards);
    else
        bim = parse_bim_(path);
        panel = load_split_by_chr_(bim, shards);
    end
end

% ---------------------------------------------------------------------------

function panel = load_sharded_(path, requested_shards)
    canonical = canonical_labels_();
    available_labels = {};
    available_paths = {};
    for i = 1:numel(canonical)
        label = canonical{i};
        bim_path = strrep(path, '@', label);
        if exist(bim_path, 'file') == 2
            available_labels{end+1} = label; %#ok<AGROW>
            available_paths{end+1} = bim_path; %#ok<AGROW>
        end
    end

    if isempty(available_labels)
        error('statgen:io', 'No BIM shards found matching template: %s', path);
    end

    selected = validate_requested_shards_(requested_shards, available_labels, 'load_reference');
    shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        path_idx = find(strcmp(available_labels, label), 1, 'first');
        bim = parse_bim_(available_paths{path_idx});
        shards{i} = statgen.ReferenceShard(label, bim.chr, bim.snp, bim.bp, bim.a1, bim.a2);
    end
    panel = statgen.ReferencePanel(shards);
end

function panel = load_split_by_chr_(bim, requested_shards)
    canonical = canonical_labels_();
    available = {};
    for i = 1:numel(canonical)
        c = canonical{i};
        if any(strcmp(bim.chr, c))
            available{end+1} = c; %#ok<AGROW>
        end
    end

    selected = validate_requested_shards_(requested_shards, available, 'load_reference');
    shards = cell(numel(selected), 1);
    for ci = 1:numel(selected)
        c = selected{ci};
        mask = strcmp(bim.chr, c);
        shards{ci} = statgen.ReferenceShard(c, ...
            bim.chr(mask), bim.snp(mask), ...
            bim.bp(mask), bim.a1(mask), bim.a2(mask));
    end

    panel = statgen.ReferencePanel(shards);
end

function bim = parse_bim_(path)
    path = char(path);
    [raw_cols, n_rows] = read_bim_tabular_(path);

    if n_rows == 0
        error('statgen:bim', 'BIM file is empty: %s', path);
    end

    chr_out = raw_cols{1};
    snp_out = raw_cols{2};
    cm_raw  = raw_cols{3};
    bp_raw  = raw_cols{4};
    a1_out  = raw_cols{5};
    a2_out  = raw_cols{6};

    bad_chr = cellfun('isempty', chr_out);
    if any(bad_chr)
        lineno = find(bad_chr, 1, 'first');
        error('statgen:bim', '%s:%d: chr must be non-empty', path, lineno);
    end

    chr_style = cellfun(@(c) strncmpi(c, 'chr', 3), chr_out);
    if any(chr_style)
        lineno = find(chr_style, 1, 'first');
        error('statgen:bim', '%s:%d: chr-style labels (e.g., chr1/chrX) are not allowed', path, lineno);
    end

    canonical = canonical_labels_();
    known = ismember(chr_out, canonical) | ismember(chr_out, {'Y', 'MT'});
    if ~all(known)
        lineno = find(~known, 1, 'first');
        error('statgen:bim', ...
            '%s:%d: unsupported chr label %s; expected 1-22, X (Y/MT are ignored)', ...
            path, lineno, chr_out{lineno});
    end

    bad_allele = cellfun('isempty', a1_out) | cellfun('isempty', a2_out);
    if any(bad_allele)
        lineno = find(bad_allele, 1, 'first');
        error('statgen:bim', '%s:%d: a1 and a2 must be non-empty', path, lineno);
    end

    bad_a1 = cellfun(@(a) isempty(regexp(a, '^[ACGT]+$', 'once')), a1_out);
    if any(bad_a1)
        lineno = find(bad_a1, 1, 'first');
        error('statgen:bim', ...
            '%s:%d: a1 must be uppercase DNA bases (A/C/G/T): %s', path, lineno, a1_out{lineno});
    end
    bad_a2 = cellfun(@(a) isempty(regexp(a, '^[ACGT]+$', 'once')), a2_out);
    if any(bad_a2)
        lineno = find(bad_a2, 1, 'first');
        error('statgen:bim', ...
            '%s:%d: a2 must be uppercase DNA bases (A/C/G/T): %s', path, lineno, a2_out{lineno});
    end

    cm_out = str2double(cm_raw);
    bad_cm = isnan(cm_out);
    if any(bad_cm)
        lineno = find(bad_cm, 1, 'first');
        error('statgen:bim', '%s:%d: cm is not a number: %s', path, lineno, cm_raw{lineno});
    end

    bp_out = str2double(bp_raw);
    bad_bp = isnan(bp_out) | (bp_out ~= floor(bp_out));
    if any(bad_bp)
        lineno = find(bad_bp, 1, 'first');
        error('statgen:bim', '%s:%d: bp is not an integer: %s', path, lineno, bp_raw{lineno});
    end

    keep = ismember(chr_out, canonical);
    line_idx = find(keep);
    validate_reference_sort_order_( ...
        chr_out(keep), bp_out(keep), a1_out(keep), a2_out(keep), path, line_idx);

    bim.chr = chr_out(keep);
    bim.snp = snp_out(keep);
    bim.bp  = bp_out(keep);
    bim.a1  = a1_out(keep);
    bim.a2  = a2_out(keep);
end

function validate_reference_sort_order_(chr_col, bp_col, a1_col, a2_col, path, line_idx)
    if isempty(chr_col)
        return
    end

    canonical = canonical_labels_();
    [is_ok_chr, chr_rank] = ismember(chr_col, canonical);
    bad_chr = ~is_ok_chr;
    if any(bad_chr)
        bad_i = find(bad_chr, 1, 'first');
        error('statgen:bim', '%s:%d: chr must use canonical labels 1-22 or X', ...
            path, line_idx(bad_i));
    end

    n = numel(chr_col);
    rank_str = arrayfun(@(x) sprintf('%03d', x), chr_rank(:), 'UniformOutput', false);
    bp_str = arrayfun(@(x) sprintf('%012d', round(x)), bp_col(:), 'UniformOutput', false);
    keys = strcat(rank_str, ':', bp_str, ':', a1_col(:), ':', a2_col(:));

    if n > 1
        sorted_keys = sort(keys);
        dup_sorted = strcmp(sorted_keys(2:end), sorted_keys(1:end-1));
        if any(dup_sorted)
            dup_pos = find(dup_sorted, 1, 'first') + 1;
            dup_key = sorted_keys{dup_pos};
            occ = find(strcmp(keys, dup_key), 2, 'first');
            if numel(occ) >= 2
                bad_line = line_idx(occ(2));
            else
                bad_line = line_idx(occ(1));
            end
            error('statgen:bim', ...
                '%s:%d: duplicate (chr, bp, a1, a2) tuple is not allowed', ...
                path, bad_line);
        end
    end

    [~, sort_idx] = sort(keys);
    expected = (1:n)';
    bad = find(sort_idx(:) ~= expected, 1, 'first');
    if ~isempty(bad)
        bad_line = line_idx(bad);
        error('statgen:bim', ...
            '%s:%d: rows must be sorted by (chr_rank, bp, a1, a2) in canonical contig order', ...
            path, bad_line);
    end
end

function labels = canonical_labels_()
    labels = [arrayfun(@num2str, 1:22, 'UniformOutput', false), {'X'}];
end

function selected = validate_requested_shards_(requested, available, where)
    if nargin < 1 || isempty(requested)
        selected = available;
        return
    end

    if ischar(requested)
        requested = {requested};
    elseif isnumeric(requested)
        error('statgen:shards', '%s: shards must be a non-empty list of unique canonical contig labels', where);
    end

    requested = requested(:)';
    if isempty(requested)
        error('statgen:shards', '%s: shards must be a non-empty list of unique canonical contig labels', where);
    end

    canonical = canonical_labels_();
    canonical_idx = zeros(1, numel(requested));
    seen = {};
    for i = 1:numel(requested)
        label = char(requested{i});
        idx = find(strcmp(canonical, label), 1, 'first');
        if isempty(idx)
            error('statgen:shards', ...
                '%s: unsupported shard label %s; expected canonical labels 1-22 or X', ...
                where, label);
        end
        if any(strcmp(seen, label))
            error('statgen:shards', '%s: duplicate shard label %s in shards', where, label);
        end
        if ~any(strcmp(available, label))
            error('statgen:shards', '%s: requested shard %s is not present', where, label);
        end
        seen{end+1} = label; %#ok<AGROW>
        canonical_idx(i) = idx;
    end

    if any(diff(canonical_idx) <= 0)
        error('statgen:shards', '%s: shards must be in canonical subsequence order', where);
    end

    selected = requested;
end

function [cols, n_rows] = read_bim_tabular_(path)
    if exist('readtable', 'file') == 2
        tbl = readtable(path, ...
            'FileType', 'text', ...
            'Delimiter', '\t', ...
            'ReadVariableNames', false, ...
            'TextType', 'char');

        if width(tbl) ~= 6
            error('statgen:bim', '%s: expected 6 tab-separated columns, got %d', ...
                path, width(tbl));
        end

        n_rows = height(tbl);
        cols = cell(1, 6);
        for c = 1:6
            cols{c} = ensure_cellstr_col_(tbl{:, c});
        end
        return
    end

    % Octave compatibility: readtable is unavailable in some builds.
    % Use textscan as a native tabular-reader equivalent.
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open BIM file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    raw_cols = textscan(fid, '%s%s%s%s%s%s', ...
        'Delimiter', '\t', ...
        'Whitespace', '', ...
        'MultipleDelimsAsOne', false, ...
        'ReturnOnError', false);
    clear cleaner;

    cols = cell(1, 6);
    for c = 1:6
        cols{c} = ensure_cellstr_col_(raw_cols{c});
    end
    n_rows = numel(cols{1});
end

function out = ensure_cellstr_col_(x)
    if iscell(x)
        out = x(:);
    elseif ischar(x)
        out = strtrim(cellstr(x));
    else
        out = strtrim(cellstr(x));
    end
end
