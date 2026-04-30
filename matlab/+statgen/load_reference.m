function panel = load_reference(path, no_shard)
% Load a reference panel from a .bim file or sharded .bim template.
%
%   panel = statgen.load_reference(path)
%   panel = statgen.load_reference(path, no_shard)
%
% If path contains '@', it is a sharded template (e.g. 'chr@.bim').
% Otherwise path is a single file split by chr column (no_shard=false)
% or loaded as one 'all' shard (no_shard=true).
    if nargin < 2, no_shard = false; end
    path = char(path);

    if ~isempty(strfind(path, '@'))
        panel = load_sharded_(path);
    else
        bim = parse_bim_(path);
        if no_shard
            shard = statgen.ReferenceShard('all', bim.chr, bim.snp, bim.bp, bim.a1, bim.a2);
            panel = statgen.ReferencePanel({shard});
        else
            panel = load_split_by_chr_(bim);
        end
    end
end

% ---------------------------------------------------------------------------

function panel = load_sharded_(path)
    at_pos = strfind(path, '@');
    at_pos = at_pos(1);
    prefix = path(1:at_pos-1);
    suffix = path(at_pos+1:end);

    % Split prefix into parent directory and name prefix
    sep_pos = max([strfind(prefix, '/'), strfind(prefix, '\')]);
    if isempty(sep_pos)
        parent_dir = '.';
        name_prefix = prefix;
    else
        parent_dir = prefix(1:sep_pos);
        name_prefix = prefix(sep_pos+1:end);
    end
    name_suffix = suffix;

    listing = dir(fullfile(parent_dir, [name_prefix, '*', name_suffix]));
    listing = listing(~[listing.isdir]);

    if isempty(listing)
        error('statgen:io', 'No BIM shards found matching template: %s', path);
    end

    % Extract shard labels from filenames
    labels = cell(numel(listing), 1);
    for i = 1:numel(listing)
        fname = listing(i).name;
        if isempty(name_suffix)
            inner = fname(length(name_prefix)+1:end);
        else
            inner = fname(length(name_prefix)+1:end-length(name_suffix));
        end
        labels{i} = inner;
    end

    % Sort by canonical chromosome order
    [labels, sort_idx] = sort_chr_labels_(labels);
    listing = listing(sort_idx);

    % Load each shard using the template path
    shards = cell(numel(listing), 1);
    for i = 1:numel(listing)
        bim_path = [prefix, labels{i}, suffix];
        bim = parse_bim_(bim_path);
        shards{i} = statgen.ReferenceShard(labels{i}, bim.chr, bim.snp, bim.bp, bim.a1, bim.a2);
    end

    panel = statgen.ReferencePanel(shards);
end

function panel = load_split_by_chr_(bim)
    % Collect unique chromosomes in first-appearance order
    chr_labels = {};
    for i = 1:numel(bim.chr)
        c = bim.chr{i};
        if ~any(strcmp(chr_labels, c))
            chr_labels{end+1} = c;
        end
    end

    shards = cell(numel(chr_labels), 1);
    for ci = 1:numel(chr_labels)
        c = chr_labels{ci};
        mask = strcmp(bim.chr, c);
        shards{ci} = statgen.ReferenceShard(c, ...
            bim.chr(mask), bim.snp(mask), ...
            bim.bp(mask), bim.a1(mask), bim.a2(mask));
    end

    panel = statgen.ReferencePanel(shards);
end

function [sorted_labels, idx] = sort_chr_labels_(labels)
    n = numel(labels);
    num_vals = nan(n, 1);
    for i = 1:n
        v = str2double(labels{i});
        if ~isnan(v)
            num_vals(i) = v;
        end
    end

    num_mask = ~isnan(num_vals);
    num_idx  = find(num_mask);
    str_idx  = find(~num_mask);

    [~, num_sort] = sort(num_vals(num_idx));
    num_idx_sorted = num_idx(num_sort);

    [~, str_sort] = sort(labels(str_idx));
    str_idx_sorted = str_idx(str_sort);

    idx = [num_idx_sorted; str_idx_sorted];
    sorted_labels = labels(idx);
end

function bim = parse_bim_(path)
    % Parse a PLINK .bim file using native tabular readers.
    % Default path uses readtable when available; Octave fallback uses
    % textscan (native table reader equivalent) to preserve compatibility.
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

    bad_allele = cellfun('isempty', a1_out) | cellfun('isempty', a2_out);
    if any(bad_allele)
        lineno = find(bad_allele, 1, 'first');
        error('statgen:bim', '%s:%d: a1 and a2 must be non-empty', path, lineno);
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

    validate_reference_sort_order_(chr_out, bp_out, path);

    bim.chr = chr_out(:);
    bim.snp = snp_out(:);
    bim.bp  = bp_out(:);
    bim.a1  = a1_out(:);
    bim.a2  = a2_out(:);
end

function validate_reference_sort_order_(chr_col, bp_col, path)
    canonical = [arrayfun(@num2str, 1:22, 'UniformOutput', false), {'X'}];

    [is_ok_chr, chr_rank] = ismember(chr_col, canonical);
    bad_chr = ~is_ok_chr;
    if any(bad_chr)
        lineno = find(bad_chr, 1, 'first');
        error('statgen:bim', '%s:%d: chr must use canonical labels 1-22 or X', path, lineno);
    end

    bad_order = (chr_rank(2:end) < chr_rank(1:end-1)) | ...
        ((chr_rank(2:end) == chr_rank(1:end-1)) & (bp_col(2:end) < bp_col(1:end-1)));
    if any(bad_order)
        lineno = find(bad_order, 1, 'first') + 1;
        error('statgen:bim', ...
            '%s:%d: rows must be sorted by canonical chromosome (1-22, X) and ascending bp', ...
            path, lineno);
    end
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
        return;
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
