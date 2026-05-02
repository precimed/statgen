function panel = load_annotations(bed_paths, reference)
% Load BED annotations and paint to reference-aligned SNP masks.
    if nargin < 2
        error('statgen:arg', 'load_annotations requires bed_paths and reference');
    end

    paths = coerce_bed_paths_(bed_paths);
    annonames = cell(numel(paths), 1);
    for i = 1:numel(paths)
        if exist(paths{i}, 'file') ~= 2
            error('statgen:io', 'BED file not found: %s', paths{i});
        end
        [~, stem, ~] = fileparts(paths{i});
        if isempty(stem)
            error('statgen:annotations', 'invalid annotation filename: %s', paths{i});
        end
        annonames{i} = stem;
    end
    if numel(unique(annonames)) ~= numel(annonames)
        error('statgen:annotations', 'duplicate annotation names derived from BED basenames');
    end

    n = reference.num_snp;
    k = numel(paths);
    if n == 0
        panel = statgen.create_annotations(reference, sparse([], [], [], 0, k), annonames);
        return
    end

    A = sparse(n, k);
    for j = 1:k
        intervals_map = parse_bed_intervals_(paths{j});
        for s = 1:numel(reference.shards)
            ref_shard = reference.shards{s};
            off = reference.shard_offsets(s);
            ix = (off.start0 + 1):off.stop0;
            mask = zeros(numel(ix), 1);
            idx = find(strcmp(intervals_map.labels, ref_shard.label), 1, 'first');
            if ~isempty(idx)
                mask = paint_mask_(ref_shard.bp, intervals_map.intervals{idx});
            end
            nz = find(mask ~= 0);
            if ~isempty(nz)
                A(ix(nz), j) = 1;
            end
        end
    end

    panel = statgen.create_annotations(reference, A, annonames);
end

function paths = coerce_bed_paths_(bed_paths)
    if ischar(bed_paths)
        bed_paths = {bed_paths};
    elseif isstring(bed_paths)
        bed_paths = cellstr(bed_paths(:));
    end
    if ~iscell(bed_paths)
        error('statgen:annotations', 'bed_paths must be a non-empty list of BED files');
    end

    paths = statgen.ensure_cell_col(bed_paths);
    if isempty(paths)
        error('statgen:annotations', 'bed_paths must be a non-empty list of BED files');
    end
    for i = 1:numel(paths)
        paths{i} = char(paths{i});
    end
end

function intervals_map = parse_bed_intervals_(path)
    [cols, n_rows] = read_bed_tabular_(path);
    if n_rows == 0
        error('statgen:annotations', '%s: BED file is empty', path);
    end
    chr_col = statgen.ensure_cell_col(cols{1});
    start_raw = statgen.ensure_cell_col(cols{2});
    end_raw = statgen.ensure_cell_col(cols{3});

    bad_chr = cellfun('isempty', chr_col);
    if any(bad_chr)
        i = find(bad_chr, 1, 'first');
        error('statgen:annotations', '%s: row %d: chromosome label must be non-empty', path, i);
    end

    start_num = str2double(start_raw);
    bad_start = isnan(start_num) | (start_num ~= floor(start_num)) | (start_num < 0);
    if any(bad_start)
        i = find(bad_start, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED start must be a non-negative integer', path, i);
    end

    end_num = str2double(end_raw);
    bad_end = isnan(end_num) | (end_num ~= floor(end_num)) | (end_num < 0);
    if any(bad_end)
        i = find(bad_end, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED end must be a non-negative integer', path, i);
    end

    bad_len = end_num < start_num;
    if any(bad_len)
        i = find(bad_len, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED interval end must be >= start', path, i);
    end

    intervals_map = struct('labels', {{}}, 'intervals', {{}});

    chr_unique = {};
    for i = 1:numel(chr_col)
        c = chr_col{i};
        if ~any(strcmp(chr_unique, c))
            chr_unique{end + 1, 1} = c; %#ok<AGROW>
        end
    end

    for i = 1:numel(chr_unique)
        mask = strcmp(chr_col, chr_unique{i});
        merged = merge_intervals_(start_num(mask), end_num(mask));
        intervals_map.labels{end + 1, 1} = chr_unique{i}; %#ok<AGROW>
        intervals_map.intervals{end + 1, 1} = merged; %#ok<AGROW>
    end
end

function merged = merge_intervals_(starts, ends)
    if isempty(starts)
        merged = zeros(0, 2);
        return
    end

    pairs = sortrows([starts(:), ends(:)], [1, 2]);
    cur_s = pairs(1, 1);
    cur_e = pairs(1, 2);
    out_s = [];
    out_e = [];

    for i = 2:size(pairs, 1)
        s = pairs(i, 1);
        e = pairs(i, 2);
        if s <= cur_e
            if e > cur_e
                cur_e = e;
            end
        else
            out_s(end + 1, 1) = cur_s; %#ok<AGROW>
            out_e(end + 1, 1) = cur_e; %#ok<AGROW>
            cur_s = s;
            cur_e = e;
        end
    end

    out_s(end + 1, 1) = cur_s;
    out_e(end + 1, 1) = cur_e;
    merged = [out_s, out_e];
end

function mask = paint_mask_(bp, intervals)
    n = numel(bp);
    mask = zeros(n, 1);
    if isempty(intervals)
        return
    end

    starts = double(intervals(:, 1));
    ends = double(intervals(:, 2));
    pos0 = double(bp(:)) - 1;

    % Vectorized candidate-interval lookup:
    % idx(k) gives the rightmost start <= pos0(k), or 0 if none.
    [~, idx] = histc(pos0, [starts; inf]);
    valid = idx > 0;
    if any(valid)
        hit = valid & (pos0 < ends(max(idx, 1)));
        mask(hit) = 1;
    end
end

function [cols, n_rows] = read_bed_tabular_(path)
    % Parse line-wise to preserve the BED header/comment contract:
    % skip blank, '#', 'track ', and 'browser ' lines before tab parsing.
    % Accept any whitespace as a field separator (tabs or spaces).
    % PERF: this exception to native table-readers is required because
    % readtable/textscan do not provide this exact multi-prefix skipping
    % behavior consistently across MATLAB and Octave.
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open BED file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    raw_lines = textscan(fid, '%s', ...
        'Delimiter', '\n', ...
        'Whitespace', '', ...
        'ReturnOnError', false);
    clear cleaner;
    lines = raw_lines{1};
    kept = {};
    for i = 1:numel(lines)
        line = char(lines{i});
        if should_skip_bed_line_(line)
            continue
        end
        kept{end + 1, 1} = line; %#ok<AGROW>
    end

    n_rows = numel(kept);
    cols = {cell(n_rows, 1), cell(n_rows, 1), cell(n_rows, 1)};
    for i = 1:n_rows
        parts = regexp(strtrim(kept{i}), '\s+', 'split');
        if numel(parts) < 3
            error('statgen:annotations', '%s: BED must have at least 3 whitespace-separated columns', path);
        end
        cols{1}{i} = parts{1};
        cols{2}{i} = parts{2};
        cols{3}{i} = parts{3};
    end
end

function tf = should_skip_bed_line_(line)
    tf = isempty(line) ...
        || strncmp(line, '#', 1) ...
        || strncmp(line, 'track ', 6) ...
        || strncmp(line, 'browser ', 8);
end
