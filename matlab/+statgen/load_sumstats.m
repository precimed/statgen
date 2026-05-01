function sumstats = load_sumstats(path, reference)
% Load sumstats TSV(.gz) and align to reference by exact chr:bp:a1:a2 key.
    if nargin < 2
        error('statgen:arg', 'load_sumstats requires path and reference');
    end
    path = char(path);
    if isempty(path)
        error('statgen:arg', 'path must be non-empty');
    end

    [tbl, cleanup_fn] = parse_sumstats_table_(path);
    cleaner = onCleanup(cleanup_fn);
    sumstats = build_sumstats_(tbl, reference, path);
    clear cleaner;
end

function sumstats = build_sumstats_(tbl, reference, path)
    if is_table_like_(tbl)
        var_names = lower(tbl.Properties.VariableNames);
    else
        var_names = lower(tbl._var_names(:)');
    end
    required = {'chr', 'bp', 'a1', 'a2', 'z', 'n'};
    for i = 1:numel(required)
        if ~any(strcmp(var_names, required{i}))
            error('statgen:sumstats', '%s: missing required column: %s', path, required{i});
        end
    end

    chr_col = ensure_cellstr_col_(get_col_(tbl, var_names, 'chr'));
    bp_raw = ensure_cellstr_col_(get_col_(tbl, var_names, 'bp'));
    a1_col = ensure_cellstr_col_(get_col_(tbl, var_names, 'a1'));
    a2_col = ensure_cellstr_col_(get_col_(tbl, var_names, 'a2'));
    z_raw = ensure_cellstr_col_(get_col_(tbl, var_names, 'z'));
    n_raw = ensure_cellstr_col_(get_col_(tbl, var_names, 'n'));

    bad_chr = cellfun('isempty', chr_col);
    if any(bad_chr), error('statgen:sumstats', '%s: row %d: chr must be non-empty', path, find(bad_chr, 1, 'first') + 1); end
    bad_a1 = cellfun('isempty', a1_col);
    if any(bad_a1), error('statgen:sumstats', '%s: row %d: a1 must be non-empty', path, find(bad_a1, 1, 'first') + 1); end
    bad_a2 = cellfun('isempty', a2_col);
    if any(bad_a2), error('statgen:sumstats', '%s: row %d: a2 must be non-empty', path, find(bad_a2, 1, 'first') + 1); end

    bp_num = str2double(bp_raw);
    bad_bp = isnan(bp_num) | (bp_num ~= floor(bp_num));
    if any(bad_bp)
        i = find(bad_bp, 1, 'first');
        error('statgen:sumstats', '%s: row %d: bp is not an integer: %s', path, i + 1, bp_raw{i});
    end

    z_num = str2double(z_raw);
    bad_z = ~isfinite(z_num);
    if any(bad_z)
        i = find(bad_z, 1, 'first');
        error('statgen:sumstats', '%s: row %d: z must be finite numeric: %s', path, i + 1, z_raw{i});
    end

    n_num = str2double(n_raw);
    bad_n = ~isfinite(n_num);
    if any(bad_n)
        i = find(bad_n, 1, 'first');
        error('statgen:sumstats', '%s: row %d: n must be finite numeric: %s', path, i + 1, n_raw{i});
    end

    optional_map = struct('p', [], 'beta', [], 'se', [], 'eaf', [], 'info', []);
    optional_names = fieldnames(optional_map);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        if any(strcmp(var_names, nm))
            raw = ensure_cellstr_col_(get_col_(tbl, var_names, nm));
            optional_map.(nm) = str2double(raw);
        end
    end

    src_key = make_key_(chr_col, bp_num, a1_col, a2_col);

    ref_chr = reference.chr;
    ref_bp = reference.bp;
    ref_a1 = reference.a1;
    ref_a2 = reference.a2;
    ref_key = make_key_(ref_chr, ref_bp, ref_a1, ref_a2);

    [has_match, loc] = ismember(ref_key, src_key);

    aligned_z = nan(numel(ref_key), 1);
    aligned_n = nan(numel(ref_key), 1);
    aligned_z(has_match) = z_num(loc(has_match));
    aligned_n(has_match) = n_num(loc(has_match));

    if isempty(optional_map.p)
        aligned_logp = nan(numel(ref_key), 1);
    else
        p_aligned = nan(numel(ref_key), 1);
        p_aligned(has_match) = optional_map.p(loc(has_match));
        aligned_logp = nan(numel(ref_key), 1);
        finite_mask = isfinite(p_aligned);
        in_range = finite_mask & p_aligned >= 0 & p_aligned <= 1;
        zero_mask = in_range & p_aligned == 0;
        pos_mask = in_range & p_aligned > 0;
        aligned_logp(zero_mask) = inf;
        aligned_logp(pos_mask) = -log10(p_aligned(pos_mask));
    end

    aligned_optional = struct('beta', [], 'se', [], 'eaf', [], 'info', []);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        vals = optional_map.(nm);
        if isempty(vals)
            aligned_optional.(nm) = [];
        else
            x = nan(numel(ref_key), 1);
            x(has_match) = vals(loc(has_match));
            aligned_optional.(nm) = x;
        end
    end

    shards = cell(numel(reference.shards), 1);
    for i = 1:numel(reference.shards)
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        beta_vec = pick_optional_(aligned_optional.beta, ix);
        se_vec = pick_optional_(aligned_optional.se, ix);
        eaf_vec = pick_optional_(aligned_optional.eaf, ix);
        info_vec = pick_optional_(aligned_optional.info, ix);

        shards{i} = statgen.SumstatsShard( ...
            s_ref.label, s_ref.checksum, ...
            aligned_z(ix), aligned_n(ix), aligned_logp(ix), ...
            beta_vec, se_vec, eaf_vec, info_vec);
    end

    sumstats = statgen.Sumstats(shards);
end

function [tbl, cleanup_fn] = parse_sumstats_table_(path)
    cleanup_fn = @() [];
    actual_path = path;
    if ends_with_(path, '.gz')
        tmpdir = tempname;
        mkdir(tmpdir);
        gunzip(path, tmpdir);
        [~, base, ~] = fileparts(path);
        actual_path = fullfile(tmpdir, base);
        cleanup_fn = @() cleanup_tmpdir_(tmpdir);
    end

    if exist('readtable', 'file') == 2
        tbl = readtable(actual_path, ...
            'FileType', 'text', ...
            'Delimiter', '\t', ...
            'ReadVariableNames', true, ...
            'TextType', 'char');
        return
    end

    fid = fopen(actual_path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open sumstats file: %s', path);
    end
    closer = onCleanup(@() fclose(fid));
    hdr = fgetl(fid);
    if ~ischar(hdr)
        error('statgen:sumstats', 'sumstats file is empty: %s', path);
    end
    names = strsplit(hdr, '\t');
    fmt = repmat('%s', 1, numel(names));
    cols = textscan(fid, fmt, 'Delimiter', '\t', 'Whitespace', '', 'MultipleDelimsAsOne', false, 'ReturnOnError', false);
    clear closer;
    S = struct();
    S._var_names = names;
    for i = 1:numel(names)
        S.(names{i}) = ensure_cellstr_col_(cols{i});
    end
    tbl = S;
end

function col = get_col_(tbl, var_names, name)
    idx = find(strcmp(var_names, name), 1, 'first');
    if isempty(idx)
        error('statgen:sumstats', 'Missing required column: %s', name);
    end
    if is_table_like_(tbl)
        col = tbl{:, idx};
    else
        field_name = tbl._var_names{idx};
        col = tbl.(field_name);
    end
end

function tf = is_table_like_(x)
    tf = false;
    if isstruct(x) && isfield(x, '_var_names')
        return
    end
    if exist('istable', 'file') == 2
        tf = istable(x);
    else
        tf = isa(x, 'table');
    end
end

function key = make_key_(chr_col, bp_col, a1_col, a2_col)
    bp_str = cellstr(num2str(double(bp_col(:)), '%d'));
    key = strcat(ensure_cellstr_col_(chr_col), {':'}, bp_str, {':'}, ensure_cellstr_col_(a1_col), {':'}, ensure_cellstr_col_(a2_col));
end

function out = pick_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end

function cleanup_tmpdir_(tmpdir)
    if exist(tmpdir, 'dir') == 7
        try
            rmdir(tmpdir, 's');
        catch
        end
    end
end

function tf = ends_with_(s, suffix)
    n = length(suffix);
    tf = length(s) >= n && strcmp(s(end - n + 1:end), suffix);
end

function out = ensure_cellstr_col_(x)
    if iscell(x)
        out = x(:);
    elseif isstring(x)
        out = cellstr(x(:));
    elseif ischar(x)
        out = strtrim(cellstr(x));
    elseif isnumeric(x) || islogical(x)
        out = arrayfun(@(v) num2str(v, '%.15g'), x(:), 'UniformOutput', false);
    else
        out = strtrim(cellstr(x));
    end
end
