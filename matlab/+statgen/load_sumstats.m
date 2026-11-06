function sumstats = load_sumstats(path, reference, varargin)
%LOAD_SUMSTATS Load summary statistics aligned to a ReferencePanel.
%
%   sumstats = statgen.load_sumstats(path, reference)
%   sumstats = statgen.load_sumstats(path, reference, 'allow_indel_swaps', true)
%
% SNV swaps are automatic; indel swaps are opt-in. Values are oriented to the
% reference panel's a1 (ALT) allele. See the guide for details:
% https://github.com/precimed/statgen/blob/HEAD/docs/ALLELE_MATCHING.md
%
% Loads one TSV or TSV.GZ summary-statistics file for one trait or source; this
% input is not sharded. Rows are aligned to reference SNPs by chromosome,
% base-pair position, and alleles. The output is a Sumstats object with
% SNP-axis vectors in reference order.
%
% Required columns: chr, bp, a1, a2, p. Optional columns include z, n, beta,
% se, eaf, and info.
% Column names are case-insensitive; POS, EffectAllele, and OtherAllele are
% accepted aliases for bp, a1, and a2. SNP may be present but is not used for
% matching.
%
% Rows with missing p (empty, NA, or NaN) are validated, then dropped before
% matching. At info verbosity the loader reports the number dropped.
% Rows are projected into reference order. Variants absent from the source have
% NaN logp values and is_present false. p == 0 is allowed and becomes Inf in
% logpvec. Optional columns absent from the source return [] in MATLAB.
%
% For large .tsv.gz inputs, set STATGEN_SCRATCH env variable to control where 
% temporary decompressed files are created.
%
% See also statgen.Sumstats, statgen.load_sumstats_cache,
% statgen.save_sumstats_cache.
    if nargin < 2
        error('statgen:arg', 'load_sumstats requires path and reference');
    end
    allow_indel_swaps = statgen.internal.parse_allow_indel_swaps(varargin{:});
    path = char(path);
    if isempty(path)
        error('statgen:arg', 'path must be non-empty');
    end

    [tbl, cleanup_fn] = parse_sumstats_table_(path);
    cleaner = onCleanup(cleanup_fn);
    sumstats = build_sumstats_(tbl, reference, path, allow_indel_swaps);
    clear cleaner;
end

function sumstats = build_sumstats_(tbl, reference, path, allow_indel_swaps)
    if is_table_like_(tbl)
        var_names = canonicalize_var_names_(tbl.Properties.VariableNames, path);
    else
        var_names = canonicalize_var_names_(tbl.statgen_var_names__(:)', path);
    end
    required = {'chr', 'bp', 'a1', 'a2', 'p'};
    for i = 1:numel(required)
        if ~any(strcmp(var_names, required{i}))
            error('statgen:sumstats', '%s: missing required column: %s', path, required{i});
        end
    end

    chr_col = ensure_cellstr_col_(get_col_(tbl, var_names, 'chr'));
    bp_num  = require_numeric_col_(get_col_(tbl, var_names, 'bp'));
    a1_col  = ensure_cellstr_col_(get_col_(tbl, var_names, 'a1'));
    a2_col  = ensure_cellstr_col_(get_col_(tbl, var_names, 'a2'));
    p_num   = require_numeric_col_(get_col_(tbl, var_names, 'p'));
    source_line = (2:(numel(chr_col) + 1))';

    keep_chr = statgen.internal.validate_variant_chr_labels(chr_col, path, 1, 'statgen:sumstats');
    chr_col = chr_col(keep_chr);
    bp_num = bp_num(keep_chr);
    a1_col = a1_col(keep_chr);
    a2_col = a2_col(keep_chr);
    p_num = p_num(keep_chr);
    source_line = source_line(keep_chr);

    [a1_lens, a2_lens] = statgen.internal.warn_allele_cap(a1_col, a2_col, path, source_line);
    [bad_a1_syntax, a1_lens] = statgen.internal.invalid_dna_allele(a1_col, a1_lens);
    [bad_a2_syntax, a2_lens] = statgen.internal.invalid_dna_allele(a2_col, a2_lens);

    bad_a1 = a1_lens == 0;
    if any(bad_a1), error('statgen:sumstats', '%s: row %d: a1 must be non-empty', path, source_line(find(bad_a1, 1, 'first'))); end
    bad_a2 = a2_lens == 0;
    if any(bad_a2), error('statgen:sumstats', '%s: row %d: a2 must be non-empty', path, source_line(find(bad_a2, 1, 'first'))); end
    if any(bad_a1_syntax)
        i = find(bad_a1_syntax, 1, 'first');
        error('statgen:sumstats', '%s: row %d: a1 must be uppercase DNA bases (A/C/G/T): %s', path, source_line(i), a1_col{i});
    end
    if any(bad_a2_syntax)
        i = find(bad_a2_syntax, 1, 'first');
        error('statgen:sumstats', '%s: row %d: a2 must be uppercase DNA bases (A/C/G/T): %s', path, source_line(i), a2_col{i});
    end

    bad_bp = ~isfinite(bp_num) | (bp_num ~= floor(bp_num));
    if any(bad_bp)
        i = find(bad_bp, 1, 'first');
        error('statgen:sumstats', '%s: row %d: bp is not an integer', path, source_line(i));
    end

    bad_p = ~isnan(p_num) & (~isfinite(p_num) | p_num < 0 | p_num > 1);
    if any(bad_p)
        i = find(bad_p, 1, 'first');
        error('statgen:sumstats', '%s: row %d: p must be finite numeric in [0, 1]', path, source_line(i));
    end

    optional_map = struct('z', [], 'n', [], 'beta', [], 'se', [], 'eaf', [], 'info', []);
    optional_names = fieldnames(optional_map);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        if any(strcmp(var_names, nm))
            vals = require_numeric_col_(get_col_(tbl, var_names, nm));
            vals = vals(keep_chr);
            vals(~isfinite(vals)) = NaN;
            optional_map.(nm) = vals;
        end
    end

    n_ref = double(reference.num_snp);
    src_a1_hash64 = statgen.internal.allele_hash64(a1_col);
    src_a2_hash64 = statgen.internal.allele_hash64(a2_col);
    same_hash = src_a1_hash64 == src_a2_hash64;
    if any(same_hash)
        same_idx = find(same_hash);
        same_exact = strcmp(a1_col(same_idx), a2_col(same_idx));
        if any(same_exact)
            i = same_idx(find(same_exact, 1, 'first'));
            error('statgen:sumstats', '%s: row %d: a1 and a2 must differ', path, source_line(i));
        end
    end

    keep_p = ~isnan(p_num);
    if any(~keep_p)
        statgen.internal.verbosity_state('info', ...
            sprintf('%s: dropped %d row(s) with missing p', path, sum(~keep_p)));
    end
    chr_col = chr_col(keep_p);
    bp_num = bp_num(keep_p);
    p_num = p_num(keep_p);
    src_a1_hash64 = src_a1_hash64(keep_p);
    src_a2_hash64 = src_a2_hash64(keep_p);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        if any(strcmp(var_names, nm))
            vals = optional_map.(nm);
            optional_map.(nm) = vals(keep_p);
        end
    end

    p_aligned = nan(n_ref, 1);

    aligned_optional = struct('z', [], 'n', [], 'beta', [], 'se', [], 'eaf', [], 'info', []);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        if any(strcmp(var_names, nm))
            aligned_optional.(nm) = nan(n_ref, 1);
        else
            aligned_optional.(nm) = [];
        end
    end

    for i = 1:numel(reference.shards)
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        src_idx = find(strcmp(chr_col, s_ref.label));
        [loc, swapped, n_strand_ambiguous, n_disabled_indels] = statgen.internal.match_shard_numeric( ...
            s_ref.bp, s_ref.a1_hash64, s_ref.a2_hash64, ...
            bp_num(src_idx), src_a1_hash64(src_idx), src_a2_hash64(src_idx), ...
            s_ref.label, 'sumstats', allow_indel_swaps);
        statgen.internal.report_variant_swaps(path, s_ref.label, sum(swapped), n_strand_ambiguous, n_disabled_indels);
        has_match = loc > 0;
        matched_ref_ix = ix(has_match);
        matched_src_ix = src_idx(loc(has_match));

        p_aligned(matched_ref_ix) = p_num(matched_src_ix);
        swapped_ix = ix(swapped);
        for j = 1:numel(optional_names)
            nm = optional_names{j};
            if ~isempty(aligned_optional.(nm))
                vals = optional_map.(nm);
                tmp = aligned_optional.(nm);
                tmp(matched_ref_ix) = vals(matched_src_ix);
                if strcmp(nm, 'z') || strcmp(nm, 'beta')
                    tmp(swapped_ix) = -tmp(swapped_ix);
                elseif strcmp(nm, 'eaf')
                    tmp(swapped_ix) = 1 - tmp(swapped_ix);
                end
                aligned_optional.(nm) = tmp;
            end
        end
        if ~any(has_match)
            warning('statgen:sumstats', '%s: no matching source rows for reference shard ''%s''', path, s_ref.label);
        end
    end

    aligned_logp = statgen.internal.sumstats_derive_logp(p_aligned);
    statgen.internal.sumstats_warn_optional_zn_completeness( ...
        aligned_optional.z, aligned_optional.n, aligned_logp, 'load_sumstats');

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
            aligned_logp(ix), ...
            pick_optional_(aligned_optional.z, ix), ...
            pick_optional_(aligned_optional.n, ix), ...
            beta_vec, se_vec, eaf_vec, info_vec);
    end

    sumstats = statgen.Sumstats(shards);
end

function var_names = canonicalize_var_names_(raw_names, path)
    var_names = lower(raw_names);
    for i = 1:numel(var_names)
        if strcmp(var_names{i}, 'pos')
            var_names{i} = 'bp';
        elseif strcmp(var_names{i}, 'effectallele')
            var_names{i} = 'a1';
        elseif strcmp(var_names{i}, 'otherallele')
            var_names{i} = 'a2';
        end
    end
    for i = 1:numel(var_names)
        if sum(strcmp(var_names, var_names{i})) > 1
            error('statgen:sumstats', '%s: duplicate columns after column normalization: %s', path, var_names{i});
        end
    end
end

function [tbl, cleanup_fn] = parse_sumstats_table_(path)
    cleanup_fn = @() [];
    actual_path = path;
    cleanup_guard = [];
    if ends_with_(path, '.gz')
        [~, base, ~] = fileparts(path);
        scratch_root = statgen_scratch_root_(path);
        tmpdir = unique_tmpdir_(scratch_root, ['sumstats_gunzip_' base]);
        mkdir(tmpdir);
        cleanup_guard = onCleanup(@() cleanup_tmpdir_(tmpdir));
        gunzip(path, tmpdir);
        actual_path = fullfile(tmpdir, base);
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
    fmt = build_col_formats_(names);
    cols = scan_sumstats_cols_(fid, fmt);
    clear closer;
    validate_col_lengths_(cols, path);
    S = struct();
    S.statgen_var_names__ = names;
    for i = 1:numel(names)
        if isnumeric(cols{i})
            S.(names{i}) = double(cols{i}(:));
        else
            S.(names{i}) = ensure_cellstr_col_(cols{i});
        end
    end
    tbl = S;
    clear cleanup_guard;
end

function col = get_col_(tbl, var_names, name)
    idx = find(strcmp(var_names, name), 1, 'first');
    if isempty(idx)
        error('statgen:sumstats', 'Missing required column: %s', name);
    end
    if is_table_like_(tbl)
        col = tbl{:, idx};
    else
        field_name = tbl.statgen_var_names__{idx};
        col = tbl.(field_name);
    end
end

function tf = is_table_like_(x)
    tf = false;
    if isstruct(x) && isfield(x, 'statgen_var_names__')
        return
    end
    if exist('istable', 'file') == 2
        tf = istable(x);
    else
        tf = isa(x, 'table');
    end
end

function out = pick_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end

function root = statgen_scratch_root_(path)
    env_root = getenv('STATGEN_SCRATCH');
    if ~isempty(env_root)
        root = char(env_root);
    else
        root = tempdir;
    end
    if exist(root, 'dir') ~= 7
        [ok, msg] = mkdir(root);
        if ~ok
            error('statgen:io', 'Cannot create STATGEN scratch directory %s: %s', root, msg);
        end
    end
end

function tmpdir = unique_tmpdir_(root, prefix)
    root = char(root);
    for i = 1:100
        suffix = sprintf('%s_%06d_%06d', prefix, round(1e6 * rand()), i);
        tmpdir = fullfile(root, suffix);
        if exist(tmpdir, 'dir') ~= 7 && exist(tmpdir, 'file') ~= 2
            return
        end
    end
    error('statgen:io', 'Could not allocate unique scratch directory under %s', root);
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

function fmt = build_col_formats_(raw_names)
    numeric = is_numeric_sumstats_col_(raw_names);
    parts = repmat({'%s'}, 1, numel(raw_names));
    for i = find(numeric)
        parts{i} = '%f';
    end
    fmt = strjoin(parts, '');
end

function cols = scan_sumstats_cols_(fid, fmt)
    cols = textscan(fid, fmt, ...
        'Delimiter', '\t', ...
        'Whitespace', '', ...
        'MultipleDelimsAsOne', false, ...
        'ReturnOnError', false, ...
        'EmptyValue', NaN, ...
        'TreatAsEmpty', {'NA', 'na', 'N/A', '.', 'NaN', 'nan', 'Inf', 'inf'});
end

function tf = is_numeric_sumstats_col_(raw_names)
    names = lower(raw_names);
    numeric_names = {'bp', 'pos', 'p', 'z', 'n', 'beta', 'se', 'eaf', 'info'};
    tf = false(1, numel(names));
    for i = 1:numel(names)
        tf(i) = any(strcmp(names{i}, numeric_names));
    end
end

function validate_col_lengths_(cols, path)
    if isempty(cols)
        return
    end
    n = numel(cols{1});
    for i = 2:numel(cols)
        if numel(cols{i}) ~= n
            error('statgen:sumstats', '%s: malformed TSV', path);
        end
    end
end

function out = require_numeric_col_(x)
    if ~isnumeric(x)
        error('statgen:sumstats', 'internal error: expected numeric sumstats column');
    end
    out = double(x(:));
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
