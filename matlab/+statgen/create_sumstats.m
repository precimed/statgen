function sumstats = create_sumstats(reference, zvec, nvec, pvec, beta_vec, se_vec, eaf_vec, info_vec)
% Create Sumstats from full-panel aligned vectors.
    if nargin < 3
        error('statgen:arg', 'create_sumstats requires reference, zvec, and nvec');
    end
    if nargin < 4, pvec = []; end
    if nargin < 5, beta_vec = []; end
    if nargin < 6, se_vec = []; end
    if nargin < 7, eaf_vec = []; end
    if nargin < 8, info_vec = []; end

    n = double(reference.num_snp);
    z_aligned = coerce_required_vec_(zvec, n, 'zvec');
    n_aligned = coerce_required_vec_(nvec, n, 'nvec');

    if isempty(pvec)
        logp_aligned = nan(n, 1);
    else
        p_aligned = coerce_optional_vec_(pvec, n, 'pvec');
        logp_aligned = derive_logp_(p_aligned);
    end

    beta_aligned = coerce_optional_vec_(beta_vec, n, 'beta_vec');
    se_aligned = coerce_optional_vec_(se_vec, n, 'se_vec');
    eaf_aligned = coerce_optional_vec_(eaf_vec, n, 'eaf_vec');
    info_aligned = coerce_optional_vec_(info_vec, n, 'info_vec');

    n_shards = numel(reference.shards);
    out_shards = cell(n_shards, 1);
    for i = 1:n_shards
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        out_shards{i} = statgen.SumstatsShard( ...
            s_ref.label, s_ref.checksum, ...
            z_aligned(ix), n_aligned(ix), logp_aligned(ix), ...
            pick_optional_(beta_aligned, ix), ...
            pick_optional_(se_aligned, ix), ...
            pick_optional_(eaf_aligned, ix), ...
            pick_optional_(info_aligned, ix));
    end
    sumstats = statgen.Sumstats(out_shards);
end

function out = coerce_required_vec_(x, n, name)
    out = coerce_vec_(x, n, name);
    bad = ~(isfinite(out) | isnan(out));
    if any(bad)
        idx = find(bad, 1, 'first');
        error('statgen:sumstats', '%s(%d) must be finite numeric or NaN', name, idx);
    end
end

function out = coerce_optional_vec_(x, n, name)
    if isempty(x)
        out = [];
        return
    end
    out = coerce_vec_(x, n, name);
    bad = ~(isfinite(out) | isnan(out));
    if any(bad)
        idx = find(bad, 1, 'first');
        error('statgen:sumstats', '%s(%d) must be finite numeric or NaN', name, idx);
    end
end

function out = coerce_vec_(x, n, name)
    try
        out = double(x);
    catch
        error('statgen:sumstats', '%s must be numeric', name);
    end
    if ~(isvector(out) || isempty(out))
        error('statgen:sumstats', '%s must be a vector with length %d', name, n);
    end
    out = out(:);
    if numel(out) ~= n
        error('statgen:sumstats', '%s length mismatch: expected %d, got %d', name, n, numel(out));
    end
end

function out = derive_logp_(pvec)
    out = nan(size(pvec));
    finite_mask = isfinite(pvec);
    in_range = finite_mask & pvec >= 0 & pvec <= 1;
    zero_mask = in_range & pvec == 0;
    pos_mask = in_range & pvec > 0;
    out(zero_mask) = inf;
    out(pos_mask) = -log10(pvec(pos_mask));
end

function out = pick_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end
