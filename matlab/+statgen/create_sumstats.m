function sumstats = create_sumstats(reference, p, z, n, beta, se, eaf, info)
%CREATE_SUMSTATS Create summary statistics from aligned vectors.
%
%   sumstats = statgen.create_sumstats(reference, p)
%   sumstats = statgen.create_sumstats(reference, p, z, n, beta, se, eaf, info)
%
% p and optional vectors must be aligned to reference and have num_snp
% elements. Optional vectors are z, n, beta, se, eaf, and info. p is
% converted to logpvec in the returned Sumstats object. Use load_sumstats for
% raw TSV input; create_sumstats expects full-panel vectors already in the
% reference panel's order and a1 (ALT) orientation. No allele swaps are applied.
% NaN p marks an absent variant: logpvec and present optional vectors are NaN
% there. Other p values must be finite in [0, 1]; p == 0 becomes Inf in logpvec.
%
% See also statgen.Sumstats, statgen.load_sumstats.
    if nargin < 2
        error('statgen:arg', 'create_sumstats requires reference and p');
    end
    if nargin < 3, z = []; end
    if nargin < 4, n = []; end
    if nargin < 5, beta = []; end
    if nargin < 6, se = []; end
    if nargin < 7, eaf = []; end
    if nargin < 8, info = []; end

    num_snp_value = double(reference.num_snp);
    p_aligned = coerce_pvec_(p, num_snp_value, 'p');
    logp_aligned = statgen.internal.sumstats_derive_logp(p_aligned);

    optional_names = {'z', 'n', 'beta', 'se', 'eaf', 'info'};
    optional_values = {z, n, beta, se, eaf, info};
    aligned_optional = struct();
    missing_p = isnan(p_aligned);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        vals = coerce_optional_vec_(optional_values{i}, num_snp_value, nm);
        if ~isempty(vals), vals(missing_p) = NaN; end
        aligned_optional.(nm) = vals;
    end

    n_shards = numel(reference.shards);
    out_shards = cell(n_shards, 1);
    for i = 1:n_shards
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        out_shards{i} = statgen.SumstatsShard( ...
            s_ref.label, s_ref.checksum, ...
            logp_aligned(ix), ...
            pick_optional_(aligned_optional.z, ix), ...
            pick_optional_(aligned_optional.n, ix), ...
            pick_optional_(aligned_optional.beta, ix), ...
            pick_optional_(aligned_optional.se, ix), ...
            pick_optional_(aligned_optional.eaf, ix), ...
            pick_optional_(aligned_optional.info, ix));
    end
    sumstats = statgen.Sumstats(out_shards);
end

function out = coerce_pvec_(x, n, name)
    out = coerce_optional_vec_(x, n, name);
    if isempty(out)
        error('statgen:sumstats', '%s is required', name);
    end
    bad = ~isnan(out) & (~isfinite(out) | out < 0 | out > 1);
    if any(bad)
        idx = find(bad, 1, 'first');
        error('statgen:sumstats', '%s(%d) must be finite numeric in [0, 1]', name, idx);
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

function out = pick_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end
