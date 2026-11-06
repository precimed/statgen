function [loc, is_allele_swapped, n_strand_ambiguous, n_disabled_indels] = match_shard_numeric(ref_bp, ref_a1_hash64, ref_a2_hash64, src_bp, src_a1_hash64, src_a2_hash64, label, context, allow_indel_swaps)
%MATCH_SHARD_NUMERIC Exact-first alignment, then swaps from unused source rows.
% loc is 1-based on the supplied source axis, with 0 for absent reference rows.
    if nargin < 8 || isempty(context), context = 'source'; end
    if nargin < 9, allow_indel_swaps = false; end

    % Reference validation belongs to load_reference; cached references are trusted.
    [loc, ~] = lookup_candidates_(ref_bp, ref_a1_hash64, ref_a2_hash64, ...
        src_bp, src_a1_hash64, src_a2_hash64, (1:numel(src_bp))', ...
        true(numel(src_bp), 1), label, context);
    is_allele_swapped = false(numel(ref_bp), 1);
    unused = true(numel(src_bp), 1);
    unused(loc(loc > 0)) = false;
    [snv, strand_ambiguous, indel] = statgen.internal.variant_hash_masks(src_a1_hash64, src_a2_hash64);
    eligible = snv | (indel & allow_indel_swaps);
    source_rows = find(unused & (snv | indel));
    missing = find(loc == 0);
    [fallback, n_disabled_indels] = lookup_candidates_( ...
        ref_bp(missing), ref_a1_hash64(missing), ref_a2_hash64(missing), ...
        src_bp(source_rows), src_a2_hash64(source_rows), src_a1_hash64(source_rows), ...
        source_rows, eligible(source_rows), label, context);
    accepted = missing(fallback > 0);
    loc(accepted) = fallback(fallback > 0);
    is_allele_swapped(accepted) = true;
    n_strand_ambiguous = sum(strand_ambiguous(loc(is_allele_swapped)));
end

function [loc, excluded_hits] = lookup_candidates_(ref_bp, ref_a1, ref_a2, src_bp, src_a1, src_a2, source_rows, eligible, label, context)
    n_ref = numel(ref_bp);
    loc = zeros(n_ref, 1);
    excluded_hits = 0;
    if n_ref == 0 || isempty(src_bp), return; end

    bp = [ref_bp(:); src_bp(:)];
    a1 = [uint64(ref_a1(:)); uint64(src_a1(:))];
    a2 = [uint64(ref_a2(:)); uint64(src_a2(:))];
    is_ref = [true(n_ref, 1); false(numel(src_bp), 1)];
    ref_index = [(1:n_ref)'; zeros(numel(src_bp), 1)];
    src_index = [zeros(n_ref, 1); source_rows(:)];
    accepted_src = [false(n_ref, 1); eligible(:)];
    order = statgen.internal.variant_key_order(bp, a1, a2);
    same_prev = [false; bp(order(2:end)) == bp(order(1:end-1)) & ...
        a1(order(2:end)) == a1(order(1:end-1)) & ...
        a2(order(2:end)) == a2(order(1:end-1))];
    group_id = cumsum(~same_prev);
    n_groups = group_id(end);
    ref_count = accumarray(group_id, double(is_ref(order)), [n_groups, 1], @sum, 0);
    src_count = accumarray(group_id, double(~is_ref(order)), [n_groups, 1], @sum, 0);
    eligible_count = accumarray(group_id, double(accepted_src(order)), [n_groups, 1], @sum, 0);
    % Eligibility is determined by the key, so a group's candidates agree.
    excluded_hits = sum(src_count(ref_count > 0 & eligible_count == 0));
    ambiguous = ref_count > 0 & eligible_count > 1;
    if any(ambiguous)
        g = find(ambiguous, 1, 'first');
        pos = order(find(group_id == g, 1, 'first'));
        error('statgen:match', ...
            'Ambiguous duplicate %s/reference matching key in shard %s: bp=%d, a1_hash64=%u, a2_hash64=%u', ...
            char(context), char(label), round(bp(pos)), a1(pos), a2(pos));
    end
    src_by_group = accumarray(group_id, src_index(order), [n_groups, 1], @max, 0);
    ref_groups = group_id(is_ref(order));
    ref_rows = ref_index(order(is_ref(order)));
    matched = eligible_count(ref_groups) == 1;
    loc(ref_rows(matched)) = src_by_group(ref_groups(matched));
end
