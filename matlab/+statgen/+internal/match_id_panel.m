function [rows, swapped] = match_id_panel(reference, sources, chromosome_scope)
    % sources: cell array of structs with context, ids, a1, a2, chr fields.
    % rows: 1-based indices into concatenated source axes; 0 denotes absent.
    n = reference.num_snp;
    rows = zeros(n, 1);
    swapped = false(n, 1);
    ref_a1 = reference.a1_hash64;
    ref_a2 = reference.a2_hash64;
    eligible = allele_eligibility_(ref_a1, ref_a2, 'reference panel');
    names = cell(numel(reference.shards) + numel(sources), 1);
    ref_rows = cell(numel(reference.shards), 1);
    ref_chr = zeros(n, 1);
    labels = cell(numel(reference.shards), 1);
    missing = 0; invalid = 0;
    for i = 1:numel(reference.shards)
        shard = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        labels{i} = shard.label;
        ref_chr(ix) = i;
        [local, names{i}, nm, ni] = shard.id_matching_names(eligible(ix));
        ref_rows{i} = local + off.start0;
        missing = missing + nm; invalid = invalid + ni;
    end
    report_names_('reference panel', missing, invalid);
    ri = vertcat(ref_rows{:});
    src_rows = cell(numel(sources), 1);
    src_a1 = cell(numel(sources), 1);
    src_a2 = cell(numel(sources), 1);
    src_chr = cell(numel(sources), 1);
    offsets = zeros(numel(sources) + 1, 1);
    for i = 1:numel(sources)
        source = sources{i};
        count = numel(source.ids);
        scope = true(count, 1);
        if chromosome_scope
            [scope, src_chr{i}] = ismember(source.chr, labels);
        else
            src_chr{i} = zeros(count, 1);
        end
        ok = false(count, 1);
        ok(scope) = allele_eligibility_(source.a1(scope), source.a2(scope), source.context);
        [local, names{numel(ref_rows) + i}, nm, ni] = statgen.internal.id_name_matrix( ...
            source.ids, count, ok, false);
        report_names_(source.context, nm, ni);
        src_rows{i} = local + offsets(i);
        src_a1{i} = source.a1;
        src_a2{i} = source.a2;
        offsets(i + 1) = offsets(i) + count;
    end
    si = vertcat(src_rows{:});
    src_a1 = vertcat(src_a1{:});
    src_a2 = vertcat(src_a2{:});
    src_chr = vertcat(src_chr{:});
    if isempty(ri), return; end
    width = 0;
    for i = 1:numel(names)
        width = max(width, size(names{i}, 2));
    end
    for i = 1:numel(names)
        names{i}(:, (size(names{i}, 2) + 1):width) = ' ';
    end
    % Exact character comparisons occur once, before all numeric key work.
    all_names = vertcat(names{:});
    if isempty(all_names), return; end
    [~, ~, codes] = unique(all_names, 'rows');
    ref_codes = codes(1:numel(ri));
    src_codes = codes(numel(ri) + 1:end);
    hashes = statgen.internal.allele_hash64({'A'; 'C'; 'G'; 'T'});
    [~, r1] = ismember(ref_a1(ri), hashes);
    [~, r2] = ismember(ref_a2(ri), hashes);
    [~, s1] = ismember(src_a1(si), hashes);
    [~, s2] = ismember(src_a2(si), hashes);
    n_codes = max(codes);
    % Packing uses only integer ID, allele-group, and shard codes, never hashes.
    % Bound the largest possible key so double arithmetic remains exact.
    if n_codes > (flintmax - 1) / (2 * (numel(labels) + 1))
        error('statgen:match', 'Too many ID codes for exact numeric matching keys');
    end
    ref_keys = ref_codes * 2 + mod(r1 + r2, 2);
    src_keys = src_codes * 2 + mod(s1 + s2, 2);
    if ~isempty(ri)
        [~, ~, group] = unique(ref_keys);
        counts = accumarray(group, 1);
        duplicate = counts(group) > 1;
        if any(duplicate)
            warning('statgen:match', 'reference panel: excluded %d slots in %d duplicate ID key groups', ...
                sum(duplicate), sum(counts > 1));
        end
        ri = ri(~duplicate);
        ref_keys = ref_keys(~duplicate);
    end
    if isempty(ri) || isempty(si), return; end
    if chromosome_scope
        ref_keys = ref_chr(ri) * (n_codes * 2) + ref_keys;
        src_keys = src_chr(si) * (n_codes * 2) + src_keys;
    end
    [hit, loc] = ismember(src_keys, ref_keys);
    matched = loc(hit);
    % Scalar logical indexing can return 0x0 in MATLAB; accumarray needs a column.
    counts = accumarray(matched(:), 1, [numel(ri), 1]);
    if any(counts > 1)
        slot = find(counts > 1, 1, 'first');
        source_row = si(find(hit & loc == slot, 1, 'first'));
        source_index = find(source_row <= offsets(2:end), 1, 'first');
        marker = sources{source_index}.ids{source_row - offsets(source_index)};
        contexts = cellfun(@(s) s.context, sources, 'UniformOutput', false);
        error('statgen:match', 'Ambiguous duplicate source/reference ID matching key for ID ''%s'' in %s', ...
            marker, strjoin(contexts, ', '));
    end
    rows(ri(loc(hit))) = si(hit);
    slots = find(rows > 0);
    [~, source_first] = ismember(src_a1(rows(slots)), hashes);
    [~, target_first] = ismember(ref_a1(slots), hashes);
    [~, target_other] = ismember(ref_a2(slots), hashes);
    swapped(slots) = source_first == target_other | 5 - source_first == target_other;
    complemented = false(n, 1);
    complemented(slots) = source_first ~= target_first & source_first ~= target_other;
    for i = 1:numel(sources)
        selected = rows > offsets(i) & rows <= offsets(i + 1);
        report_count_(sources{i}.context, sum(selected & swapped), 'matched row(s) requiring inversion');
        report_count_(sources{i}.context, sum(selected & complemented), 'matched row(s) requiring strand complementation');
    end
end

function ok = allele_eligibility_(a1, a2, context)
    [snv, palindromic] = statgen.internal.variant_hash_masks(a1, a2);
    report_count_(context, sum(palindromic), 'palindromic SNV(s)', true);
    report_count_(context, sum(~snv), 'multi-base variant(s)', true);
    ok = snv & ~palindromic & a1 ~= a2;
end

function report_count_(context, count, label, excluded)
    if count == 0, return; end
    if nargin < 4, excluded = false; end
    prefix = '';
    if excluded, prefix = 'excluded '; end
    statgen.internal.verbosity_state('info', sprintf('%s: %s%d %s', context, prefix, count, label));
end

function report_names_(context, missing, invalid)
    if missing > 0
        warning('statgen:match', '%s: excluded %d unusable missing/empty/. ID(s)', context, missing);
    end
    if invalid > 0
        warning('statgen:match', '%s: excluded %d unusable overlength or invalid-character ID(s)', context, invalid);
    end
end
