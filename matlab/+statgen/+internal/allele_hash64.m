function out = allele_hash64(alleles, workspace_bytes)
%ALLELE_HASH64 Exact polynomial hash; source loaders own cap warnings.
    alleles = statgen.internal.ensure_cell_col(alleles);
    full_lens = cellfun('length', alleles);
    [cap, default_workspace] = statgen.internal.allele_processing_limits();
    if nargin < 2, workspace_bytes = default_workspace; end
    lens = min(full_lens, cap);
    out = zeros(numel(alleles), 1, 'uint64');
    if isempty(alleles), return; end
    p = 2147483647;
    bases = [257, 263];
    single = lens == 1;
    base_chars = {'A'; 'C'; 'G'; 'T'};
    for k = 1:4
        x = double(base_chars{k}) + 1;
        h = mod(bases + x, p);
        mask = single & strcmp(alleles, base_chars{k});
        out(mask) = bitshift(uint64(h(1)), 32) + uint64(h(2));
        single(mask) = false;
    end
    idx = find(lens ~= 1 | single);
    if isempty(idx), return; end
    longest = max(lens(idx));
    powers = ones(longest + 1, 2);
    % PERF: exact small-base recurrence over power positions, once per call.
    for j = 1:longest
        powers(j + 1, :) = mod(powers(j, :) .* bases, p);
    end
    % Gather from contiguous vectors rather than repeatedly indexing a column
    % of the two-lane matrix for every byte in each batch.
    power1 = powers(:, 1);
    power2 = powers(:, 2);
    clear powers;
    % Each reduction has at most 8192 terms, keeping double sums exact.
    block_chars = min(8192, max(1, floor(workspace_bytes / 64)));
    % PERF: bucket dispatch and bounded batches replace whole-vector scans.
    for k = 0:floor(log2(max(1, longest)))
        upper = 2^(k + 1) - 1;
        bucket = idx(lens(idx) >= 2^k & lens(idx) <= upper);
        capacity = max(1, floor(workspace_bytes / (64 * min(upper, block_chars))));
        for first = 1:capacity:numel(bucket)
            rows = bucket(first:min(first + capacity - 1, numel(bucket)));
            lengths = lens(rows);
            width = max(lengths);
            batch = alleles(rows);
            % Long alleles are sliced directly into capped byte blocks below;
            % avoid copying a complete million-base prefix for every batch row.
            h = [power1(lengths + 1), power2(lengths + 1)];
            for offset = 0:block_chars:width-1
                last = min(width, offset + block_chars);
                if offset == 0 && width <= block_chars
                    chars = char(batch);
                else
                    pieces = cellfun(@(x) x(offset+1:min(numel(x), last)), batch, 'UniformOutput', false);
                    chars = char(pieces);
                end
                cols = offset + (1:size(chars, 2));
                active = bsxfun(@le, cols, lengths);
                x = (double(uint8(chars)) + 1) .* active;
                positions = max(0, bsxfun(@minus, lengths, cols)) + 1;
                weights = reshape(power1(positions(:)), size(positions));
                h(:, 1) = mod(h(:, 1) + sum(x .* weights, 2), p);
                weights = reshape(power2(positions(:)), size(positions));
                h(:, 2) = mod(h(:, 2) + sum(x .* weights, 2), p);
            end
            out(rows) = bitshift(uint64(h(:, 1)), 32) + uint64(h(:, 2));
        end
    end
end
