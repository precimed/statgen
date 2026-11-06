function [bad, lens] = invalid_dna_allele(alleles, lens, workspace_bytes)
%INVALID_DNA_ALLELE Validate complete strings using bounded length buckets.
    alleles = statgen.internal.ensure_cell_col(alleles);
    if nargin < 2, lens = cellfun('length', alleles); end
    bad = lens == 0;
    if isempty(alleles), return; end
    [~, default_workspace] = statgen.internal.allele_processing_limits();
    if nargin < 3, workspace_bytes = default_workspace; end
    % Conservative allowance for char, masks, positions and numeric temporaries.
    max_chars = max(1, floor(workspace_bytes / 32));
    % Single-base validation needs neither padding masks nor a result merge.
    single_idx = find(lens == 1);
    for first = 1:max_chars:numel(single_idx)
        rows = single_idx(first:min(first + max_chars - 1, numel(single_idx)));
        chars = char(alleles(rows));
        bad(rows) = chars ~= 'A' & chars ~= 'C' & chars ~= 'G' & chars ~= 'T';
    end
    if numel(single_idx) == numel(alleles), return; end
    % PERF: geometrical buckets keep padded work linear in sequence volume.
    for k = 1:floor(log2(max([1; lens(:)])))
        upper = 2^(k + 1) - 1;
        idx = find(lens >= 2^k & lens <= upper);
        capacity = max(1, floor(max_chars / upper));
        for first = 1:capacity:numel(idx)
            batch_idx = idx(first:min(first + capacity - 1, numel(idx)));
            batch_lens = lens(batch_idx);
            width = max(batch_lens);
            for offset = 0:max_chars:width-1
                last = min(width, offset + max_chars);
                if offset == 0 && width <= max_chars
                    chars = char(alleles(batch_idx));
                else
                    % PERF: only oversized alleles require character blocks.
                    pieces = cellfun(@(x) x(offset+1:min(numel(x), last)), ...
                        alleles(batch_idx), 'UniformOutput', false);
                    chars = char(pieces);
                end
                bad(batch_idx) = bad(batch_idx) | padded_bad_(chars, batch_lens - offset);
            end
        end
    end
end

function bad = padded_bad_(chars, lens)
% Reusable bounded padded kernel, including the single-base bucket.
    invalid = chars ~= 'A' & chars ~= 'C' & chars ~= 'G' & chars ~= 'T';
    valid_position = bsxfun(@le, 1:size(chars, 2), lens);
    bad = any(invalid & valid_position, 2);
end
