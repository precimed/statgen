function [a1_lens, a2_lens] = warn_allele_cap(a1, a2, path, rows)
% Source-boundary diagnostic, once across both allele columns.
    [cap, ~] = statgen.internal.allele_processing_limits();
    a1_lens = cellfun('length', a1);
    a2_lens = cellfun('length', a2);
    affected = [a1_lens(:) > cap, a2_lens(:) > cap];
    count = sum(affected(:));
    if count
        row = find(any(affected, 2), 1);
        col = find(affected(row, :), 1);
        lens = [a1_lens(row), a2_lens(row)];
        warning('statgen:allele_hash', ...
            '%s: %d allele(s) exceed hash cap %d; first affected row %d, column a%d, length %d; hashing uses first %d bases', ...
            path, count, cap, rows(row), col, lens(col), cap);
    end
end
