function [rows, names, n_missing, n_invalid] = id_name_matrix(values, n, eligible, packed)
    % Convert eligible names only; packed payloads never become per-ID cells.
    if packed
        payload = values(:);
        if n == 0
            starts = zeros(0, 1); lengths = zeros(0, 1);
        else
            breaks = find(payload == char(10));
            starts = [1; breaks + 1];
            lengths = [breaks; numel(payload) + 1] - starts;
        end
        if numel(lengths) ~= n
            error('statgen:cache', 'Invalid reference cache: SNP ID payload length mismatch');
        end
    else
        lengths = cellfun('length', values(:));
    end
    width = min(32, max([0; lengths]));
    missing = lengths == 0;
    candidates = find(eligible & lengths >= 1 & lengths <= 32);
    names = repmat(' ', numel(candidates), width);
    valid = false(numel(candidates), 1);
    placeholders = false(n, 1);
    % PERF: bounded batches cap index workspace at 65536-by-32, including
    % packed IDs. Each batch uses native indexing and character operations.
    for first = 1:65536:numel(candidates)
        ix = first:min(first + 65535, numel(candidates));
        selected = candidates(ix);
        if packed
            columns = 0:(width - 1);
            inside = bsxfun(@lt, columns, lengths(selected));
            indices = bsxfun(@plus, starts(selected), columns);
            block = repmat(' ', numel(selected), width);
            block(inside) = payload(indices(inside));
        else
            block = char(values(selected));
            block(:, (size(block, 2) + 1):width) = ' ';
            inside = bsxfun(@lt, 0:(width - 1), lengths(selected));
        end
        placeholder = lengths(selected) == 1 & block(:, 1) == '.';
        placeholders(selected) = placeholder;
        good = ~any(inside & (block < char(33) | block > char(126)), 2) & ~placeholder;
        % Invalid IDs never enter the exact-comparison matrix.
        names(ix(good), :) = block(good, :);
        valid(ix) = good;
    end
    missing = missing | placeholders;
    n_missing = sum(eligible & missing);
    n_invalid = sum(eligible) - n_missing - sum(valid);
    rows = candidates(valid);
    names = names(valid, :);
end
