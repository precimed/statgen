function order = variant_key_order(bp, a1_hash64, a2_hash64)
%VARIANT_KEY_ORDER Sort numeric keys without converting allele hashes to double.
    order = (1:numel(bp))';
    if isempty(order)
        return
    end
    order = stable_sort_by_(order, a2_hash64);
    order = stable_sort_by_(order, a1_hash64);
    order = stable_sort_by_(order, bp);
end

function order = stable_sort_by_(order, values)
    current_rank = (1:numel(order))';
    sorted_values = values(order);
    [sorted_values, p] = sort(sorted_values, 'ascend');
    order = order(p);
    sorted_rank = current_rank(p);

    same_prev = [false; sorted_values(2:end) == sorted_values(1:end-1)];
    group_id = cumsum(~same_prev);
    [~, q] = sortrows([double(group_id(:)), double(sorted_rank(:))], [1, 2]);
    order = order(q);
end
