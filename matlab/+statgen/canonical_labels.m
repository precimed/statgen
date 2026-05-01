function labels = canonical_labels()
% Canonical contig labels for shard order and validation.
    labels = [arrayfun(@num2str, 1:22, 'UniformOutput', false), {'X'}];
end
