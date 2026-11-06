function panel = load_reference(path, shards)
%LOAD_REFERENCE Load a ReferencePanel from BIM metadata.
%
%   panel = statgen.load_reference(path)
%   panel = statgen.load_reference(path, shards)
%
% Loads SNP identifiers, chromosome labels, base-pair positions, and alleles
% from PLINK .bim input. path may be a single non-sharded .bim file or an
% @-sharded template such as 'chr@.bim'. The optional shards argument restricts
% loading to selected canonical shard labels, for example {'21', '22'}.
%
% Supported contig labels are 1-22 and X. Labels are not normalized: inputs
% such as chr1 or chrX are rejected and should be fixed upstream.
% BIM a1 is ALT and a2 is REF; the source loaders orient values to this panel.
%
% See also statgen.ReferencePanel, statgen.load_reference_cache.
    if nargin < 2
        shards = [];
    end

    path = char(path);
    if ~isempty(strfind(path, '@'))
        panel = load_sharded_(path, shards);
    else
        bim = reference_bim_rows_(statgen.internal.bfile_parse_bim(path));
        panel = load_split_by_chr_(bim, shards);
    end
end

% ---------------------------------------------------------------------------

function panel = load_sharded_(path, requested_shards)
    canonical = statgen.internal.canonical_labels();
    available_labels = {};
    available_paths = {};
    for i = 1:numel(canonical)
        label = canonical{i};
        bim_path = strrep(path, '@', label);
        if exist(bim_path, 'file') == 2
            available_labels{end+1} = label; %#ok<AGROW>
            available_paths{end+1} = bim_path; %#ok<AGROW>
        end
    end

    if isempty(available_labels)
        error('statgen:io', 'No BIM shards found matching template: %s', path);
    end

    selected = statgen.internal.validate_requested_shards(requested_shards, available_labels, 'load_reference');
    shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        path_idx = find(strcmp(available_labels, label), 1, 'first');
        bim = reference_bim_rows_(statgen.internal.bfile_parse_bim(available_paths{path_idx}));
        shards{i} = statgen.ReferenceShard(label, bim.chr, bim.snp, bim.bp, bim.a1, bim.a2);
        validate_reciprocal_snvs_(shards{i});
    end
    panel = statgen.ReferencePanel(shards);
end

function panel = load_split_by_chr_(bim, requested_shards)
    canonical = statgen.internal.canonical_labels();
    available = {};
    for i = 1:numel(canonical)
        c = canonical{i};
        if any(strcmp(bim.chr, c))
            available{end+1} = c; %#ok<AGROW>
        end
    end

    selected = statgen.internal.validate_requested_shards(requested_shards, available, 'load_reference');
    shards = cell(numel(selected), 1);
    for ci = 1:numel(selected)
        c = selected{ci};
        mask = strcmp(bim.chr, c);
        shards{ci} = statgen.ReferenceShard(c, ...
            bim.chr(mask), bim.snp(mask), ...
            bim.bp(mask), bim.a1(mask), bim.a2(mask));
        validate_reciprocal_snvs_(shards{ci});
    end

    panel = statgen.ReferencePanel(shards);
end

function out = reference_bim_rows_(bim)
    canonical = statgen.internal.canonical_labels();
    keep = ismember(bim.chr, canonical);
    out.chr = bim.chr(keep);
    out.snp = bim.snp(keep);
    out.bp = bim.bp(keep);
    out.a1 = bim.a1(keep);
    out.a2 = bim.a2(keep);
end

function validate_reciprocal_snvs_(shard)
    snv = shard.is_single_nucleotide_variant;
    bp = shard.bp(snv);
    a1 = shard.a1_hash64(snv);
    a2 = shard.a2_hash64(snv);
    lo = min(a1, a2);
    hi = max(a1, a2);
    order = statgen.internal.variant_key_order(bp, lo, hi);
    duplicate = bp(order(2:end)) == bp(order(1:end-1)) & ...
        lo(order(2:end)) == lo(order(1:end-1)) & hi(order(2:end)) == hi(order(1:end-1));
    if any(duplicate)
        pos = order(find(duplicate, 1, 'first'));
        error('statgen:reference', 'reference shard %s: reciprocal SNV keys are not allowed at bp=%d', shard.label, bp(pos));
    end
end
