function bim = bfile_parse_bim(path)
%BFILE_PARSE_BIM Parse and validate a PLINK BIM sidecar.
%
% Returns all input rows in input order with fields:
%   chr, snp, cm, bp, a1, a2
% Validation ignores Y/MT for supported-row sort checks but preserves those rows.
    path = char(path);
    [raw_cols, n_rows] = read_bim_tabular_(path);

    if n_rows == 0
        error('statgen:bim', 'BIM file is empty: %s', path);
    end

    bim.chr = raw_cols{1};
    bim.snp = raw_cols{2};
    bim.cm  = raw_cols{3};
    bim.bp  = raw_cols{4};
    bim.a1  = raw_cols{5};
    bim.a2  = raw_cols{6};

    canonical = statgen.internal.canonical_labels();
    statgen.internal.validate_variant_chr_labels(bim.chr, path, 0, 'statgen:bim');

    [a1_lens, a2_lens] = statgen.internal.warn_allele_cap(bim.a1, bim.a2, path, (1:n_rows)');
    [bad_a1, a1_lens] = statgen.internal.invalid_dna_allele(bim.a1, a1_lens);
    [bad_a2, a2_lens] = statgen.internal.invalid_dna_allele(bim.a2, a2_lens);

    bad_allele = (a1_lens == 0) | (a2_lens == 0);
    if any(bad_allele)
        lineno = find(bad_allele, 1, 'first');
        error('statgen:bim', '%s:%d: a1 and a2 must be non-empty', path, lineno);
    end

    if any(bad_a1)
        lineno = find(bad_a1, 1, 'first');
        error('statgen:bim', ...
            '%s:%d: a1 must be uppercase DNA bases (A/C/G/T): %s', path, lineno, bim.a1{lineno});
    end
    if any(bad_a2)
        lineno = find(bad_a2, 1, 'first');
        error('statgen:bim', ...
            '%s:%d: a2 must be uppercase DNA bases (A/C/G/T): %s', path, lineno, bim.a2{lineno});
    end
    same_allele = strcmp(bim.a1, bim.a2);
    if any(same_allele)
        lineno = find(same_allele, 1, 'first');
        error('statgen:bim', '%s:%d: a1 and a2 must differ', path, lineno);
    end

    bad_cm = isnan(bim.cm);
    if any(bad_cm)
        lineno = find(bad_cm, 1, 'first');
        error('statgen:bim', '%s:%d: cm is not a number', path, lineno);
    end

    bad_bp = isnan(bim.bp) | (bim.bp ~= floor(bim.bp));
    if any(bad_bp)
        lineno = find(bad_bp, 1, 'first');
        error('statgen:bim', '%s:%d: bp is not an integer', path, lineno);
    end

    keep = ismember(bim.chr, canonical);
    statgen.internal.bfile_validate_source_sort_order(subset_bim_(bim, keep), path, find(keep));
end

function out = subset_bim_(bim, mask)
    out.chr = bim.chr(mask);
    out.snp = bim.snp(mask);
    out.cm = bim.cm(mask);
    out.bp = bim.bp(mask);
    out.a1 = bim.a1(mask);
    out.a2 = bim.a2(mask);
end

function [cols, n_rows] = read_bim_tabular_(path)
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open BIM file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    raw_cols = textscan(fid, '%s%s%f%f%s%s', ...
        'Delimiter', sprintf(' \t'), ...
        'MultipleDelimsAsOne', true, ...
        'ReturnOnError', false);
    clear cleaner;

    cols = cell(1, 6);
    for c = [1 2 5 6]
        cols{c} = statgen.internal.ensure_cell_col(raw_cols{c});
    end
    for c = [3 4]
        cols{c} = double(raw_cols{c}(:));
    end
    n_rows = numel(cols{1});
    for c = 2:6
        if numel(cols{c}) ~= n_rows
            error('statgen:bim', '%s: expected 6 whitespace-delimited columns', path);
        end
    end
end
