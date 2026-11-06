function panel = load_genotype(bfile_prefix, reference, varargin)
%LOAD_GENOTYPE Load PLINK genotype metadata aligned to a ReferencePanel.
%
%   genotype = statgen.load_genotype(bfile_prefix, reference)
%   genotype = statgen.load_genotype(bfile_prefix, reference, 'allow_indel_swaps', true)
%
% SNV swaps are automatic; indel swaps are opt-in. Fetched counts are oriented
% to reference a1 (ALT); is_allele_swapped records the matched orientation.
%
% bfile_prefix may be either a single non-sharded PLINK bed/bim/fam prefix,
% such as 'genotypes/all', or an @-sharded template, such as 'genotypes/chr@'.
% In the sharded form, @ is replaced by each reference shard label. The
% returned GenotypePanel stores metadata and source BED paths; genotype calls
% are read lazily with fetch_genotypes.
%
% For sharded genotype input, autosomal FAM files must match exactly and define
% the panel-level sample axis. The chrX FAM file may be a subset of that axis;
% chrX subjects are matched by fid/iid, and absent panel-level subjects receive
% missing calls when chrX genotypes are fetched.
%
% A .ploidy sidecar may be present next to each bfile prefix, with columns
% male_ploidy and female_ploidy in BIM row order. Missing .ploidy defaults
% autosomal rows to (2, 2). Missing chrX .ploidy defaults chrX rows to (1, 2)
% and warns because PAR/non-PAR mixtures and nonstandard encodings should
% provide explicit ploidy metadata.
%
% See also statgen.GenotypePanel, statgen.load_genotype_cache,
% statgen.save_genotype_cache.
    if nargin < 2 || isempty(reference)
        error('statgen:genotype', 'load_genotype requires bfile_prefix and reference');
    end

    allow_indel_swaps = statgen.internal.parse_allow_indel_swaps(varargin{:});
    ref_shards = validate_reference_(reference);
    prefix = char(bfile_prefix);
    ref_labels = cell(numel(ref_shards), 1);
    for i = 1:numel(ref_shards)
        ref_labels{i} = ref_shards{i}.label;
    end

    if ~isempty(strfind(prefix, '@'))
        source_layout = 'sharded';
        sources = struct();
        for i = 1:numel(ref_labels)
            label = ref_labels{i};
            sources.(matlab.lang.makeValidName(['chr_' label])) = load_source_record_(strrep(prefix, '@', label), label);
        end
    else
        source_layout = 'non_sharded';
        shared = load_source_record_(prefix, []);
        if any(strcmp(ref_labels, 'X')) && any(strcmp(shared.bim.chr, 'X')) && exist([prefix '.ploidy'], 'file') ~= 2
            warning('statgen:genotype', ...
                '%s: chrX genotype source has no .ploidy sidecar; defaulting chrX rows to male/female ploidy (1, 2)', ...
                prefix);
        end
        sources = struct();
        for i = 1:numel(ref_labels)
            label = ref_labels{i};
            sources.(matlab.lang.makeValidName(['chr_' label])) = shared;
        end
    end

    subject_maps = struct();
    if strcmp(source_layout, 'non_sharded')
        panel_fam = get_source_(sources, ref_labels{1}).fam;
        for i = 1:numel(ref_labels)
            label = ref_labels{i};
            source = get_source_(sources, label);
            if ~any(strcmp(source.bim.chr, label))
                warning('statgen:genotype', ...
                    '%s: no source BIM rows for requested reference shard %s; marking shard absent', ...
                    prefix, label);
            end
            [subject_present, source_subject_row0] = all_subjects_present_(source.source_num_sample);
            subject_maps.(matlab.lang.makeValidName(['chr_' label])).subject_present = subject_present;
            subject_maps.(matlab.lang.makeValidName(['chr_' label])).source_subject_row0 = source_subject_row0;
        end
    else
        autosomal_labels = ref_labels(~strcmp(ref_labels, 'X'));
        if ~isempty(autosomal_labels)
            panel_fam = get_source_(sources, autosomal_labels{1}).fam;
            for i = 2:numel(autosomal_labels)
                label = autosomal_labels{i};
                if ~fam_equal_(get_source_(sources, label).fam, panel_fam)
                    error('statgen:genotype', 'load_genotype: autosomal FAM mismatch for shard %s', label);
                end
            end
            for i = 1:numel(autosomal_labels)
                label = autosomal_labels{i};
                source = get_source_(sources, label);
                [subject_present, source_subject_row0] = all_subjects_present_(source.source_num_sample);
                subject_maps.(matlab.lang.makeValidName(['chr_' label])).subject_present = subject_present;
                subject_maps.(matlab.lang.makeValidName(['chr_' label])).source_subject_row0 = source_subject_row0;
            end
            if isfield(sources, matlab.lang.makeValidName('chr_X'))
                [subject_present, source_subject_row0] = map_chrx_subjects_(get_source_(sources, 'X').fam, panel_fam, 'load_genotype');
                subject_maps.(matlab.lang.makeValidName('chr_X')).subject_present = subject_present;
                subject_maps.(matlab.lang.makeValidName('chr_X')).source_subject_row0 = source_subject_row0;
            end
        else
            panel_fam = get_source_(sources, 'X').fam;
            [subject_present, source_subject_row0] = all_subjects_present_(get_source_(sources, 'X').source_num_sample);
            subject_maps.(matlab.lang.makeValidName('chr_X')).subject_present = subject_present;
            subject_maps.(matlab.lang.makeValidName('chr_X')).source_subject_row0 = source_subject_row0;
        end
    end

    out_shards = cell(numel(ref_shards), 1);
    for i = 1:numel(ref_shards)
        label = ref_shards{i}.label;
        m = subject_maps.(matlab.lang.makeValidName(['chr_' label]));
        out_shards{i} = build_shard_(ref_shards{i}, get_source_(sources, label), ...
            m.subject_present, m.source_subject_row0, allow_indel_swaps);
    end

    panel = statgen.GenotypePanel(out_shards, panel_fam.fid, panel_fam.iid, ...
        panel_fam.father_id, panel_fam.mother_id, panel_fam.sex, source_layout);
end

function ref_shards = validate_reference_(reference)
    try
        ref_shards = reference.shards;
    catch
        error('statgen:genotype', 'load_genotype requires a ReferencePanel-like object with shards');
    end
    if isempty(ref_shards)
        error('statgen:genotype', 'load_genotype requires a reference with at least one shard');
    end
    labels = cell(numel(ref_shards), 1);
    for i = 1:numel(ref_shards)
        labels{i} = ref_shards{i}.label;
    end
    statgen.internal.validate_requested_shards(labels, labels, 'load_genotype reference');
end

function source = get_source_(sources, label)
    source = sources.(matlab.lang.makeValidName(['chr_' label]));
end

function source = load_source_record_(prefix, label)
    [bed_path, bim_path, fam_path, ploidy_path] = require_source_paths_(prefix, label);
    source_bim = statgen.internal.bfile_parse_bim(bim_path);
    source_num_snp = numel(source_bim.chr);
    source_bim.line = (1:source_num_snp)';
    source_bim.source_row0 = (0:(source_num_snp - 1))';
    if ~isempty(label)
        wrong = ~strcmp(source_bim.chr, label);
        if any(wrong)
            idx = find(wrong, 1, 'first');
            error('statgen:genotype', ...
                '%s:%d: sharded genotype BIM for %s contains chr %s', ...
                bim_path, source_bim.line(idx), char(label), source_bim.chr{idx});
        end
    end
    source_bim.a1_hash64 = statgen.internal.allele_hash64(source_bim.a1);
    source_bim.a2_hash64 = statgen.internal.allele_hash64(source_bim.a2);

    source_fam = statgen.internal.bfile_parse_fam(fam_path);
    if isequal(label, 'X') && isempty(ploidy_path)
        warning('statgen:genotype', ...
            '%s: chrX genotype source has no .ploidy sidecar; defaulting chrX rows to male/female ploidy (1, 2)', ...
            prefix);
    end
    [ploidy_male, ploidy_female] = statgen.internal.bfile_parse_ploidy( ...
        ploidy_path, source_num_snp, source_bim.chr);
    bed_file_size = statgen.internal.bfile_validate_bed(bed_path, numel(source_fam.fid), source_num_snp);

    source.bed_path = bed_path;
    source.bim_path = bim_path;
    source.bed_file_size = bed_file_size;
    source.source_num_snp = source_num_snp;
    source.source_num_sample = numel(source_fam.fid);
    source.bim = source_bim;
    source.fam = source_fam;
    source.ploidy_male = ploidy_male;
    source.ploidy_female = ploidy_female;
end

function [bed_path, bim_path, fam_path, ploidy_path] = require_source_paths_(prefix, label)
    bed_path = [prefix '.bed'];
    bim_path = [prefix '.bim'];
    fam_path = [prefix '.fam'];
    paths = {bed_path, bim_path, fam_path};
    for i = 1:numel(paths)
        if exist(paths{i}, 'file') ~= 2
            if isempty(label)
                error('statgen:io', 'Missing genotype source: %s', paths{i});
            else
                error('statgen:io', 'Missing genotype source for requested shard %s: %s', char(label), paths{i});
            end
        end
    end
    candidate = [prefix '.ploidy'];
    if exist(candidate, 'file') == 2
        ploidy_path = candidate;
    else
        ploidy_path = [];
    end
end

function ok = fam_equal_(lhs, rhs)
    ok = isequal(lhs.fid, rhs.fid) && isequal(lhs.iid, rhs.iid) && ...
        isequal(lhs.father_id, rhs.father_id) && isequal(lhs.mother_id, rhs.mother_id) && ...
        isequal(lhs.sex, rhs.sex);
end

function [subject_present, source_subject_row0] = all_subjects_present_(source_num_sample)
    subject_present = true(source_num_sample, 1);
    source_subject_row0 = (0:(source_num_sample - 1))';
end

function [subject_present, source_subject_row0] = map_chrx_subjects_(chrx_fam, panel_fam, where)
    num_sample = numel(panel_fam.fid);
    if fam_equal_(chrx_fam, panel_fam)
        [subject_present, source_subject_row0] = all_subjects_present_(num_sample);
        return
    end

    panel_keys = strcat(panel_fam.fid, char(9), panel_fam.iid);
    src_keys = strcat(chrx_fam.fid, char(9), chrx_fam.iid);
    [found, panel_rows] = ismember(src_keys, panel_keys);
    if any(~found)
        idx = find(~found, 1, 'first');
        error('statgen:genotype', ...
            '%s: chrX FAM subject (%s, %s) is absent from the autosomal sample axis', ...
            where, chrx_fam.fid{idx}, chrx_fam.iid{idx});
    end

    bad = ~strcmp(chrx_fam.father_id, panel_fam.father_id(panel_rows)) | ...
        ~strcmp(chrx_fam.mother_id, panel_fam.mother_id(panel_rows)) | ...
        chrx_fam.sex ~= panel_fam.sex(panel_rows);
    if any(bad)
        idx = find(bad, 1, 'first');
        error('statgen:genotype', ...
            '%s: chrX FAM metadata mismatch for subject (%s, %s)', ...
            where, chrx_fam.fid{idx}, chrx_fam.iid{idx});
    end

    subject_present = false(num_sample, 1);
    source_subject_row0 = -ones(num_sample, 1);
    source_rows0 = (0:(numel(chrx_fam.fid) - 1))';
    subject_present(panel_rows) = true;
    source_subject_row0(panel_rows) = source_rows0;
end

function shard = build_shard_(ref_shard, source, subject_present, source_subject_row0, allow_indel_swaps)
    source_mask = strcmp(source.bim.chr, ref_shard.label);
    source_bim = subset_bim_(source.bim, source_mask);
    is_allele_swapped = false(ref_shard.num_snp, 1);
    if isempty(source_bim.bp)
        local_match = zeros(ref_shard.num_snp, 1);
    else
        [local_match, is_allele_swapped, n_strand_ambiguous, n_disabled_indels] = statgen.internal.match_shard_numeric( ...
            ref_shard.bp, ref_shard.a1_hash64, ref_shard.a2_hash64, ...
            source_bim.bp, source_bim.a1_hash64, source_bim.a2_hash64, ...
            ref_shard.label, 'genotype', allow_indel_swaps);
        statgen.internal.report_variant_swaps(source.bim_path, ref_shard.label, sum(is_allele_swapped), n_strand_ambiguous, n_disabled_indels);
    end

    is_present = local_match > 0;
    source_row0 = -ones(ref_shard.num_snp, 1);
    ploidy_male = nan(ref_shard.num_snp, 1);
    ploidy_female = nan(ref_shard.num_snp, 1);
    if any(is_present)
        source_rows = source_bim.source_row0(local_match(is_present));
        source_row0(is_present) = source_rows;
        ploidy_male(is_present) = source.ploidy_male(source_rows + 1);
        ploidy_female(is_present) = source.ploidy_female(source_rows + 1);
    end

    shard = statgen.GenotypeShard(ref_shard.label, source.bed_path, ...
        source.bed_file_size, source.source_num_snp, source.source_num_sample, ...
        source_row0, subject_present, source_subject_row0, is_present, is_allele_swapped, ...
        ploidy_male, ploidy_female, ref_shard.checksum);
end

function out = subset_bim_(bim, mask)
    out.chr = bim.chr(mask);
    out.snp = bim.snp(mask);
    out.cm = bim.cm(mask);
    out.bp = bim.bp(mask);
    out.a1 = bim.a1(mask);
    out.a2 = bim.a2(mask);
    out.line = bim.line(mask);
    out.source_row0 = bim.source_row0(mask);
    out.a1_hash64 = bim.a1_hash64(mask);
    out.a2_hash64 = bim.a2_hash64(mask);
end
