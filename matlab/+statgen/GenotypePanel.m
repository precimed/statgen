classdef GenotypePanel
%STATGEN.GENOTYPEPANEL Reference-aligned PLINK genotype metadata.
%
%   genotype = statgen.load_genotype(bfile_prefix, reference)
%   genotype = statgen.load_genotype_cache(path)
%
% A GenotypePanel stores sample metadata and SNP mappings aligned to a
% ReferencePanel. Genotype calls are read lazily from PLINK BED files with
% fetch_genotypes or fetch_genotypes_int8. Calls count reference a1 (ALT).
% is_allele_swapped marks swaps and is false for exact matches and absent variants.
%
% FAM rows define the panel-level sample axis. For sharded genotype input,
% loaded autosomal shards must have identical FAM rows and define that axis.
% chrX is the exception: its FAM may be a subset of the panel sample axis.
% genotype.is_subject_present('X') reports which panel-level subjects are
% present in the chrX source; fetched chrX genotypes are expanded back to the
% panel sample order, with missing calls for absent chrX subjects.
%
% Common properties:
%   num_snp      Number of reference-axis SNPs.
%   num_sample   Number of samples.
%   fid, iid     Sample identifiers from FAM metadata.
%   is_present   Logical vector marking SNPs present in the genotype source.
%   is_allele_swapped  Logical vector marking swaps into reference orientation.
%   ploidy_male, ploidy_female  Per-SNP ploidy metadata.
%
% Common methods:
%   is_subject_present   Return sample-presence mask for one shard.
%   select_shards        Restrict the panel to selected shards.
%   fetch_genotypes      Read selected genotype calls as double.
%   fetch_genotypes_int8 Read selected genotype calls as int8.
%   save_cache           Save metadata to a MATLAB .mat cache.
%
% See also statgen.load_genotype, statgen.load_genotype_cache,
% statgen.save_genotype_cache.
    properties (SetAccess = private)
        num_snp
        num_sample
        shard_offsets
        shards
        fid
        iid
        father_id
        mother_id
        sex
        is_male
        is_female
        source_layout
    end
    properties (Dependent)
        is_present
        is_allele_swapped
        ploidy_male
        ploidy_female
        source_row0
    end

    methods
        function obj = GenotypePanel(shards_cell, fid, iid, father_id, mother_id, sex, source_layout)
            if nargin == 0, return; end
            source_layout = char(source_layout);
            if ~strcmp(source_layout, 'non_sharded') && ~strcmp(source_layout, 'sharded')
                error('statgen:genotype', 'source_layout must be ''non_sharded'' or ''sharded''');
            end
            if isempty(shards_cell)
                error('statgen:genotype', 'GenotypePanel requires at least one shard');
            end
            obj.source_layout = source_layout;
            obj.shards = shards_cell;
            obj.fid = statgen.internal.ensure_cell_col(fid);
            obj.iid = statgen.internal.ensure_cell_col(iid);
            obj.father_id = statgen.internal.ensure_cell_col(father_id);
            obj.mother_id = statgen.internal.ensure_cell_col(mother_id);
            obj.sex = double(sex(:));
            obj.num_sample = numel(obj.fid);
            if numel(obj.iid) ~= obj.num_sample || numel(obj.father_id) ~= obj.num_sample || numel(obj.mother_id) ~= obj.num_sample || numel(obj.sex) ~= obj.num_sample
                error('statgen:genotype', 'GenotypePanel FAM vector lengths must match');
            end
            if any(obj.sex ~= floor(obj.sex) | ~ismember(obj.sex, [0; 1; 2]))
                error('statgen:genotype', 'GenotypePanel sex values must be 0, 1, or 2');
            end
            obj.is_male = obj.sex == 1;
            obj.is_female = obj.sex == 2;

            n_shards = numel(shards_cell);
            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});
            total = 0;
            for i = 1:n_shards
                s = shards_cell{i};
                if numel(s.subject_present) ~= obj.num_sample
                    error('statgen:genotype', ...
                        'Genotype shard %s: subject_present length does not match panel sample axis', ...
                        s.label);
                end
                offsets(i).shard_label = s.label;
                offsets(i).start0 = total;
                offsets(i).stop0 = total + s.num_snp;
                total = total + s.num_snp;
            end
            obj.num_snp = total;
            obj.shard_offsets = offsets;
        end

        function display(obj)
            name = inputname(1);
            if ~isempty(name)
                fprintf('%s =\n\n', name);
            end
            disp(obj);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  statgen.GenotypePanel array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            fprintf('  statgen.GenotypePanel object\n\n');
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    num_sample: %d\n', obj.num_sample);
            fprintf('    source_layout: %s\n', obj.source_layout);
            fprintf('    shards: %d\n', numel(obj.shards));
            fprintf('    shard_labels: %s\n', shard_labels_string_(obj.shards));
            fprintf('    present_snps: %d\n', count_present_(obj.shards));
            fprintf('    samples: male=%d, female=%d, unknown=%d\n', ...
                sum(obj.is_male), sum(obj.is_female), sum(obj.sex == 0));
            fprintf('    fid, iid: %d-by-1 cell\n', obj.num_sample);
            fprintf('    genotype_calls: lazy BED fetch via fetch_genotypes/fetch_genotypes_int8\n');
        end

        function out = get.is_present(obj)
            out = concat_required_(obj.shards, 'is_present');
        end

        function out = get.is_allele_swapped(obj)
            out = concat_required_(obj.shards, 'is_allele_swapped');
        end

        function out = get.ploidy_male(obj)
            out = concat_required_(obj.shards, 'ploidy_male');
        end

        function out = get.ploidy_female(obj)
            out = concat_required_(obj.shards, 'ploidy_female');
        end

        function out = get.source_row0(obj)
            out = concat_required_(obj.shards, 'source_row0');
        end

        function out = is_subject_present(obj, shard)
        %IS_SUBJECT_PRESENT Return sample-presence mask for one shard.
        %
        %   mask = genotype.is_subject_present(shard)
        %
        % Returns a num_sample-by-1 logical vector for the requested shard
        % label, for example '21' or 'X'.
            label = char(shard);
            for i = 1:numel(obj.shards)
                if strcmp(obj.shards{i}.label, label)
                    out = obj.shards{i}.subject_present;
                    return
                end
            end
            error('statgen:genotype', ...
                'GenotypePanel.is_subject_present: requested shard %s is not loaded', label);
        end

        function out = select_shards(obj, shards)
        %SELECT_SHARDS Return genotype metadata restricted to selected shards.
        %
        %   out = genotype.select_shards(shards)
        %
        % shards is a cell array or string array of canonical shard labels, for
        % example {'21', '22'}. The returned GenotypePanel preserves the
        % requested shard order.
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.internal.validate_requested_shards( ...
                shards, available, 'GenotypePanel.select_shards');

            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.GenotypePanel(out_shards, obj.fid, obj.iid, ...
                obj.father_id, obj.mother_id, obj.sex, obj.source_layout);
        end

        function out = fetch_genotypes_int8(obj, snp_indices, varargin)
        %FETCH_GENOTYPES_INT8 Read selected genotype calls as int8.
        %
        %   G = genotype.fetch_genotypes_int8(snp_indices)
        %   G = genotype.fetch_genotypes_int8(snp_indices, bed_path)
        %
        % Returns a num_sample-by-numel(snp_indices) matrix. Calls are encoded
        % as reference-a1 counts 0, 1, or 2; missing calls are -1. snp_indices are
        % one-based panel SNP indices in MATLAB. bed_path overrides cached
        % source BED paths and may contain '@' for sharded panels. Requesting a
        % reference SNP absent from the genotype source is an error.
            bed_path = parse_bed_path_(varargin{:});
            indices0 = normalize_snp_indices_(snp_indices, obj.num_snp);
            out = int8(-ones(obj.num_sample, numel(indices0)));
            if isempty(indices0)
                return
            end

            addressed = addressed_shards_(obj, indices0);
            bed_paths = resolve_bed_paths_(obj, addressed, bed_path);
            validate_requested_present_(obj, addressed);

            for i = 1:numel(addressed)
                a = addressed{i};
                shard = obj.shards{a.shard_index};
                effective_bed = bed_paths{i};
                statgen.internal.bfile_validate_bed( ...
                    effective_bed, shard.source_num_sample, shard.source_num_snp);
                source_rows0 = shard.source_row0(a.local_indices0 + 1);
                source_geno = statgen.internal.bfile_read_bed_rows_int8( ...
                    effective_bed, source_rows0, shard.source_num_sample);
                swapped_columns = find(shard.is_allele_swapped(a.local_indices0 + 1));
                if ~isempty(swapped_columns)
                    % Keep lookup indices int8 to avoid a double-sized block temporary.
                    flip_lut = int8([-1; 2; 1; 0]);
                    calls = source_geno(:, swapped_columns);
                    source_geno(:, swapped_columns) = reshape(flip_lut(calls(:) + int8(2)), size(calls));
                end
                panel_rows = find(shard.subject_present);
                if ~isempty(panel_rows)
                    source_subject_rows = shard.source_subject_row0(panel_rows) + 1;
                    out(panel_rows, a.columns) = source_geno(source_subject_rows, :);
                end
            end
        end

        function out = fetch_genotypes(obj, snp_indices, varargin)
        %FETCH_GENOTYPES Read selected genotype calls as double.
        %
        %   G = genotype.fetch_genotypes(snp_indices)
        %   G = genotype.fetch_genotypes(snp_indices, bed_path)
        %   G = genotype.fetch_genotypes(snp_indices, haploid_mode)
        %   G = genotype.fetch_genotypes(snp_indices, bed_path, haploid_mode)
        %
        % Returns a num_sample-by-numel(snp_indices) matrix. Missing calls are
        % returned as NaN. haploid_mode is 'raw' by default. Use
        % 'ploidy_scaled' to map PLINK diploid-style hardcalls onto declared
        % biological ploidy: male chrX with ploidy 1 maps 0/2 to 0/1, while
        % diploid calls stay 0/1/2.
        %
        % snp_indices are one-based panel SNP indices in MATLAB. bed_path
        % overrides cached source BED paths and may contain '@' for sharded
        % panels. Requesting a reference SNP absent from the genotype source is
        % an error.
        %
        % See also statgen.GenotypePanel.fetch_genotypes_int8.
            [bed_path, haploid_mode] = parse_fetch_genotypes_args_(varargin{:});
            geno_int8 = obj.fetch_genotypes_int8(snp_indices, bed_path);
            out = double(geno_int8);
            out(geno_int8 == int8(-1)) = NaN;
            if strcmp(haploid_mode, 'ploidy_scaled')
                indices0 = normalize_snp_indices_(snp_indices, obj.num_snp);
                out = apply_ploidy_scaled_(obj, out, indices0);
            end
        end

        function save_cache(obj, path, varargin)
        %SAVE_CACHE Save genotype metadata to a MATLAB .mat cache.
        %
        %   genotype.save_cache(path)
        %   genotype.save_cache(path, 'format', format)
        %
        % Saves metadata and source BED paths, not dense genotype calls.
        %
        % Default format is 'v7.3' in MATLAB and 'v7' in Octave.
        %
        % See also statgen.save_genotype_cache, statgen.load_genotype_cache.
            statgen.save_genotype_cache(obj, path, varargin{:});
        end
    end
end

function [bed_path, haploid_mode] = parse_fetch_genotypes_args_(varargin)
    bed_path = [];
    haploid_mode = 'raw';
    if nargin == 0
        return
    elseif nargin == 1
        if is_haploid_mode_(varargin{1})
            haploid_mode = normalize_haploid_mode_(varargin{1});
        else
            bed_path = parse_bed_path_(varargin{1});
        end
    elseif nargin == 2
        bed_path = parse_bed_path_(varargin{1});
        haploid_mode = normalize_haploid_mode_(varargin{2});
    else
        error('statgen:genotype', ...
            'GenotypePanel.fetch_genotypes accepts at most bed_path and haploid_mode optional arguments');
    end
end

function tf = is_haploid_mode_(value)
    if isstring(value)
        tf = isscalar(value) && any(strcmp(char(value), {'raw', 'ploidy_scaled'}));
    elseif ischar(value)
        tf = any(strcmp(value, {'raw', 'ploidy_scaled'}));
    else
        tf = false;
    end
end

function haploid_mode = normalize_haploid_mode_(value)
    if isstring(value)
        if ~isscalar(value)
            error('statgen:genotype', 'GenotypePanel.fetch_genotypes: haploid_mode must be raw or ploidy_scaled');
        end
        value = char(value);
    end
    if ~ischar(value) || ~any(strcmp(value, {'raw', 'ploidy_scaled'}))
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: haploid_mode must be raw or ploidy_scaled');
    end
    haploid_mode = value;
end

function out = apply_ploidy_scaled_(obj, out, indices0)
    if isempty(indices0)
        return
    end
    male_rows = find(obj.sex == 1);
    female_rows = find(obj.sex == 2);
    unknown_rows = find(obj.sex == 0);
    addressed = addressed_shards_(obj, indices0);
    for i = 1:numel(addressed)
        a = addressed{i};
        shard = obj.shards{a.shard_index};
        ploidy_male = shard.ploidy_male(a.local_indices0 + 1);
        ploidy_female = shard.ploidy_female(a.local_indices0 + 1);
        differs = ploidy_male ~= ploidy_female;
        unknown_present = unknown_rows(shard.subject_present(unknown_rows));
        if ~isempty(unknown_present) && any(differs) && ...
                any(isfinite(reshape(out(unknown_present, a.columns(differs)), [], 1)))
            error('statgen:genotype', ...
                ['GenotypePanel.fetch_genotypes: haploid_mode=''ploidy_scaled'' ' ...
                 'requires known FAM sex when male and female ploidy differ']);
        end
        if ~isempty(male_rows)
            out(male_rows, a.columns) = bsxfun(@times, out(male_rows, a.columns), ploidy_male(:)' ./ 2);
        end
        if ~isempty(female_rows)
            out(female_rows, a.columns) = bsxfun(@times, out(female_rows, a.columns), ploidy_female(:)' ./ 2);
        end
        if ~isempty(unknown_rows)
            out(unknown_rows, a.columns) = bsxfun(@times, out(unknown_rows, a.columns), ploidy_male(:)' ./ 2);
        end
    end
end

function out = concat_required_(shards, field_name)
    if isempty(shards)
        out = [];
        return
    end
    vals = cell(numel(shards), 1);
    for i = 1:numel(shards)
        vals{i} = shards{i}.(field_name);
    end
    out = vertcat(vals{:});
end

function out = shard_labels_string_(shards)
    labels = cell(numel(shards), 1);
    for i = 1:numel(shards)
        labels{i} = shards{i}.label;
    end
    out = statgen.internal.display_join_strings(labels);
end

function n = count_present_(shards)
    n = 0;
    for i = 1:numel(shards)
        n = n + sum(shards{i}.is_present);
    end
end

function bed_path = parse_bed_path_(varargin)
    if nargin == 0
        bed_path = [];
    elseif nargin == 1
        bed_path = varargin{1};
    else
        error('statgen:genotype', ...
            'GenotypePanel.fetch_genotypes accepts at most one optional bed_path argument');
    end
    if isstring(bed_path)
        if ~isscalar(bed_path)
            error('statgen:genotype', 'GenotypePanel.fetch_genotypes: bed_path must be a string scalar');
        end
        bed_path = char(bed_path);
    elseif ~isempty(bed_path) && ~ischar(bed_path)
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: bed_path must be a character vector');
    end
end

function indices0 = normalize_snp_indices_(snp_indices, num_snp)
    if isempty(snp_indices)
        indices0 = zeros(0, 1);
        return
    end
    if ~isnumeric(snp_indices)
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: snp_indices must be integer indices');
    end
    indices = double(snp_indices(:));
    if any(~isfinite(indices)) || any(indices ~= floor(indices))
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: snp_indices must be integer indices');
    end
    if any(indices < 1) || any(indices > num_snp)
        bad = indices(find(indices < 1 | indices > num_snp, 1, 'first'));
        error('statgen:genotype', ...
            'GenotypePanel.fetch_genotypes: SNP index %.0f is out of bounds for num_snp=%.0f', ...
            bad, num_snp);
    end
    indices0 = indices - 1;
end

function addressed = addressed_shards_(obj, indices0)
    addressed = {};
    for i = 1:numel(obj.shards)
        offset = obj.shard_offsets(i);
        cols = find(indices0 >= offset.start0 & indices0 < offset.stop0);
        if ~isempty(cols)
            a.shard_index = i;
            a.columns = cols(:)';
            a.local_indices0 = indices0(cols) - offset.start0;
            addressed{end + 1, 1} = a; %#ok<AGROW>
        end
    end
end

function bed_paths = resolve_bed_paths_(obj, addressed, bed_path)
    bed_paths = cell(numel(addressed), 1);
    if isempty(bed_path)
        for i = 1:numel(addressed)
            bed_paths{i} = obj.shards{addressed{i}.shard_index}.bed_path;
        end
        return
    end

    if ~isempty(strfind(bed_path, '@'))
        if strcmp(obj.source_layout, 'non_sharded')
            error('statgen:genotype', ...
                'GenotypePanel.fetch_genotypes: @ override incompatible with non-sharded panel metadata');
        end
        for i = 1:numel(addressed)
            shard = obj.shards{addressed{i}.shard_index};
            bed_paths{i} = strrep(bed_path, '@', shard.label);
        end
        return
    end

    if strcmp(obj.source_layout, 'sharded') && ~isempty(addressed)
        source_num_snp = zeros(numel(addressed), 1);
        source_num_sample = zeros(numel(addressed), 1);
        for i = 1:numel(addressed)
            shard = obj.shards{addressed{i}.shard_index};
            source_num_snp(i) = shard.source_num_snp;
            source_num_sample(i) = shard.source_num_sample;
        end
        if numel(unique(source_num_snp)) ~= 1 || numel(unique(source_num_sample)) ~= 1
            error('statgen:genotype', ...
                ['GenotypePanel.fetch_genotypes: flat bed_path override requires equal ' ...
                 'source_num_snp and source_num_sample across addressed shards']);
        end
    end

    for i = 1:numel(addressed)
        bed_paths{i} = bed_path;
    end
end

function validate_requested_present_(obj, addressed)
    for i = 1:numel(addressed)
        a = addressed{i};
        shard = obj.shards{a.shard_index};
        present = shard.is_present(a.local_indices0 + 1);
        if ~all(present)
            j = find(~present, 1, 'first');
            panel_index = obj.shard_offsets(a.shard_index).start0 + a.local_indices0(j) + 1;
            error('statgen:genotype', ...
                'GenotypePanel.fetch_genotypes: requested SNP %.0f in shard %s is not present in the genotype source', ...
                panel_index, shard.label);
        end
    end
end
