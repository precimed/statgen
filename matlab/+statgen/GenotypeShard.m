classdef GenotypeShard
% Immutable metadata/accessor state for one reference-aligned genotype shard.
    properties (SetAccess = private)
        label
        chr
        num_snp
        bed_path
        bed_file_size
        source_num_snp
        source_num_sample
        source_row0
        subject_present
        source_subject_row0
        is_present
        is_allele_swapped
        ploidy_male
        ploidy_female
        reference_checksum
    end
    methods
        function obj = GenotypeShard(label, bed_path, bed_file_size, source_num_snp, source_num_sample, source_row0, subject_present, source_subject_row0, is_present, is_allele_swapped, ploidy_male, ploidy_female, reference_checksum)
            if nargin == 0, return; end

            obj.label = char(label);
            obj.chr = obj.label;
            obj.bed_path = char(bed_path);
            obj.bed_file_size = double(bed_file_size);
            obj.source_num_snp = double(source_num_snp);
            obj.source_num_sample = double(source_num_sample);
            if obj.source_num_snp < 0 || obj.source_num_sample < 0 || obj.source_num_snp ~= floor(obj.source_num_snp) || obj.source_num_sample ~= floor(obj.source_num_sample)
                error('statgen:genotype', 'GenotypeShard %s: source_num_snp and source_num_sample must be non-negative integers', obj.label);
            end
            expected = statgen.internal.bfile_expected_bed_size(obj.source_num_sample, obj.source_num_snp);
            if obj.bed_file_size ~= expected
                error('statgen:genotype', ...
                    'Invalid genotype metadata for shard %s: bed_file_size %.0f does not match expected %.0f', ...
                    obj.label, obj.bed_file_size, expected);
            end

            obj.source_row0 = double(source_row0(:));
            obj.is_present = logical(is_present(:));
            obj.is_allele_swapped = logical(is_allele_swapped(:));
            obj.ploidy_male = double(ploidy_male(:));
            obj.ploidy_female = double(ploidy_female(:));
            obj.num_snp = numel(obj.is_present);
            if numel(obj.is_allele_swapped) ~= obj.num_snp || numel(obj.source_row0) ~= obj.num_snp || numel(obj.ploidy_male) ~= obj.num_snp || numel(obj.ploidy_female) ~= obj.num_snp
                error('statgen:genotype', 'GenotypeShard %s: SNP-axis vector lengths must match', obj.label);
            end
            if ~isequal(obj.is_present, obj.source_row0 >= 0)
                error('statgen:genotype', 'GenotypeShard %s: is_present must match source_row0 >= 0', obj.label);
            end
            bad_source = obj.source_row0 < -1 | obj.source_row0 >= obj.source_num_snp | obj.source_row0 ~= floor(obj.source_row0);
            if any(bad_source)
                error('statgen:genotype', 'GenotypeShard %s: source_row0 out of source BIM bounds', obj.label);
            end
            if any(obj.is_allele_swapped & ~obj.is_present)
                error('statgen:genotype', 'GenotypeShard %s: is_allele_swapped requires is_present', obj.label);
            end
            absent = ~obj.is_present;
            if any(~isnan(obj.ploidy_male(absent))) || any(~isnan(obj.ploidy_female(absent)))
                error('statgen:genotype', 'GenotypeShard %s: absent SNPs must have NaN ploidy', obj.label);
            end
            present_ploidy = [obj.ploidy_male(obj.is_present); obj.ploidy_female(obj.is_present)];
            if any(~ismember(present_ploidy, [0; 1; 2]))
                error('statgen:genotype', 'GenotypeShard %s: present ploidy values must be 0, 1, or 2', obj.label);
            end

            obj.subject_present = logical(subject_present(:));
            obj.source_subject_row0 = double(source_subject_row0(:));
            if numel(obj.subject_present) ~= numel(obj.source_subject_row0)
                error('statgen:genotype', 'GenotypeShard %s: sample-axis vector lengths must match', obj.label);
            end
            if ~isequal(obj.subject_present, obj.source_subject_row0 >= 0)
                error('statgen:genotype', 'GenotypeShard %s: subject_present must match source_subject_row0 >= 0', obj.label);
            end
            bad_subject = obj.source_subject_row0 < -1 | obj.source_subject_row0 >= obj.source_num_sample | obj.source_subject_row0 ~= floor(obj.source_subject_row0);
            if any(bad_subject)
                error('statgen:genotype', 'GenotypeShard %s: source_subject_row0 out of source FAM bounds', obj.label);
            end

            obj.reference_checksum = char(reference_checksum);
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
                fprintf('  statgen.GenotypeShard array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            fprintf('  statgen.GenotypeShard object\n\n');
            fprintf('    label: %s\n', obj.label);
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    present_snps: %d\n', sum(obj.is_present));
            fprintf('    source_num_snp: %d\n', obj.source_num_snp);
            fprintf('    source_num_sample: %d\n', obj.source_num_sample);
            fprintf('    present_subjects: %d\n', sum(obj.subject_present));
            fprintf('    bed_path: %s\n', obj.bed_path);
            fprintf('    reference_checksum: %s\n', obj.reference_checksum);
        end
    end
end
