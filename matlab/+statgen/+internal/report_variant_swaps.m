function report_variant_swaps(path, label, swaps, strand_ambiguous, disabled_indels)
%REPORT_VARIANT_SWAPS Report nonzero source-alignment counts at info verbosity.
    if ~strcmp(statgen.get_verbosity(), 'info'), return; end
    if swaps > 0
        fprintf(2, '%s: shard %s: %d successful allele swap(s)\n', path, char(label), swaps);
    end
    if strand_ambiguous > 0
        fprintf(2, '%s: shard %s: %d strand-ambiguous A/T or C/G swap(s)\n', path, char(label), strand_ambiguous);
    end
    if disabled_indels > 0
        fprintf(2, '%s: shard %s: %d indel row(s) unmatched because indel swaps are disabled\n', path, char(label), disabled_indels);
    end
end
