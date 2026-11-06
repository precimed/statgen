function [cap, workspace_bytes] = allele_processing_limits()
% Internal limits; workspace accounts for simultaneous matrices and masks.
    cap = 1000000;
    workspace_bytes = 16 * 1024 * 1024;
end
