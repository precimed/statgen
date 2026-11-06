function save_cache_mat(path, payload, save_arg)
% Save cache variables, failing if the selected format cannot store them.
    warning_id = 'MATLAB:save:sizeTooBigForMATFile';
    previous = warning('query', warning_id);
    restore_warning = onCleanup(@() warning(previous));
    warning('error', warning_id);
    save(path, '-struct', 'payload', save_arg);
end
