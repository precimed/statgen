function format = default_cache_mat_format()
% Choose a cache MAT-file format supported by the current runtime.
    if statgen.internal.is_octave()
        format = 'v7';
    else
        format = 'v7.3';
    end
end
