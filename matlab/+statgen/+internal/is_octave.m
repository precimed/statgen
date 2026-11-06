function tf = is_octave()
% Return true when running under GNU Octave.
    persistent cached;
    if isempty(cached)
        cached = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    end
    tf = cached;
end
