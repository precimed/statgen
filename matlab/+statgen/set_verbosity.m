function set_verbosity(level)
%SET_VERBOSITY Set statgen runtime verbosity.
%
%   statgen.set_verbosity(level)
%
% level must be one of 'quiet' or 'info'. The default is 'info', which reports
% LD loading progress and nonzero source swap/disabled-indel counts per shard.
% 'quiet' suppresses these reports while preserving validation warnings.
%
% See also statgen.get_verbosity.
    statgen.internal.verbosity_state('set', level);
end
