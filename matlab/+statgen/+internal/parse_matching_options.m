function [allow_indel_swaps, match_by] = parse_matching_options(varargin)
    allow_indel_swaps = false;
    match_by = 'genomic_position';
    if mod(numel(varargin), 2) ~= 0
        error('statgen:arg', 'Expected matching name-value options');
    end
    seen = {};
    for i = 1:2:numel(varargin)
        name = varargin{i};
        if ~(ischar(name) && isrow(name)) && ~(isstring(name) && isscalar(name))
            error('statgen:arg', 'Expected matching name-value options');
        end
        name = lower(char(name));
        if any(strcmp(seen, name))
            error('statgen:arg', 'Duplicate matching option: %s', name);
        end
        seen{end + 1} = name;
        value = varargin{i + 1};
        switch name
            case 'allow_indel_swaps'
                allow_indel_swaps = statgen.internal.parse_allow_indel_swaps(name, value);
            case 'match_by'
                if ~((ischar(value) && isrow(value)) || (isstring(value) && isscalar(value))) || ...
                        ~any(strcmp(char(value), {'genomic_position', 'id'}))
                    error('statgen:arg', 'match_by must be ''genomic_position'' or ''id''');
                end
                match_by = char(value);
            otherwise
                error('statgen:arg', 'Unknown matching option: %s', name);
        end
    end
    if strcmp(match_by, 'id') && allow_indel_swaps
        error('statgen:arg', 'allow_indel_swaps=true is incompatible with match_by=''id''');
    end
end
