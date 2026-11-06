function enabled = parse_allow_indel_swaps(varargin)
%PARSE_ALLOW_INDEL_SWAPS Parse the shared source-loader name-value option.
    enabled = false;
    if isempty(varargin), return; end
    if numel(varargin) ~= 2 || ...
            ~(ischar(varargin{1}) || isstring(varargin{1})) || ...
            ~strcmpi(char(varargin{1}), 'allow_indel_swaps')
        error('statgen:arg', 'Expected name-value option ''allow_indel_swaps'', true or false');
    end
    if ~islogical(varargin{2}) || ~isscalar(varargin{2})
        error('statgen:arg', 'allow_indel_swaps must be a logical scalar');
    end
    enabled = varargin{2};
end
