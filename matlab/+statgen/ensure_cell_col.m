function out = ensure_cell_col(x)
% Convert input to an n-by-1 cell array of trimmed strings.
    if iscell(x)
        out = x(:);
    elseif isstring(x)
        out = cellstr(x(:));
    elseif ischar(x)
        out = strtrim(cellstr(x));
        out = out(:);
    elseif isnumeric(x) || islogical(x)
        out = arrayfun(@(v) num2str(v, '%.15g'), x(:), 'UniformOutput', false);
    else
        out = strtrim(cellstr(x));
        out = out(:);
    end
end
