function varargout = statgen_load_reference(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_reference(varargin{:});
end
