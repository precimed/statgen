function varargout = statgen_version(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.version(varargin{:});
end
