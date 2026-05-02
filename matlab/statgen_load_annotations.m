function varargout = statgen_load_annotations(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_annotations(varargin{:});
end
