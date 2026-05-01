function varargout = statgen_load_sumstats(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_sumstats(varargin{:});
end
