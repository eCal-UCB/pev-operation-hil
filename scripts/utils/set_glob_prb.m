function varargout = set_glob_prb(var)
global prb
prb = var;

if nargout == 0
    % do nothing
elseif nargout == 1
    varargout = {};
    varargout{1} = var;
else
    error('[ ERROR] number of output must be either 0 or 1');
end
end