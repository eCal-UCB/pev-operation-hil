function varargout = set_glob_par(var)
global par
par = var;

if nargout == 0
    % do nothing
elseif nargout == 1
    varargout = {};
    varargout{1} = var;
else
    error('[ ERROR] number of output must be either 0 or 1');
end
end