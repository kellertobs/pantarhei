function [xplt, zplt, varargout] = domainsubset (boxind, x, z, varargin)

% check if want to plot a subset of domain

if isempty(boxind)
    % do not take a subset
    xplt = x; zplt = z;
    varargout = varargin;
else
    if size(boxind,1)==1
        % if only one set of indices are provided
        xind = boxind(1):boxind(2);
        zind = xind;
        
    elseif size(boxind,1)==2
        % if both x, z, indices provided
        xind = boxind(2,1):boxind(2,2);
        zind = boxind(1,1):boxind(1,2);
        
    end
    
    xplt = x(xind);
    zplt = z(zind);
    
    varargout = cell(length(varargin),1);
    for vi = 1:length(varargin)
        varargout{vi} = varargin{vi}(:,zind,xind,:,:);
    end
        
end
end