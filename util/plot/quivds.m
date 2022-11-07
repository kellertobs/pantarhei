function [Xquiv, Zquiv, varargout] = quivds (x, z, opt, scl, varargin)
% downsample and rescale quiver depending on options


% downsample u, w for quiver
xqind = unique(round(linspace(opt.bcind*length(x), (1-opt.bcind)*length(x), opt.Nquiv)));
zqind = unique(round(linspace(opt.bcind*length(z), (1-opt.bcind)*length(z), opt.Nquiv)));
xqind(xqind==0) = [];
zqind(zqind==0) = [];

% make grid for quiver
[Xquiv,Zquiv] = meshgrid(x(xqind),z(zqind));


v = varargin;
for vi = 1:length(varargin)
    v{vi} = varargin{vi}(:,zqind,xqind,:,:);
end



% rescale u, w?
if opt.uquiv
    D = max(z(end)-z(1), x(end)-x(1));
    
    for vi = 1:length(varargin)
        v{vi} = v{vi}./scl*D/length(zqind);
    end
end

varargout = v;

end