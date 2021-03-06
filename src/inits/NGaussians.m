function [f] = NGaussians (X, Z, D, f0, df, gw)
% 
% [f] = NGaussians (X, Z, D, f0, df, zfrac, smopt)
% 
% this function makes an initial condition multiple Gaussians, with the
% number of Gaussians defined by df
% All the Gaussians are centered at z=0 at regular x intervals
% 
% INPUTS
% X, Z      meshgrid values of grid cells [NPHS x Nz x Nz]
% D         domain size [scalar]
% f0        background phase fraction, below z0 [NPHS x 1]
% df        phase fraction step change, above z0 [NPHS x Ng]
% gw        Gaussian width, fraction of domain size [scalar]
% 
% OUTPUTS
% f         phase fraction field [NPHS x Nz x Nz]
% 
% YQW, 24 Feb 2021


Ng = size(df,2);                    % number of Gaussian peaks
x0 = X(1) + (1:2:(2*Ng))*0.5*D/Ng;  % centroids of the Gaussians
f  = f0.*ones(size(X));             % initialize phase fractions

for gi = 1:Ng
    gsn = exp(-(X-x0(gi)).^2./(gw*D).^2).*exp(-Z.^2./(gw*D).^2);
    f = f + df(:,gi).*gsn;
end




end