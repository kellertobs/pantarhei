function [f] = Gaussian1D (X, Z, D, f0, df, gw)
% 
% [f] = Gaussian1D (X, Z, D, f0, df, zfrac, smopt)
% 
% this function makes an initial condition with a 1D Gaussian, varying only
% in z and constant in x
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
z0 = Z(1) + (1:2:(2*Ng))*0.5*D/Ng;  % centroids of the Gaussians

if length(gw)==1, gw = gw.*ones(Ng,1); end

f  = f0.*ones(size(X));             % initialize phase fractions

for gi = 1:Ng
    gsn = exp(-(Z-z0(gi)).^2./(gw(gi)*D).^2);
    f = f + df(:,gi).*gsn;
end





end