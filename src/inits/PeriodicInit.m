function [f] = PeriodicInit (X, Z, D, f0, df, Nw)
% 
% [f] = PeriodicInit (X, Z, D, f0, df, Nw)
% 
% this function makes an initial condition with periodic phase fraction 
% (BC should also be periodic)
% 
% INPUTS
% X, Z      meshgrid values of grid cells [NPHS x Nz x Nz]
% D         domain size [scalar]
% f0        background phase fraction, below z0 [NPHS x 1]
% df        phase fraction step change, above z0 [NPHS x 1]
% Nw        number of cycles per domain [1x1 or 2x1]
% 
% OUTPUTS
% f         phase fraction field [NPHS x Nz x Nz]

% wave number
k = 2*pi*Nw/D;

if length(k)==1
    % cycles only in z direction
    f = f0 + df.*sin(k*Z);
elseif length(k)==2
    % cycles in both x and z directions
    f = f0 + df.*sin(k*Z).*sin(k*X);
end

end