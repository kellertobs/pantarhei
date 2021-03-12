
PHS  = {'olv','plg','bas'};
NPHS = length(PHS);

rho0 = [ 3200; 2400; 2700]; % pure-phase densities
eta0 = [1e+18;1e+15;1e+01]; % pure-phase viscosities
d0   = [ 5e-3; 5e-3; 5e-3]; % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A  =  [ 0.60, 0.30, 0.25; ...
        0.40, 0.25, 0.20; ...
        0.30, 0.25, 0.20; ];  % permission slopes
    
B  =  [ 0.45, 0.35, 0.20; ...
        0.35, 0.30, 0.35; ...
        0.45, 0.54, 0.01; ];  % permission step locations
    
C  =  [ 0.40, 0.40, 0.20; ...
        0.40, 0.40, 0.20; ...
        0.60, 0.20, 0.60; ];  % permission step widths
    