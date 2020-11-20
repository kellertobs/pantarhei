
clear variables;

% set model parameters
RunID  = 'olvplgbas';
NPHS   = 3;                 % number of phases
N      = 100;               % number of grid points in each direction
D      = 10;                % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
BC     = 'open';            % boundary conditions open (periodic flow) or closed (no flow)
tend   = 1e24;              % model run time [s]

nop    = 1;                 % plot and store output every [nop] time step
nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 5e-4;              % residual tolerance for convergence of iterative solver
rtol   = 1e-2;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 2e3;               % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.00;              % second-order iterative step size (reduce if not converging)
cfl    = 0.125;             % Courant number to limit physical time step size
flim   = 1e-3;              % limit to phase fractions in coefficient closures
thtlim = 1e+3;              % limit to phase-internal permission contrasts
cfflim = 1e12;              % limit to coefficient contrasts between phases

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [0.25; 0.25; 0.50];  % initial background phase fractions (unity sum!)
dfg  = [0.00; 0.00; 0.00];  % initial guassian peak amplitude (unity sum!)
dfr  = [0.05; 0.05; 0.05];  % initial random perturbation amplitude (unity sum!)
smth = (D/20)^2/h^2;        % smoothing parameter for random perturbation field

rho0 = [3500 ;2400; 2800];  % pure-phase densities
eta0 = [1e+18;1e+16;1e+1];  % pure-phase viscosities
d0   = [3e-3 ;1e-2 ;1e-3];  % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A  =  [ 0.60, 0.30, 0.25; 0.30, 0.60, 0.25; 0.20, 0.20, 0.20; ];  % permission slopes
B  =  [ 0.40, 0.45, 0.15; 0.39, 0.35, 0.26; 0.48, 0.51, 0.01; ];  % permission step locations
C  =  [ 0.40, 0.20, 0.40; 0.20, 0.25, 0.25; 0.60, 0.20, 0.60; ];  % permission step widths
