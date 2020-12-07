% model parameters for use in manufactured solution plg_dac_mvp in ../usr/
% YQW, 25 Nov 2020


NPHS   = 3;                 % number of phases
N      = 40;               % number of grid points in each direction
D      = 80;               % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
tend   = 1e5;
dt     = 0.05*tend;                 % initial time step [s]

grav = [-9.81,0];           % gravity in vertical and horizontal direction
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e16;1e2;1e-3];   % pure-phase viscosities
d0   = [5e-3 ;5e-3;5e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths


%% parameters for numerical solution

nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.99;              % first-order iterative step size (reduce if not converging)
beta   = 0.90;              % second-order iterative step size (reduce if not converging)
cfl    = 0.75;              % Courant number to limit physical time step size
flim   = 1e-16;             % limit phase fractions in coefficient closures
thtlim = 1e+16;             % limit phase-internal permission contrasts
cfflim = 1e+16;              % limit inter-phase coefficient contrasts
