
clear variables;

% set model parameters
RunID  = 'twophase_melt';
nop    = 1;              % plot and store output every [nop] time step
svop   = 0;                 % save output

NPHS   = 2;              % number of phases
N      = 100;            % number of grid points in each direction
D      = 10e3;           % domain dimension in each direction [m]
h      = D/N;            % grid spacing [m]
BC     = 'open';         % boundary conditions periodic, open, or closed
tend   = 1e20;           % model run time [s]
dt     = 1e9;            % initial time step [s]

nupd   = 50;             % update residual and permissions every [nupd] iterations
atol   = 1e-5;           % residual tolerance for convergence of iterative solver
rtol   = 1e-3;           % residual tolerance for convergence of iterative solver
minits = 500;            % minimum iteration count for iterative solver
maxits = 2.5e3;          % maximum iteration count for iterative solver
alpha  = 0.95;           % first-order iterative step size (reduce if not converging)
beta   = 0.50;           % second-order iterative step size (reduce if not converging)
cfl    = 0.50;           % Courant number to limit physical time step size
flim   = 1e-3;           % limit to phase fractions in coefficient closures
thtlim = 1e+4;           % limit to phase-internal permission contrasts
cfflim = 1e+9;           % limit to coefficient contrasts between phases

grav = [-9.81,0];        % gravity in vertical and horizontal direction
f0   = [ 0.95; 0.05];    % initial background phase fractions (unity sum!)
dfg  = [-0.00; 0.00];    % initial guassian peak amplitude (unity sum!)
dfr  = [ 0.00;-0.00];    % initial random perturbation amplitude (unity sum!)
Gmg  = [1;-1].*2e-9;     % impose gaussian mass transfer rate (unity sum!)
smth = (D/20)^2/h^2;     % smoothing parameter for random perturbation field

rho0 = [3000;2500];      % pure-phase densities
eta0 = [1e+18;1e+2];     % pure-phase viscosities
d0   = [1e-3;1e-3];      % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.6945, 0.1832; 0.5360, 0.1834;];  % permission slopes
B = [ 0.6906, 0.3094; 0.9993, 0.0007;];  % permission step locations
C = [ 0.6889, 0.1750; 0.8154, 1.5642;];  % permission step widths

