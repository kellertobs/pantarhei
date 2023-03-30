% pantarhei run script 
% calibrated for plagioclase + dacitic melt + magmatic volatile phase

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'plgdacmvp';
outdir = '../out/';         % directory to save output files
pltits = false;             % whether to plot residuals 
nop    = 5;                 % plot and store output every [nop] time step
svop   = 1;                 % save output and print figures
restart= 0;

NPHS   = 3;                 % number of phases
N      = 250;               % number of grid points in each direction
Lfac   = [100,100];         % domain dimension in each direction [delta0]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 1e3;               % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.50;              % second-order iterative step size (reduce if not converging)
dmp    = 0;                 % damping parameter, acts as numerical bulk viscosity
cfl    = 0.90;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.10; 0.85; 0.05]; % initial background phase fractions (unity sum!)
dfg  = [-0.01; 0.01; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00; 0.00; 0.00]; % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)
Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [   0.60, 0.25, 0.30; 
        0.20, 0.20, 0.20; 
        0.20, 0.20, 0.20; ];  % permission slopes
B = [   0.30, 0.15, 0.55; 
        0.48, 0.02, 0.50;   
        0.80, 0.08, 0.12; ];  % permission step locations
C = [   0.20, 0.20, 0.20; 
        0.60, 0.60, 0.12; 
        0.20, 0.25, 0.50; ];  % permission step widths

% run model
run('../src/pantarhei');
