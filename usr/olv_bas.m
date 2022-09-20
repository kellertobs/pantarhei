% pantarhei run script 
% calibrated for olivine + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'test';         % identifier for this run
outdir = '../out/';         % directory to save output files
nop    = 1;               % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = 1000;               % number of grid points in each direction
D      = 50;                % domain dimension in each direction [delta0]
BC     = 'closed';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = nop*500;           % maximum number of time steps
tend   = 1e16;              % model run time [s]

nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.7;               % second-order iterative step size (reduce if not converging)
cfl    = 0.9;               % Courant number to limit physical time step size
flim   = 1e-4;              % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.90; 0.10];       % initial background phase fractions (unity sum!)
dfg  = [-0.00;+0.00];       % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00;+0.00];       % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;        % impose gaussian mass transfer rate (unity sum!)

rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.6945, 0.1832; 0.5360, 0.1834;];  % permission slopes
B = [ 0.6906, 0.3094; 0.9993, 0.0007;];  % permission step locations
C = [ 0.6889, 0.1750; 0.8154, 1.5642;];  % permission step widths

% check problem scales (returns delta0, w0)
scales;

% reset domain depth to multiple of max segr-comp-length
D  = D.*max(delta0(:));
h  = D/N;

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));

% additional functions to calculate initial f
% fInit = @(X,Z) PermBarrier(X, Z, BC, f0, 0.004*[1;-1], [0.25,0.75], 0.1);
% fInit = @(X,Z) Gaussian2D(X, Z, D, f0, [0.01,0.02].*[-1;1], 1/32);
% fInit = @(X,Z) Gaussian1D(X, Z, D, f0, dfg, 0.05);
% fInit = @(X,Z) InitFromMatFile (X, Z, f0, 0.4*[-1;1], 2, 'fsoltest2');
% fInit = @(X,Z) RichardSolWave (X, Z, f0, 2, delta0);

% run model
run('../src/pantarhei');
