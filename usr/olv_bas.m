% pantarhei run script 
% calibrated for olivine + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'olvbas';          % identifier for this run
outdir = '../out/';         % directory to save output files
nop    = 1;                 % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = 200;               % number of grid points in each direction
Lfac   = 50;                % domain dimension in each direction [delta0]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = nop*500;           % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.4;               % second-order iterative step size (reduce if not converging)
cfl    = 0.5;               % Courant number to limit physical time step size
flim   = 1e-6;             % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.90; 0.10];       % initial background phase fractions (unity sum!)
dfg  = [-0.00;+0.00];       % initial guassian peak amplitude (unity sum!)
dfr  = [-0.01;+0.01];       % initial random perturbation amplitude (unity sum!)
smth = (N/20)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;        % impose gaussian mass transfer rate (unity sum!)
Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.5989, 0.1772; 0.0397, 0.1182 ];  % permission slopes
B = [ 0.6870, 0.3130; 0.9998, 0.0002 ];  % permission step locations
C = [ 9.0105, 0.1592; 0.7249, 3.5524 ];  % permission step widths

% run model
run('../src/pantarhei');
