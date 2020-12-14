% pantarhei run script 
% calibrated for olivine + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'olv40_bas60';
nop    = 1;              % plot and store output every [nop] time step
svop   = 0;              % save output

NPHS   = 2;              % number of phases
N      = 100;            % number of grid points in each direction
D      = 0.5;            % domain dimension in each direction [m]
h      = D/N;            % grid spacing [m]
BC     = 'closed';       % boundary conditions: 'open', 'closed', 'periodic'
tend   = 1e16;           % model run time [s]
dt     = 1e0;            % initial time step [s]

nupd   = 50;             % update residual and permissions every [nupd] iterations
npr    = 1000;           % print iteration diagnostics every [npr] iterations
atol   = 1e-6;           % residual tolerance for convergence of iterative solver
rtol   = 1e-4;           % residual tolerance for convergence of iterative solver
minits = 500;            % minimum iteration count for iterative solver
maxits = 4000;           % maximum iteration count for iterative solver
alpha  = 0.80;           % first-order iterative step size (reduce if not converging)
beta   = 0.40;           % second-order iterative step size (reduce if not converging)
cfl    = 1.00;           % Courant number to limit physical time step size
flim   = 1e-4;           % limit phase fractions in coefficient closures
thtlim = 1e+4;           % limit phase-internal permission contrasts
cfflim = 1e+4;           % limit inter-phase coefficient contrasts

grav = [-9.81,0];        % gravity in vertical and horizontal direction
f0   = [ 0.40; 0.60];    % initial background phase fractions (unity sum!)
dfg  = [-0.00; 0.00];    % initial guassian peak amplitude (unity sum!)
dfr  = [ 0.02;-0.02];    % initial random perturbation amplitude (unity sum!)
smth = (D/30)^2/h^2;     % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;     % impose gaussian mass transfer rate (unity sum!)

rho0 = [ 3000; 2500];    % pure-phase densities
eta0 = [1e+18;1e+02];    % pure-phase viscosities
d0   = [ 5e-3; 5e-3];    % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.6945, 0.1832; 0.5360, 0.1834;];  % permission slopes
B = [ 0.6906, 0.3094; 0.9993, 0.0007;];  % permission step locations
C = [ 0.6889, 0.1750; 0.8154, 1.5642;];  % permission step widths

% run model
run('../src/pantarhei');