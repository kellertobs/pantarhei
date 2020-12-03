% pantarhei run script 
% calibrated for olivine + plagioclase + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'olv30_plg10_bas60';
nop    = 5;                 % plot and store output every [nop] time step
svop   = 1;                 % save output and print figures

NPHS   = 3;                 % number of phases
N      = 100;               % number of grid points in each direction
D      = 0.25;              % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
tend   = 1e16;              % model run time [s]
dt     = 1e0;               % initial time step [s]

nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 1e-6;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 3000;              % maximum iteration count for iterative solver
alpha  = 0.90;              % first-order iterative step size (reduce if not converging)
beta   = 0.60;              % second-order iterative step size (reduce if not converging)
cfl    = 1.00;              % Courant number to limit physical time step size
flim   = 1e-9;              % limit phase fractions in coefficient closures
thtlim = 1e+9;              % limit phase-internal permission contrasts
cfflim = 1e+9;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.30; 0.10; 0.60]; % initial background phase fractions (unity sum!)
dfg  = [ 0.10;-0.10; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00; 0.00; 0.00]; % initial random perturbation amplitude (unity sum!)
smth = (D/30)^2/h^2;        % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [ 3300; 2400; 2800]; % pure-phase densities
eta0 = [1e+18;1e+16;1e+01]; % pure-phase viscosities
d0   = [ 5e-3; 5e-3; 5e-3]; % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A  =  [ 0.60, 0.30, 0.25; ...
        0.60, 0.30, 0.25; ...
        0.30, 0.25, 0.20; ];  % permission slopes
    
B  =  [ 0.45, 0.35, 0.20; ...
        0.35, 0.30, 0.35; ...
        0.45, 0.54, 0.01; ];  % permission step locations
    
C  =  [ 0.40, 0.40, 0.20; ...
        0.40, 0.40, 0.20; ...
        0.60, 0.20, 0.60; ];  % permission step widths

% run model
run('../src/pantarhei');
