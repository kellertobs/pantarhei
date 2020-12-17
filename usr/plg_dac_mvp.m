% pantarhei run script 
% calibrated for plagioclase + dacitic melt + magmatic volatile phase

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'plg60_dac20_mvp20';
nop    = 0;                 % plot and store output every [nop] time step
svop   = 0;                 % save output and print figures

NPHS   = 3;                 % number of phases
N      = 100;               % number of grid points in each direction
D      = 20;                % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
tend   = 1e9;               % model run time [s]
dt     = 1e2;               % initial time step [s]

nupd   = 50;                % update residual and permissions every [nupd] iterations
npr    = 50;                % print iteration diagnostics every [npr] iterations
atol   = 1e-6;              % absolute residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % relative residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.60;              % second-order iterative step size (reduce if not converging)
cfl    = 1.00;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+4;              % limit phase-internal permission contrasts
cfflim = 1e+4;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [0.60; 0.20; 0.20];  % initial background phase fractions (unity sum!)
dfg  = [-0.04; 0.02; 0.02]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00; 0.00; 0.00]; % initial random perturbation amplitude (unity sum!)
smth = (D/30)^2/h^2;        % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [5e-3 ;5e-3;5e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% run model
run('../src/pantarhei');
