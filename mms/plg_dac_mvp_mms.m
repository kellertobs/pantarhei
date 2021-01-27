% pantarhei run script 
% calibrated for plagioclase + dacitic melt + magmatic volatile phase
% here to verify code using MMS
% YQW, 8 Dec 2020

% clear workspace
clear; close all; 

% set model parameters
RunID  = 'plg60_dac20_mvp20_mms';
nop    = 1;                 % plot and store output every [nop] time step
svop   = 0;                 % save output and print figures

NPHS   = 3;                 % number of phases
N      = 25;               % number of grid points in each direction
D      = 100;                % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
tend   = 0.99;              % model run time [s]
dt     = 1;                 % initial time step [s]
NtMax  = 1;                 % max number of time steps

nupd   = 100;               % update residual and permissions every [nupd] iterations
npr    = 100;               % print iteration diagnostics every [npr] iterations
atol   = 1e-6;              % absolute residual tolerance for convergence of iterative solver
rtol   = 1e-10;             % relative residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 1e5;               % maximum iteration count for iterative solver
alpha  = 0.99;              % first-order iterative step size (reduce if not converging)
beta   = 0.50;              % second-order iterative step size (reduce if not converging)
cfl    = 0.9;               % Courant number to limit physical time step size
flim   = 1e-16;             % limit phase fractions in coefficient closures
thtlim = 1e+16;             % limit phase-internal permission contrasts
cfflim = 1e+5;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.60; 0.20; 0.20]; % initial background phase fractions (unity sum!)
dfg  = [-0.00; 0.00; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.04; 0.02; 0.02]; % initial random perturbation amplitude (unity sum!)
smth = (D/30)^2/h^2;        % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [5e-3 ;5e-3;5e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

%% manufactured solution

run('mms_utils/mms_params_periodicBC.m');

%%

% run model
run('../src/pantarhei');