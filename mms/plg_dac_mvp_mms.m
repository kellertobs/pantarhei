% pantarhei run script 
% calibrated for plagioclase + dacitic melt + magmatic volatile phase
% here to verify code using MMS
% YQW, 8 Dec 2020

% clear workspace
clear; close all; 

% set model parameters
RunID  = 'plg30_dac60_mvp10_mms';
outdir = '../out/';         % directory to save output files
nop    = 1;               % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 3;                 % number of phases
N      = 100;               % number of grid points in each direction
D      = 10;                % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
tend   = 0.99;              % model run time [s]
NtMax  = 0;                 % max number of time steps

nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % absolute residual tolerance for convergence of iterative solver
rtol   = 1e-10;             % relative residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 1e5;               % maximum iteration count for iterative solver
alpha  = 0.90;              % first-order iterative step size (reduce if not converging)
beta   = 0.80;              % second-order iterative step size (reduce if not converging)
cfl    = 0.9;               % Courant number to limit physical time step size
flim   = 1e-16;             % limit phase fractions in coefficient closures
thtlim = 1e+16;             % limit phase-internal permission contrasts
cfflim = 1e+16;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.30; 0.60; 0.10]; % initial background phase fractions (unity sum!)
dfg  = [-0.00; 0.00; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [ 0.02;-0.04; 0.02]; % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [5e-3 ;5e-3;5e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% check problem scales (returns delta0, w0)
run('../usr/scales.m');

% reset domain depth to multiple of max segr-comp-length
D  = D.*max(delta0(:));
h  = D/N;

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));


%% manufactured solution

run('mms_utils/mms_params_periodicBC.m');

%%

% run model
run('../src/pantarhei');
