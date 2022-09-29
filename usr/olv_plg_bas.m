% pantarhei run script 
% calibrated for olivine + plagioclase + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'olvplgbas';
outdir = '../out/';         % directory to save output files
nop    = 1;                 % plot and store output every [nop] time step
svop   = 0;                 % save output and print figures
restart= 0;

NPHS   = 3;                 % number of phases
N      = 200;               % number of grid points in each direction
D      = 100;               % domain dimension in each direction [delta0]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 1e3;               % maximum number of time steps
tend   = 1e16;              % model run time [s]

nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 1e-6;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.99;              % first-order iterative step size (reduce if not converging)
beta   = 0.85;              % second-order iterative step size (reduce if not converging)
cfl    = 1.00;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.20; 0.20; 0.60]; % initial background phase fractions (unity sum!)
dfg  = [ 0.00;-0.00; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.02; 0.02; 0.00]; % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [0;-0;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [ 3200; 2400; 2700]; % pure-phase densities
eta0 = [1e+18;1e+15;1e+01]; % pure-phase viscosities
d0   = [ 5e-3; 5e-3; 5e-3]; % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A  =  [ 0.60, 0.30, 0.25; ...
        0.40, 0.25, 0.20; ...
        0.30, 0.25, 0.20; ];  % permission slopes
    
B  =  [ 0.45, 0.35, 0.20; ...
        0.35, 0.30, 0.35; ...
        0.45, 0.54, 0.01; ];  % permission step locations
    
C  =  [ 0.40, 0.40, 0.20; ...
        0.40, 0.40, 0.20; ...
        0.60, 0.20, 0.60; ];  % permission step widths
    
% check problem scales (returns delta0, w0)
scales;

% reset domain depth to multiple of max segr-comp-length
D  = D.*max(delta0(:));
h  = D/N;

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));

% run model
run('../src/pantarhei');
