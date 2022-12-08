% pantarhei run script 
% calibrated for olivine + basaltic melt
% here to verify code using MMS

% clear workspace
clear; close all; clc;

% run method of manufactured solutions for time step 0 to check consistent
% f, p, u, w
mms = true;

% set model parameters
RunID  = 'olvbas_mms';
outdir = '../out/';         % directory to save output files
nop    = 1;                 % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;                 

NPHS   = 2;                 % number of phases
N      = 100;               % number of grid points in each direction
Lfac   = 40;                % domain dimension in each direction [delta0]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 1;                 % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'fromm', 'weno5', 'tvdim'
nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 1e-6;              % residual tolerance for convergence of iterative solver
rtol   = 1e-10;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.90;              % first-order iterative step size (reduce if not converging)
beta   = 0.50;              % second-order iterative step size (reduce if not converging)
cfl    = 0.50;              % Courant number to limit physical time step size
flim   = 1e-16;              % limit phase fractions in coefficient closures
thtlim = 1e+16;              % limit phase-internal permission contrasts
cfflim = 1e+16;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.10; 0.90];       % initial background phase fractions (unity sum!)
dfg  = [-0.01; 0.01];       % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00; 0.00];       % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;        % impose gaussian mass transfer rate (unity sum!)
Pu   = 0;                   % pure shear strain rate [/s]
Si   = 0;                   % simple shear strain rate [/s]

rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents 

% set permission weight parameters for coefficient closure model
A = [ 0.5989, 0.1772; 0.0397, 0.1182 ];  % permission slopes
B = [ 0.6870, 0.3130; 0.9998, 0.0002 ];  % permission step locations
C = [ 9.0105, 0.1592; 0.7249, 3.5524 ];  % permission step widths

% run model
run('../src/pantarhei');
