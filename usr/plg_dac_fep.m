% pantarhei run script 
% calibrated for plagioclase + dacitic melt + Fe-P-rich melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'plg50_dac40_fep10';
nop    = 1;                  % plot and store output every [nop] time step
svop   = 0;                 % save output and print figures

NPHS   = 3;                 % number of phases
N      = 100;               % number of grid points in each direction
D      = 100;               % domain dimension in each direction [delta0]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 1e3;               % maximum number of time steps
tend   = 1e16;              % model run time [s]

nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 1e-6;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.50;              % second-order iterative step size (reduce if not converging)
cfl    = 1.00;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [0.50; 0.40; 0.10];  % initial background phase fractions (unity sum!)
dfg  = [ 0.00; 0.00; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.01;-0.01; 0.02]; % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)

rho0 = [3000 ;2500;4000];   % pure-phase densities
eta0 = [1e+16;1e+3;1e+0];   % pure-phase viscosities
d0   = [1e-2 ;1e-2;1e-2];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% check problem scales (returns delta0, w0)
scales;

% reset domain depth to multiple of max segr-comp-length
D  = D.*max(delta0(:));
h  = D/N;

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));

% run model
run('../src/pantarhei');
