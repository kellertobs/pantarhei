% pantarhei run script 
% calibrated for plagioclase + dacitic melt + Fe-P-rich melt

% clear workspace
clear; close all; clc;

% set model parameters
% detect system
RunID  = 'NEW1P2';
outdir = '../out/';         % directory to save output files
% outdir = '/media/43TB_RAID_Array/tchai/out/';
nop    = 10;                % plot and store output every [nop] time step
svop   = 1;                 % save output and print figures
restart= 0;

NPHS   = 3;                 % number of phases
N      = 1000;               % number of grid points in each direction
Lfac   = 10;                % domain dimension in each direction [delta0]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 1e3;               % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 1e-8;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 2e4;               % maximum iteration count for iterative solver
alpha  = 0.80;              % first-order iterative step size (reduce if not converging)
beta   = 0.30;              % second-order iterative step size (reduce if not converging)
cfl    = 0.10;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+4;              % limit phase-internal permission contrasts
cfflim = 1e+9;              % limit inter-phase coefficient contrasts
dmp    = 0.5;               % damping parameter, acts as numerical bulk viscosity

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.65; 0.25; 0.10]; % initial background phase fractions (unity sum!)
dfg  = [ 0.00; 0.00; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [ 0.00; 0.00; 0.00]; % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)
Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

rho0 = [2700 ; 2400; 4000]; % pure-phase densities
eta0 = [1e+16; 1e+5;    1]; % pure-phase viscosities
d0   = [1e-3 ; 1e-3; 1e-3]; % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A  =  [ 0.25, 0.25, 0.25; ...
        0.25, 0.25, 0.25; ...
        0.25, 0.25, 0.25; ];  % permission slopes
B  =  [ 0.44, 0.18, 0.38; ...
        0.60, 0.02, 0.38; ...
        0.72, 0.25, 0.03; ];  % permission step locations
C  =  [ 0.30, 0.30, 0.30; ...
        0.60, 0.60, 0.12; ...
        0.60, 0.12, 0.60; ];  % permission step widths



% run model
run('../src/pantarhei');
