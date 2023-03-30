% pantarhei run script 
% calibrated for olivine + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'plgdac_blayer';   % identifier for this run
outdir = '../out/';         % directory to save output files
nop    = 10;                % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = 200;               % number of grid points in each direction
Lfac   = 40;                % domain dimension in each direction [delta0]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = nop*0;             % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-5;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 50000;             % maximum iteration count for iterative solver
alpha  = 0.90;              % first-order iterative step size (reduce if not converging)
beta   = 0.70;              % second-order iterative step size (reduce if not converging)
dmp    = 0;                 % damping parameter, acts as numerical bulk viscosity
cfl    = 0.50;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.05; 0.95];       % initial background phase fractions (unity sum!)
dfg  = [-0.00;+0.00];       % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00;+0.00];       % initial random perturbation amplitude (unity sum!)
smth = (N/20)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;        % impose gaussian mass transfer rate (unity sum!)
Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

rho0 = [ 2600; 2300];       % pure-phase densities
eta0 = [1e+16;1e+05];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.9516, 0.2426; 0.0199, 0.2005 ];  % permission slopes
B = [ 0.5978, 0.4022; 0.9933, 0.0067 ];  % permission step locations
C = [ 0.1370, 9.5783; 0.7336, 1.8839 ];  % permission step widths

wInit = @(w, w0) wbl(w, w0, f0);

% run model
run('../src/pantarhei');
run('../src/checkmasscons.m');
%%


function [w,wstar] = wbl(win, w0, f0)

[Nz,Nx] = size(win,[2,3]);

wstar = -sum(w0(:)).*ones(size(win));
w  = (w0./f0 + wstar).*ones(size(win));
w(:,[1,end],:) = 0;

end