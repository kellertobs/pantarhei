% pantarhei run script 
% calibrated for olivine + basaltic melt
% 
% This script runs a simulation where concentration shock-rarefaction waves
% develop during particle settling, because particle settling speeds
% decrease with increasing particle fraction. 
% 
% This script demonstrates 3 things:
%   1. 1D mode of pantarhei is called by setting dfg = dfr = 0;
%   2. how to define an fInit function to be called within pantarhei to get
%      an arbitrary initial condition on phase fraction
%   3. performance of the advection scheme. Since shocks form, the
%      non-oscillatory schemes (weno3,weno5,tvdim) show improved performance
%      over the simpler schemes (central, quick, fromm).

% clear workspace
clear; close all; %clc;

% set model parameters
RunID  = 'olv05_bas95_d01_w10';
outdir = '../out/';         % directory to save output files
nop    = 10;                % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = 1000;              % number of grid points in each direction
D      = 500;               % domain dimension in each direction [delta0]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = nop*100;           % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'fromm', 'weno5', 'tvdim'
nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 8000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.4;               % second-order iterative step size (reduce if not converging)
cfl    = 0.90;              % Courant number to limit physical time step size
flim   = 1e-16;             % limit phase fractions in coefficient closures
thtlim = 1e+16;             % limit phase-internal permission contrasts
cfflim = 1e+16;             % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.05; 0.95];       % initial background phase fractions (unity sum!)
dfg  = [ 0.00; 0.00];       % initial guassian peak amplitude (unity sum!)
dfr  = [ 0.00; 0.00];       % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;        % impose gaussian mass transfer rate (unity sum!)
Pu   = 0;                   % pure shear strain rate [/s]
Si   = 0;                   % simple shear strain rate [/s]

rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 1e-3; 1e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.5989, 0.1772; 0.0397, 0.1182 ];  % permission slopes
B = [ 0.6870, 0.3130; 0.9998, 0.0002 ];  % permission step locations
C = [ 9.0105, 0.1592; 0.7249, 3.5524 ];  % permission step widths

% special function for initialization
fInit = @(X,Z,Ddim) Gaussian1D(X, Z, Ddim, f0, [+0.10;-0.10], 0.1);

% run model
run('../src/pantarhei');

%%

function [f] = Gaussian1D (X, Z, D, f0, df, gw)
% 
% [f] = Gaussian1D (X, Z, D, f0, df, zfrac, smopt)
% 
% this function makes an initial condition with a 1D Gaussian, varying only
% in z and constant in x
% 
% INPUTS
% X, Z      meshgrid values of grid cells [NPHS x Nz x Nz]
% D         domain size in dimensional units [scalar or 2x1, meters]
% f0        background phase fraction, below z0 [NPHS x 1]
% df        phase fraction step change, above z0 [NPHS x Ng]
% gw        Gaussian width, fraction of domain size [scalar]
% 
% OUTPUTS
% f         phase fraction field [NPHS x Nz x Nz]
% 
% YQW, 24 Feb 2021

f = f0 + df.*exp(-(Z./(gw*D(1))).^2);

end
