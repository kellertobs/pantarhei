% pantarhei run script 
% calibrated for olivine + basaltic melt

% clear workspace
clear; close all; %clc;

% set model parameters
IO.RunID   = 'olvbas';          % identifier for this run
IO.outdir  = '../out/';         % directory to save output files
IO.pltits  = true;              % whether to plot residuals 
IO.figvis = 'on';               % toggle for figure visibility on/off
IO.nop     = 1;                 % plot and store output every [nop] time step
IO.svop    = 1;                 % save output
IO.restart = 0;

NUM.NPHS   = 2;                 % number of phases
NUM.N      = 2^8;               % number of grid points in each direction [power of 2]
NUM.Lfac   = 20;                % domain dimension in each direction [delta0]
NUM.BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NUM.NtMax  = IO.nop*10;         % maximum number of time steps
NUM.tend   = 1e16;              % model run time [s]

NUM.advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
NUM.smth   = 30;                % number of smoothing iterations at each MG level
NUM.atol   = 1e-5;              % residual tolerance for convergence of iterative solver
NUM.rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
NUM.minits = 2;                 % minimum iteration count for iterative solver
NUM.maxits = 100;               % maximum iteration count for iterative solver
NUM.alpha  = 0.80;              % first-order iterative step size (reduce if not converging)
NUM.beta   = 0.00;              % second-order iterative step size (reduce if not converging)
NUM.dmp    = 1;                 % damping parameter, acts as numerical bulk viscosity
NUM.minlvl = 2;                 % minimum size MG level [2^minlvl]
NUM.minlvl_tol  = 1e-3;         % residual tolerance for at minimum MG level
NUM.minlvl_smth = 1e3;          % residual tolerance for at minimum MG level
NUM.cfl    = 0.25;              % Courant number to limit physical time step size
NUM.flim   = 1e-6;              % limit phase fractions in coefficient closures
NUM.thtlim = 1e+6;              % limit phase-internal permission contrasts
NUM.cfflim = 1e+9;              % limit inter-phase coefficient contrasts

PHS.grav = [-9.81,0];           % gravity in vertical and horizontal direction
PHS.f0   = [ 0.10; 0.90];       % initial background phase fractions (unity sum!)
PHS.dfg  = [-0.00;+0.00];       % initial guassian peak amplitude (unity sum!)
PHS.dfr  = [-0.00;+0.00];       % initial random perturbation amplitude (unity sum!)
PHS.smth = (NUM.N/10)^2;        % smoothing parameter for random perturbation field
PHS.Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
PHS.Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

PHS.rho0 = [ 3200; 2700];       % pure-phase densities
PHS.eta0 = [1e+18;1e+02];       % pure-phase viscosities
PHS.d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
PHS.A = [ 0.5989, 0.1772; 0.0397, 0.1182 ];  % permission slopes
PHS.B = [ 0.6870, 0.3130; 0.9998, 0.0002 ];  % permission step locations
PHS.C = [ 9.0105, 0.1592; 0.7249, 3.5524 ];  % permission step widths

% run model
run('../src/pantarhei');
