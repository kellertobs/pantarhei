% pantarhei run script 
% calibrated for plagioclase + dacitic melt + Fe-P-rich melt

% clear workspace
clear; close all; %clc;

% set model parameters
% detect system
IO.RunID  = 'P5direct';
IO.outdir = '../out/';         % directory to save output files
% IO.outdir = '/media/43TB_RAID_Array/tchai/out/';
IO.pltits  = false;            % whether to plot residuals 
IO.figvis  = 'on';             % toggle for figure visibility on/off
IO.nop     = 50;               % plot and store output every [nop] time step
IO.svop    = 1;                % save output
IO.restart = 0;

NUM.NPHS   = 3;                 % number of phases
NUM.N      = 5000;              % number of grid points in each direction
NUM.Lfac   = 8;                 % domain dimension in each direction [delta0]
NUM.BC     = 'opentop';         % boundary conditions: 'open', 'closed', 'periodic'
NUM.NtMax  = 1e4;               % maximum number of time steps
NUM.tend   = 1e16;              % model run time [s]

NUM.direct = true;              % select direct solver (false: multigrid solver)
NUM.tint   = 'bd2si';           % time integration scheme. 'be1im', 'bd2im', 'cn2si', 'bd2si'
NUM.advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
NUM.smth   = 10;                % update residual and permissions every [nupd] iterations
NUM.atol   = 1e-9;              % residual tolerance for convergence of iterative solver
NUM.rtol   = 1e-6;              % residual tolerance for convergence of iterative solver
NUM.minits = 1;                 % minimum iteration count for iterative solver
NUM.maxits = 100;               % maximum iteration count for iterative solver
NUM.cfl    = 0.25;              % Courant number to limit physical time step size
NUM.flim   = 1e-6;              % limit phase fractions in coefficient closures
NUM.thtlim = 1e+6;              % limit phase-internal permission contrasts
NUM.cfflim = 1e+9;              % limit inter-phase coefficient contrasts

PHS.grav = [-9.81,0];           % gravity in vertical and horizontal direction
PHS.f0   = [ 0.57; 0.38; 0.05]; % initial background phase fractions (unity sum!)
PHS.dfg  = [ 0.00; 0.00; 0.00]; % initial guassian peak amplitude (unity sum!)
PHS.dfr  = [ 0.00; 0.00; 0.00]; % initial random perturbation amplitude (unity sum!)
PHS.smth = (NUM.N/40)^2;        % smoothing parameter for random perturbation field
PHS.Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)
PHS.Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
PHS.Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

PHS.rho0 = [2700 ; 2400; 4000]; % pure-phase densities
PHS.eta0 = [1e+16; 1e+5;    1]; % pure-phase viscosities
PHS.d0   = [1e-3 ; 1e-3; 1e-3]; % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
PHS.A  =  [ 0.25, 0.25, 0.25; ...
            0.25, 0.25, 0.25; ...
            0.25, 0.25, 0.25; ];  % permission slopes
PHS.B  =  [ 0.44, 0.18, 0.38; ...
            0.60, 0.02, 0.38; ...
            0.72, 0.25, 0.03; ];  % permission step locations
PHS.C  =  [ 0.30, 0.30, 0.30; ...
            0.60, 0.60, 0.12; ...
            0.60, 0.12, 0.60; ];  % permission step widths


% run model
run('../src/pantarhei');
