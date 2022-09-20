% pantarhei run script 
% calibrated for plagioclase + dacitic melt + magmatic volatile phase

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'pdm_conv100';
outdir = '../out/';         % directory to save output files
nop    = 5;                 % plot and store output every [nop] time step
svop   = 1;                 % save output and print figures
restart= 0;

NPHS   = 3;                 % number of phases
N      = 250;               % number of grid points in each direction
D      = [100,100];               % domain dimension in each direction [delta0]
BC     = 'closed';          % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 1e3;               % maximum number of time steps
tend   = 1e16;              % model run time [s]

nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % residual tolerance for convergence of iterative solver
rtol   = 1e-4;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.50;              % second-order iterative step size (reduce if not converging)
cfl    = 0.90;              % Courant number to limit physical time step size
flim   = 1e-6;              % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.10; 0.85; 0.05]; % initial background phase fractions (unity sum!)
dfg  = [-0.09; 0.09; 0.00]; % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00; 0.00; 0.00]; % initial random perturbation amplitude (unity sum!)
smth = (N/40)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1;0].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)
Pu   = 0;                   % pure shear strain rate [/s]
Si   = 0;                   % simple shear strain rate [/s]

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% check problem scales (returns delta0, w0)
scales;

% reset domain depth to multiple of max segr-comp-length
D  = D.*max(delta0(:));
h  = D(1)/N;

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));

fInit = @(X,Z) initcond(X, Z, f0, dfg, dfr, D, h, 0.5, 2, smth);

% run model
run('../src/pantarhei');


%% initial condition function

function [f] = initcond (X, Z, f0, dfg, dfr, D, h, lambda, dA, smth)

% initialize
[NPHS, Nz, Nx] = size(Z);

% make a sine wave 
x = permute(X(1,1,:),[1,3,2]);
zboundary = repmat( cos(2*pi.*x./(lambda*D(2)))*dA*h, Nz, 1);

inds     = squeeze(Z(1,:,:)) > (0.4*D(1) + zboundary);
df       = dfr(1)*randn(1,Nz,Nx);   % random noise
df(inds) = df(inds) + dfg(1);       % add the top/bottom phase fraction diff

% smooth out the random noise
rng(15);
icz = [1,1:Nz,Nz];
icx = [1,1:Nx,Nx];

for i = 1:smth
    df = df + diff(df(:,icz,:),2,2)./8 + diff(df(:,:,icx),2,3)./8;
end
df(df<-f0(1)) = -f0(1);
f = f0 + [df;-df;zeros(size(df))];

end