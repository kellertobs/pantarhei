% pantarhei run script 
% calibrated for olivine + basaltic melt
% here to verify code using MMS

% clear workspace
clear; close all; clc;

% run method of manufactured solutions for time step 0 to check consistent
% f, p, u, w
mms = true;
mms_regime = "olv_bas_suspension"

% set model parameters
RunID  = 'olvbas_mms';
outdir = '../out/';         % directory to save output files
pltits = true;             % whether to plot residuals 
nop    = 1;                 % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;                 

NPHS   = 2;                 % number of phases
N      = 120;               % number of grid points in each direction
Lfac   = 40;                % domain dimension in each direction [delta0]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 0;                 % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'fromm', 'weno5', 'tvdim'
nupd   = 100;                % update residual and permissions every [nupd] iterations
atol   = 1e-6;              % residual tolerance for convergence of iterative solver
rtol   = 1e-5;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
beta   = 0.70;              % second-order iterative step size (reduce if not converging)
dmp    = 0;                 % damping parameter, acts as numerical bulk viscosity
cfl    = 0.50;              % Courant number to limit physical time step size
flim   = 1e-8;              % limit phase fractions in coefficient closures
thtlim = 1e+8;              % limit phase-internal permission contrasts
cfflim = 1e+8;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.05; 0.95];       % initial background phase fractions (unity sum!)
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


% special function for initialization
fInit = @(X,Z,Ldim) phasefracsine(X, Z, Ldim, f0, dfg);

% run model
run('../src/pantarhei'); 

%%

if ~exist('psrc','var'), psrc = zeros(size(p)); end

disp('max div_qf error');

norm( sum(Div_qf-psrc,1), 'fro' ) / ( norm(Div_qf-psrc, 'fro') + TINY )
figure; plot(Div_qf(:,:,21)-psrc(:,:,21), z, sum(Div_qf(:,:,21)-psrc(:,:,21)), z); axis tight; grid on;

dqffig = figure; 
subplot(311); imagesc(x, z, squeeze(Div_qf(1,:,:)-psrc(1,:,:))); axis xy equal tight; colormap(ocean); colorbar;
subplot(312); imagesc(x, z, squeeze(Div_qf(2,:,:)-psrc(2,:,:))); axis xy equal tight; colormap(ocean); colorbar;
subplot(313); imagesc(x, z, squeeze( sum(Div_qf  -psrc   , 1))); axis xy equal tight; colormap(ocean); colorbar;

%%
wmean = mean(w, [2,3])
max(w-wmean, [],[2,3])
max(u      , [],[2,3])
max(p      , [],[2,3])

%%

function [f] = phasefracsine (X, Z, L, f0, df)

NPHS = length(f0);

Xmf = ones(NPHS,1)./(2*pi)*L(2);      % wavelength of variation in x direction
Zmf = ones(NPHS,1)./(2*pi)*L(1);      % wavelength of variation in z direction

xfrac    = X./Xmf;
zfrac    = Z./Zmf;
f        = f0 + df.*cos(xfrac).*cos(zfrac);

end
