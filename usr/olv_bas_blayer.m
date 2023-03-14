% pantarhei run script 
  % calibrated for olivine + basaltic melt

% clear workspace
clear; close all; clc;

% set model parameters
RunID  = 'olvbas_blayer';          % identifier for this run
outdir = '../out/';         % directory to save output files
pltits = true;
nop    = 4;                % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = 320;               % number of grid points in each direction
Lfac   = 80;                % domain dimension in each direction [delta0]
BC     = 'closed';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = nop*0;           % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'weno5'
nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-6;              % residual tolerance for convergence of iterative solver
rtol   = 1e-5;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 50000;              % maximum iteration count for iterative solver
alpha  = 0.5;              % first-order iterative step size (reduce if not converging)
beta   = 0.00;              % second-order iterative step size (reduce if not converging)
dmp    = 0;
cfl    = 0.5;               % Courant number to limit physical time step size
flim   = 1e-6;             % limit phase fractions in coefficient closures
thtlim = 1e+6;              % limit phase-internal permission contrasts
cfflim = 1e+6;              % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
f0   = [ 0.95; 0.05];       % initial background phase fractions (unity sum!)
dfg  = [-0.00;+0.00];       % initial guassian peak amplitude (unity sum!)
dfr  = [-0.00;+0.00];       % initial random perturbation amplitude (unity sum!)
smth = (N/20)^2;            % smoothing parameter for random perturbation field
Gmg  = [1;-1].*0e-9;        % impose gaussian mass transfer rate (unity sum!)
Pu   = 0;                   %   pure shear strain rate [multiple of max segr speed]
Si   = 0;                   % simple shear strain rate [multiple of max segr speed]

rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.5989, 0.1772; 0.0397, 0.1182 ];  % permission slopes
B = [ 0.6870, 0.3130; 0.9998, 0.0002 ];  % permission step locations
C = [ 9.0105, 0.1592; 0.7249, 3.5524 ];  % permission step widths

wInit = @(w, w0) wbl(w, w0, f0);
% pInit = @(p, p0) pbl(p, p0, f0);

% run model
run('../src/pantarhei');

%%
kfx = - Kfx.*(diff(p(:,:,icx),1,3)./h - Gx_pstar);
kfz = - Kfz.*(diff(p(:,icz,:),1,2)./h - Gz_pstar);
div_kf = diff(kfx,1,3)./h + diff(kfz,1,2)./h;

disp('max div_qf error'); max(sum(Div_qf+TINY-TINY))./max(Div_qf(:)+TINY)
disp('max     Gf error'); max(sum(    Gf+TINY-TINY))./max(    Gf(:)+TINY)
disp('max div_kf error'); max(sum(div_kf+TINY-TINY))./max(div_kf(:)+TINY)

figure(21); 
subplot(131); plot(Div_qf, z, sum(Div_qf), z); grid on; axis tight;
subplot(132); plot(    Gf, z, sum(    Gf), z); grid on; axis tight;
subplot(133); plot(div_kf, z, sum(div_kf), z); grid on; axis tight;

%%

function [w,wstar] = wbl(win, w0, f0)

[Nz,Nx] = size(win,[2,3]);

wstar = -sum(w0(:)).*ones(size(win));
w  = (w0./f0 + wstar).*ones(size(win));
w(:,[1,end],:) = 0;

end

function p = pbl(p, p0, f0)

p(1,[1,end],:) = p0.*[1,-1];

end
