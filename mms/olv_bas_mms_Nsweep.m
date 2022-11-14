% pantarhei run script
% calibrated for olivine + basaltic melt
% here to verify code using MMS

% clear workspace
clear; close all; clc;

% run method of manufactured solutions for time step 0 to check consistent
% f, p, u, w

% Nvec = [400,240,200,150,120,100,80,50,40,25];
Nvec = 100;
beta = 0.70;

for ni = 1:length(Nvec)
    beta = RunModel(Nvec(ni), beta);
end


%%
function [beta] = RunModel (Nin, beta)

mms = true;

% set model parameters
RunID  = ['olvbas_mms_' num2str(Nin,'%03d')];
% RunID = 'tmp';
outdir = '../out/';         % directory to save output files
nop    = 1;                 % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = Nin;          % number of grid points in each direction
Lfac   = 40;                % domain dimension in each direction [delta0]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 0;                 % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'fromm', 'weno5', 'tvdim'
nupd   = 50;                % update residual and permissions every [nupd] iterations
atol   = 1e-8;              % residual tolerance for convergence of iterative solver
rtol   = 1e-10;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 5000;              % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
% beta   = 0.80;              % second-order iterative step size (reduce if not converging)
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

flag = 0;
% run model; use while loop in case initial guess for beta is too large
while flag==0 && beta>=0
    fprintf(1,'N = %d, beta = %.1f. \n', Nin, beta);
    try
        run('../src/pantarhei');
        
        % assign flags to know why pantarhei finished running
        if res<atol, flag = 1;              % abs tol reached
        elseif it>=maxits, flag = 0;        % max iterations reached
        end
        
        NormErr = [fNormErr, pNormErr, uNormErr, wNormErr];
        MaxmErr = [fMaxErr , pMaxErr , uMaxErr , wMaxErr ];
        
        name = [outdir,RunID,'/',RunID,'_mmserr.mat'];
        save(name, 'NormErr', 'MaxmErr');
        
    catch exception
        % print the error
        getReport(exception, 'basic')
        flag = 0;
    end
    beta = beta-0.1;
end
end
