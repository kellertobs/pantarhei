% pantarhei run script
% calibrated for olivine + basaltic melt
% here to verify code using MMS

% clear workspace
clear; close all; clc;

% run method of manufactured solutions for time step 0 to check consistent
% f, p, u, w

Nvec = [120, 160, 200, 240, 280];
Nvecfrac = (Nvec(1)./Nvec).^2;
beta = 0.7;
NormErr = nan(2,4,length(Nvec)); MaxmErr = NormErr;

%%

for ni = 1:length(Nvec)

    [NormErr(:,:,ni), MaxmErr(:,:,ni),beta] = RunModel(Nvec(ni), Nvecfrac(ni), beta);
    beta = min(beta+0.1, 0.9);

    ne = squeeze(sum(NormErr,1));
    ne./ne(:,1)
    Nvec(1).^2./Nvec.^2

end

%%
function [NormErr,MaxmErr,beta] = RunModel (Nin, Nvecfrac, beta)

mms = true;
mms_regime = "olv_bas_suspension"

% set model parameters
RunID  = ['olvbas_mms_susp_' num2str(Nin,'%03d')];
% RunID = 'tmp';
outdir = '../out/';         % directory to save output files
pltits = false;              % whether to plot residuals 
nop    = 1;                 % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

NPHS   = 2;                 % number of phases
N      = Nin;               % number of grid points in each direction
Lfac   = 40;                % domain dimension in each direction [delta0]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
NtMax  = 0;                 % maximum number of time steps
tend   = 1e16;              % model run time [s]

advn   = 'weno5';           % advection scheme. best ones: 'quick', 'fromm', 'weno5', 'tvdim'
nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5*Nvecfrac;     % residual tolerance for convergence of iterative solver
rtol   = 1e-5;              % residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 10000;             % maximum iteration count for iterative solver
alpha  = 0.95;              % first-order iterative step size (reduce if not converging)
% beta   = 0.6;             % second-order iterative step size (reduce if not converging)
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
    close all;

    beta = beta-0.1;
end


disp('max div_qf error');


norm( sum(Div_qf-psrc,1), 'fro' ) / ( norm(Div_qf-psrc, 'fro') + TINY )
% figure; plot(Div_qf(:,:,21)-psrc(:,:,21), z, sum(Div_qf(:,:,21)-psrc(:,:,21)), z); axis tight; grid on;
% 
% dqffig = figure; 
% subplot(311); imagesc(x, z, squeeze(Div_qf(1,:,:)-psrc(1,:,:))); axis xy equal tight; colormap(ocean); colorbar;
% subplot(312); imagesc(x, z, squeeze(Div_qf(2,:,:)-psrc(2,:,:))); axis xy equal tight; colormap(ocean); colorbar;
% subplot(313); imagesc(x, z, squeeze( sum(Div_qf  -psrc   , 1))); axis xy equal tight; colormap(ocean); colorbar;
% 
% name = [outdir,RunID,'/',RunID,'_divqfpsrc_',num2str(step/nop)];
% print(dqffig,'-dpdf','-r200','-opengl',name,'-loose');

end
