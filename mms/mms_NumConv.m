% pantarhei run script to verify code using MMS with a numerical
% convergence test for different grid sizes.
% calibrated for olivine rich rock and basaltic melt
% YQW, 9 Dec 2020

clear all; close all

mpScript = 'olv_bas_params';

f0   = [ 0.10; 0.90];   % initial background phase fractions (unity sum!)
dfr  = [ 0.01;-0.01];   % initial random perturbation amplitude (unity sum!)

Nvec = [25,40,50,60,80,100,120,150,200,400];
Nn   = length(Nvec);
Din  = 40;

% initialize output matrices
NormErr = nan(4, Nn);
MaxErr  = nan(4, Nn);
betaOut = nan(1, Nn);
flag    = nan(1, Nn);
Nit     = nan(1, Nn);

% loop over N to get error between true and num solution
for ni = 1:Nn
    beta = 0.80;
    
    % loop to adjust beta to achieve convergence
    while flag(1,ni)~=1 && beta>0
        
        [NormErrOut, MaxErrOut, Nit(ni), flag(ni)] = ...
            RunSolver(mpScript, Din, Nvec(ni), f0, dfr, beta);
        
        NormErr(:,ni) = NormErrOut;
        MaxErr (:,ni) = MaxErrOut ;
        betaOut(1,ni) = beta;
                
        beta = beta - 0.1;
    end
end


%% function that runs the solver for given N, alpha

function [NormErr, MaxErr, it, flag] = RunSolver (mpScript, Din, Nin, f0, dfr, betaIn)

% initialize output error vectors
NormErr = nan(4,1);
MaxErr  = nan(4,1);
flag    = nan;

% load phase parameters
run(mpScript);
RunID = [strjoin(strcat(PHS(:),num2str(round(f0*100), '%02d'))', '_') '_mms_N' num2str(Nin)];


% load solver params and reassign based on inputs
solver_params;
N    = Nin;
beta = betaIn;

% check problem scales (returns delta0, w0)
run('../usr/scales.m');

% reset domain depth to multiple of max segr-comp-length
D    = Din.*max(delta0(:));
h    = D/N;
smth = (N/40)^2;            % smoothing parameter for random perturbation field

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));

% get properties for manufactured solution
switch BC
    case 'periodic'
        run('../mms/mms_utils/mms_params_periodicBC.m');
        
    case 'closed'
        % this is not tested yet
        run('../mms/mms_utils/mms_params_closedBC.m');
end

fprintf(1, '    N = %d, beta = %.2f.', N, beta);

try
    % run model
    run('../src/pantarhei');
        
    % output 2norm and maximum errors
    NormErr = [fNormErr; pNormErr; uNormErr; wNormErr];
    MaxErr  = [fMaxErr ; pMaxErr ; uMaxErr ; wMaxErr ];
    
    % assign flags to know why pantarhei finished running
    if res<atol, flag = 1;              % abs tol reached
    elseif it>=maxits, flag = 2;        % max iterations reached
    end
    
catch exception
    
    % print the error
    getReport(exception, 'basic')
end
end


