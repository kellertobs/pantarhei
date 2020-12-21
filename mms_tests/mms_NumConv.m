% pantarhei run script to verify code using MMS with a numerical
% convergence test for different grid sizes.
% calibrated for plagioclase + dacitic melt + magmatic volatile phase
% YQW, 9 Dec 2020

clear all; close all

Nvec = [10,15,20,25,30,40,50,60,80,100];
Nn   = length(Nvec);

% initialize output matrices
NormErr = nan(4, Nn);
MaxErr  = nan(4, Nn);
betaOut = nan(1, Nn);
flag    = nan(1, Nn);
Nit     = nan(1, Nn);

% choose one fixed alpha
alphaIn = 0.90;

% loop over N to get error between true and num solution
for ni = 1:Nn
    beta = alphaIn - 0.1;
    
    % loop to adjust beta to achieve convergence
    while flag(1,ni)~=1 && beta>0.4
        
        [NormErrOut, MaxErrOut, Nit(ni), flag(ni)] = ...
            RunSolver(Nvec(ni), alphaIn, beta);
        
        NormErr(:,ni) = NormErrOut;
        MaxErr (:,ni) = MaxErrOut ;
        betaOut(1,ni) = beta;
        
        if flag(ni)==2, keyboard; end
        
        beta = beta - 0.1;
    end
end

%% save outputs

clear Nn ni beta NormErrOut MaxErrOut ans

% load parameters used in solution to save together with errors
plg_dac_mvp_params;
alpha = alphaIn;
run('../mms/mms_periodic_params');

FileNameDefault = '../out/plg60_dac20_mvp20_mms/NumConvTest';

% if file exists, adjust final number of file to avoid overwriting
NfInDir = length(dir([FileNameDefault '*.mat']))+1;
FileName = [FileNameDefault '_' num2str(NfInDir) '.mat'];

save(FileName);


%% function that runs the solver for given N, alpha

function [NormErr, MaxErr, it, flag] = RunSolver (Nin, alphaIn, betaIn)

% initialize output error vectors
NormErr = nan(4,1);
MaxErr  = nan(4,1);
flag    = nan;

% load parameters
plg_dac_mvp_params;

% reassign number of grid points and grid spacing
N = Nin;
h = D/N;

% reassign convergence params
alpha = alphaIn;
beta  = betaIn;

fprintf(1, '    N = %d, alpha = %.2f, beta = %.2f.', N, alpha, beta);

% manufactured solution
run('../mms/mms_periodic_params');

try
    % run model
    run('../src/pantarhei');
    run('../mms/mms_results.m');
        
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


