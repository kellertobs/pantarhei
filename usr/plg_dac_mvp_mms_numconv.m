% pantarhei run script to verify code using MMS with a numerical
% convergence test for different grid sizes.
% calibrated for plagioclase + dacitic melt + magmatic volatile phase
% YQW, 9 Dec 2020

clear all; close all

Nvec = [10,15,20,25,30,35,40,50,60];
avec = [0.6:0.1:0.7];
% test different alphas because iteration convergence is affected by N.

Nn   = length(Nvec);
Na   = length(avec);

NormErr = zeros(4, Nn, Na);
MaxErr  = zeros(4, Nn, Na);

% loop over N and alpha to get error between true and num solution
for ai = 1:Na
    for ni = 1:Nn
        
        [NormErrOut, MaxErrOut] = RunSolver(Nvec(ni), avec(ai));
        NormErr(:,ni,ai) = NormErrOut;
        MaxErr (:,ni,ai) = MaxErrOut ;
    end
end

%% save outputs

clear Nn Na ai ni NormErrOut MaxErrOut ans

% load parameters used in solution to save together with errors
plg_dac_mvp_params;
run('../mms/mms_periodic_params');

FileNameDefault = '../out/plg60_dac20_mvp20_mms/NumConvTest';

% if file exists, adjust final number of file to avoid overwriting
NfInDir = length(dir([FileNameDefault '*.mat']))+1;
FileName = [FileNameDefault '_' num2str(NfInDir) '.mat'];

save(FileName);


%% function that runs the solver for given N, alpha

function [NormErr, MaxErr] = RunSolver (Nin, alphaIn)

% initialize output error vectors
NormErr = nan(4,1);
MaxErr  = nan(4,1);

% load parameters
plg_dac_mvp_params;

% reassign number of grid points and grid spacing
N = Nin;
h = D/N;

% reassign convergence params
alpha = alphaIn;
beta  = alpha - 0.1;

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
    
catch exception
    
    % print the error
    getReport(exception, 'basic')
end
end


