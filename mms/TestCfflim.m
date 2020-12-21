% pantarhei run script to test mms with different cfflim
% using calibration for plagioclase + dacitic melt + magmatic volatile phase
% YQW, 15 Dec 2020

clear all; close all

Nvec    = [10,20,30,50];       
cffvec  = 10.^[2,3,5,8];

Nn = length(Nvec);
Nc = length(cffvec);

NormErr = zeros(4, Nn, Nc);
MaxErr  = zeros(4, Nn, Nc);
betaOut = nan  (1, Nn, Nc);
flag    = nan  (1, Nn, Nc);

% loop over N, alpha, beta to get error between true and num solution
for ni = 1:Nn
    for ci = 1:Nc
        beta = 0.80;
        
        while isnan(flag(1,ni,ci)) && beta>=0.4
        
            [NormErrOut, MaxErrOut, betaOut(1,ni,ci), flag(1,ni,ci)] = ...
                RunSolver(Nvec(ni), beta, cffvec(ci));
            
            NormErr(:,ni,ci) = NormErrOut;
            MaxErr (:,ni,ci) = MaxErrOut ;
            
            beta = beta - 0.1;
            
        end
    end
end

%% save outputs

clear Nn Nc ci ni beta NormErrOut MaxErrOut ans

% load parameters used in solution to save together with errors
plg_dac_mvp_params;
run('mms_utils/mms_periodic_params');

FileNameDefault = '../out/plg60_dac20_mvp20_mms/TestCfflim';

% if file exists, adjust final number of file to avoid overwriting
NfInDir = length(dir([FileNameDefault '*.mat']))+1;
FileName = [FileNameDefault '_' num2str(NfInDir) '.mat'];

save(FileName);


%% function that runs the solver for given N, cfflim

function [NormErr, MaxErr, beta, flag] = RunSolver (Nin, betaIn, cfflimIn)

% initialize output error vectors
NormErr = nan(4,1);
MaxErr  = nan(4,1);
flag    = nan;

% load parameters
plg_dac_mvp_params;

% reassign number of grid points and grid spacing
N = Nin;
h = D/N;

% reassign cfflim
cfflim = cfflimIn; 

% hardcode alpha, allow adjustable beta
alpha = 0.9;
beta = betaIn;

fprintf(1, '    N = %d, beta = %.2f, cfflim = %.1e.\n', N, beta, cfflim);

% manufactured solution
run('mms_utils/mms_periodic_params');

try
    % run model
    run('../src/pantarhei');
    run('mms_utils/mms_results.m');
    
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

