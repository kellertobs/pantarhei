function [f] = PermBarrier (X, Z, BC, f0, df, zfrac, smopt)
% 
% [f] = PermBarrier (X, Z, BC, f0, df, zfrac, smopt)
% 
% this function makes an initial condition with a step change in phase
% fraction
% 
% INPUTS
% X, Z      meshgrid values of grid cells [NPHS x Nz x Nz]
% BC        type of boundary condition [closed or periodic]
% f0        background phase fraction, below z0 [NPHS x 1]
% df        phase fraction step change, above z0 [NPHS x 1]
% zfrac     step change location, as fraction of domain size [scalar or 2x1]
% smopt     whether to smooth and over how many grid cells 
%           [0 if no smoothing, scalar for number of grid cells before and
%           after step change]
% 
% OUTPUTS
% f         phase fraction field [NPHS x Nz x Nz]

% initialize
N = size(Z,2);

% check what kind of BCs
switch BC
    case 'closed'
        % just one step
        zInd = round(zfrac(1)*N):N;
    case 'periodic'
        % has to go up and back down to maintain periodic bc
        zInd = round(zfrac(1)*N):round(zfrac(2)*N);
end


% check smoothing
if smopt == 0   
    % no smoothing
    f           = f0.*ones(1,N,N);
    f(:,zInd,:) = f(:,zInd,:) + df;
else
    
    % smoothing over smopt cells
    Nmat = permute(repmat(1:N,N,1),[3,2,1]);
    f = f0 + df./(1+exp(-smopt*(Nmat - zInd(1))));
    
    if strcmp(BC, 'periodic')
        f = f  - df./(1+exp(-smopt*(Nmat - zInd(end))));
    end
end


end