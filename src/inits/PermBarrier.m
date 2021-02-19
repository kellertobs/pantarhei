function [f] = PermBarrier (BC, f0, df, zfrac, N)
% [f] = PermBarrier (f0, df, z0, zMat)
% 
% this function makes an initial condition with a step change in phase
% fraction
% 
% INPUTS
% f0        background phase fraction, below z0 [NPHS x 1]
% df        phase fraction step change, above z0 [NPHS x 1]
% zfrac     step change location, as fraction of domain size [scalar or 2x1]
% N         number of grid points [scalar]
% 
% OUTPUTS
% f         phase fraction field [NPHS x Nz x Nz]

switch BC
    case 'closed'
        % just one step
        zInd = round(zfrac(1)*N):N;
    case 'periodic'
        % has to go up and back down to maintain periodic bc
        zInd = round(zfrac(1)*N):round(zfrac(2)*N);
end

f           = f0.*ones(1,N,N);
f(:,zInd,:) = f(:,zInd,:) + df;



end