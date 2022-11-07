function [vf] = GetCellFaceVals (v, dim, BC)
% [vc] = GetCellCenterVals (v, dim)
% 
% use this function to interpolate the phase fraction/pressure fields, 
% which are defined at the cell centers, to get values at cell faces where
% velocity is defined
% 
% INPUTS
% v     field (phase fraction, pressure, etc), defined at cell centers
% dim   which dimension (2 or 3) to interpolate over
% 
% OUTPUTS
% vf    field, defined at cell faces
% 
% YQW, 24 June 2021

N = size(v,2);

switch BC
    case 'periodic'
        im = [N,1:N]; ip = [1:N,1];
    case {'open','closed'}
        im = [1,1:N]; ip = [1:N,N];
end


if dim==2
    vf = (v(:,im,:,:) + v(:,ip,:,:))./2;
    
elseif dim==3
    vf = (v(:,:,im,:) + v(:,:,ip,:))./2;
    
else
    vf = [];
    disp('!!! Wrong dimension defined !!!');
end



end