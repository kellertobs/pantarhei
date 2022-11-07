function [vc] = GetCellCenterVals (v, dim)
% [vc] = GetCellCenterVals (v, dim)
% 
% use this function to interpolate the velocity fields, which are defined
% at the cell faces, to get the velocities at cell centers
% 
% INPUTS
% v     velocity field, defined at cell faces
% dim   which dimension (2 or 3) to interpolate over
% 
% OUTPUTS
% vc    velocity field, defined at cell centers
% 
% YQW, 9 Feb 2021


if dim==2
    vc = 0.5*(v(:,1:end-1,:,:,:) + v(:,2:end,:,:,:));
    
elseif dim==3
    vc = 0.5*(v(:,:,1:end-1,:,:) + v(:,:,2:end,:,:));
    
else
    vc = [];
    disp('!!! Wrong dimension defined !!!');
end



end