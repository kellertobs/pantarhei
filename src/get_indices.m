% update z-index arrays   
NUM.izv = 1:NUM.Nz+1;
if strcmp(NUM.BC{1},'periodic')
    NUM.icz = [NUM.Nz,1:NUM.Nz,1]; NUM.imz = [NUM.Nz,1:NUM.Nz]; NUM.ipz = [1:NUM.Nz,1];
elseif strcmp(NUM.BC{1}(1:4),'open') || strcmp(NUM.BC{1},'closed')
    NUM.icz = [1,1:NUM.Nz,NUM.Nz]; NUM.imz = [1,1:NUM.Nz]; NUM.ipz = [1:NUM.Nz,NUM.Nz];
end

% update x-index arrays    
NUM.ixv = 1:NUM.Nx+1;
if strcmp(NUM.BC{2},'periodic')
    NUM.icx = [NUM.Nx,1:NUM.Nx,1]; NUM.imx = [NUM.Nx,1:NUM.Nx]; NUM.ipx = [1:NUM.Nx,1];
elseif strcmp(NUM.BC{2}(1:4),'open') || strcmp(NUM.BC{2},'closed')
    NUM.icx = [1,1:NUM.Nx,NUM.Nx]; NUM.imx = [1,1:NUM.Nx]; NUM.ipx = [1:NUM.Nx,NUM.Nx];
end