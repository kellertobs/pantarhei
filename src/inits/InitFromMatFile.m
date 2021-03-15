function [f] = InitFromMatFile (X, Z, f0, A, w, MatFile)
% 
% [f] = InitFromMatFile (X, Z, MatFile)
% 
% this function gets the form of the initial condition from a mat file.
% 
% INPUTS
% X, Z      meshgrid values of grid cells [NPHS x Nz x Nz]
% f0        background phase fraction [NPHS x 1]
% A         amplitude stretch of the wave [NPHS x 1]
% w         stretching of the wave in z direction [scalar]
% MatFile   name of mat file where initial condition form is stored
% 
% OUTPUTS
% f         phase fraction field [NPHS x Nz x Nz]
% 
% YQW, 24 Feb 2021

load(MatFile,'f','x');

fOut = interp1(x, f, Z*w, 'pchip', 'extrap');
f    = f0 + A.*fOut;

end