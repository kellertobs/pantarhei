function [x, z, zunit] = RescaleLengths (x, z)
% 
% rescale the units of length into something that makes sense for the viewer
% 
% YQW, 9 Nov 2022

zunitvec = { 'mm'; 'cm'; 'm'; 'km'};
zscales  = [ 1e-3; 1e-2;  1 ;  1e3];


[~,zi] = min(abs(log10(max(z)./zscales)));

x     = x./zscales(zi);
z     = z./zscales(zi);
zunit = zunitvec{zi};


end