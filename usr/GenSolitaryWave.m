function [fc, z, c] = GenSolitaryWave (A0, fc0, z0, zlim, N)
% [fctmp, z] = GenSolitaryWave (c, z0, zlim, N)
% 
% This function generates a 1D solitary wave solution. Under the right 
% conditions (constant viscosity, permeability is of power law 2), the 
% wave travels at speed c with its shape unchanged.
% Solution from Simpson and Spiegelman (2011). Since they offer an implicit
% equation, we first calculate z in terms of fc, then interpolate to
% desired z vector.
% 
% INPUTS 
% A0        wave amplitude, relative to fc0
% fc0       reference porosity at z = inf [scalar, dimensionless]
% z0        location of wave peak [scalar, length]
% zlim      domain limits, either
%               1) [2x1 vector, length]
%               2) [Nx1 vector, length] - in this case N should not be defined
% N         number of points
% NB: c, z0, zlim should have the same length unit
% 
% OUTPUTS
% fc        wave [Nx1, dimensionless %]
% z         positions where wave is calculated [Nx1, length]
% c         wave speed [length/time]
% 
% YQW, 17 Feb 2021

z = []; fc = [];

% check format of zlim, N
if length(zlim)==2
    z = linspace(zlim(1), zlim(2), N);
elseif length(zlim)>2
    z = zlim;
else
    error('!!! Input zlim has the wrong size !!!');
end

% transform wave amplitude to reference f_infty = 1
A = A0/fc0;

% precalc for speed since constant
sA1 = sqrt(A-1);

% implicit equation relating position and phi_c, assuming z0 = 0
% eqn 8 from Simpson and Spiegelman (2011)
zfunc = @(phi_c) (A+0.5).*(-2*sqrt(A-phi_c) + 1./sA1.*log( (sA1-sqrt(A-phi_c)) ./ (sA1 + sqrt(A-phi_c)) ) ).^2; 

% calculate z in terms of phi_c
fctmp = [linspace(A,1.02,101),linspace(1.01,1,1e5)];
ztmp  = sqrt(zfunc(fctmp));

% reassign ztmp(fctmp == 1) to a finite value
ztmp(end) = interp1(fctmp(end-10:end-1),ztmp(end-10:end-1),fctmp(end),'pchip','extrap');

% remove all points greater than the max zlim desired
zmax  = max(abs(zlim-z0));
fctmp(ztmp>zmax) = []; ztmp(ztmp>zmax) = []; 

% now set fc(z = large value) == 1
zmaxmax  = 10*max(abs(zlim-z0));
fcmaxmax = 1;

% assemble left and right sides of ztmp, fctmp
% move peak to z0
z1 = [-zmaxmax,  -fliplr( ztmp(2:end)),  ztmp,  zmaxmax] + z0;
f1 = [ fcmaxmax,  fliplr(fctmp(2:end)), fctmp, fcmaxmax];

% interpolate to desired z vector
fc = fc0*interp1(z1,f1,z,'linear','extrap');
c  = (2*A + 1).*fc0;
end