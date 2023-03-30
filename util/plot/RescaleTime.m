function [trscl, tunit, tscl] = RescaleTime (t)
% 
% rescale the units of time into something that makes sense for the viewer
% 
% YQW, 9 Nov 2022

tunitvec =         {'s'; 'hr'; 'day';  'yr' ; 'kyr'};
tscales  = cumprod([ 1 ; 3600;   24 ; 365.25;  1e3 ]);


[~,ti] = min(abs(log10(max(t)./tscales)));

tscl  = tscales(ti);
trscl = t./tscales(ti);
tunit = tunitvec{ti};


end