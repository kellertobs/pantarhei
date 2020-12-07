function [dsc] = SegCompLength (f, eta0, d0, A, B, C)
% calculates the segregation-compaction length between two phases.
% 
% order of dsc for three phases:
% [0            f2segf1comp     f3segf1comp;
%  f1segf2comp  0               f3segf2comp;
%  f1segf3comp  f2segf3comp     0]

thtlim  = inf;
flim    = 0;
cfflim  = inf;

fmean = mean(mean(f,3),2);
[~, ~, Cv, Cf] = MMSsource('CalcPermissions', ...
    fmean, eta0, d0, A, B, C, thtlim, cfflim);

dsc = fmean.'.*fmean./sqrt(Cv.'.*Cf);
dsc = dsc - diag(diag(dsc));


end