% params for MMS for periodic boundary conditions
% YQW, 8 Dec 2020

mms = true;

% manufactured solution for [f, p, u, w];
% equation form is Amf + dmf*cos(t/Tmf)*cos(x/Xmf)*cos(z/Zmf)
Tmf = 1/(2*pi)*1e5*ones(NPHS,4);
Xmf = 0.25/(2*pi)*D ./ones(NPHS,4);
Zmf = 0.25/(2*pi)*D ./ones(NPHS,4);

Amf = [f0, zeros(NPHS,3)];

% sine wave amplitude
dp  = 1e+1*sign(dfr);
du  = 1e-4*sign(dfr);
dw  = 1e-4*sign(dfr);
dmf = [dfr, dp, du, dw];

% dmf = [dfr, [1e2,1e-3,1e-3].*ones(NPHS,3)];

