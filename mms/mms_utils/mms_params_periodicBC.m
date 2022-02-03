% params for MMS for periodic boundary conditions
% YQW, 8 Dec 2020

mms = true;

piterm = deg2rad(180);

% manufactured solution for [f, p, u, w];
% equation form is Amf + dmf*cos(t/Tmf)*cos(x/Xmf)*cos(z/Zmf)
Amf = [f0, zeros(NPHS,3)];
Tmf =    1/(2*piterm)*1e5*ones(NPHS,4);
Xmf = 0.25/(2*piterm)*D ./ones(NPHS,4);
Zmf = 0.25/(2*piterm)*D ./ones(NPHS,4);

% sine wave amplitude
dp  = 1e+0*sign(dfr);
du  = 1e-5*sign(dfr);
dw  = 1e-5*sign(dfr);
dmf = [dfr, dp, du, dw];

% dmf = [dfr, [1e2,1e-3,1e-3].*ones(NPHS,3)];

