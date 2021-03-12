% params for MMS with closed boundary conditions
% YQW, 8 Dec 2020

mms = true;

% manufactured solution for [f, p, u, w];
% equation form is Amf + dmf*cos(t/Tmf)*cos(x/Xmf)*cos(z/Zmf)
Tmf = 1/(2*pi)*1e5*ones(NPHS,4);
Xmf = 1/(2*pi)*D ./ones(NPHS,4);
Zmf = 1/(2*pi)*D ./ones(NPHS,4);

Amf = [f0, zeros(NPHS,3)];

dp  = 1e+1*sign(dfr);
du  = 1e-4*sign(dfr);
dw  = 1e-4*sign(dfr);
dmf = [dfr, dp, du, dw];

