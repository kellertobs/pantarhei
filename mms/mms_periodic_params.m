% params for MMS
% YQW, 8 Dec 2020

mms = true;

% manufactured solution for [f, p, u, w];
% equation form is Amf + dmf*cos(t/Tmf)*cos(x/Xmf)*cos(z/Zmf)
Amf = [f0, zeros(NPHS,3)];
dmf = [dfr, [1e4,1e-6,1e-6].*ones(3,3)];
Tmf = 1/(2*pi)*1e5*ones(NPHS,4);
Xmf = 1/(2*pi)*D ./ones(NPHS,4);
Zmf = 1/(2*pi)*D ./ones(NPHS,4);
