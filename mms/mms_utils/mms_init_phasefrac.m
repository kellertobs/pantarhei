
%% parameters for the manufactured solution for [f, p, u, w]

% general equation form (may be phase shifted or sign flipped)
% Amf + dmf*cos(t/Tmf)*cos(x/Xmf)*cos(z/Zmf)
%
% f = f0 + df0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
% p = P0 + dP0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
% u = U0 + dU0.*cos(tfrac).*sin(xfrac).*sin(zfrac);
% w = W0 + dW0.*cos(tfrac).*cos(xfrac).*cos(zfrac);

Amf = [f0, [0;0], [0;0], 6.6e-4.*[-1;1]];           % baseline amplitude of variation
Tmf =   1*ones(NPHS,4)./(2*deg2rad(180))*1e5;       %  timescale of variation
Xmf = 0.5*ones(NPHS,4)./(2*deg2rad(180))*L(2);      % wavelength of variation in x direction
Zmf = 0.5*ones(NPHS,4)./(2*deg2rad(180))*L(1);      % wavelength of variation in z direction

% sine/cosine wave amplitude
dp  = [  0.5 ; 0.5 ];
du  = [  4e-5; 4e-5];
dw  = [  8e-5; 5e-5].*sign(dfg);
dmf = [dfg, dp, du, dw];

%% now initialize the phase fractions and additional position matrices

% define x and z positions of field variables
XuGrid = cat(3,X-h/2,X(:,:,end)+h/2);   ZuGrid = cat(3,Z,Z(:,:,end));
XwGrid = cat(2,X,X(:,end,:));           ZwGrid = cat(2,Z-h/2,Z(:,end,:)+h/2);

% generate manufactured solutions
f = MMSsource('Calc_f', 0, X, Z, Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
 