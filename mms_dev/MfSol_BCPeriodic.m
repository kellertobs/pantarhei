% generates parameters for manufactured solution
% assuming periodic boundary conditions
% 
% YQW, 18 Nov 2020

BC = 'periodic';

% manufactured solution properties ----------------------------------------
Tf  = 1/(2*pi)*1e5*[1;1;1];     % period of variation for phi [s]
f0  = [0.55;0.2;0.25];            % mean value for phase fraction Gaussian
df0 = [-0.04;0.2;0.02];          % scale of variation of phase fraction
Xf  = 1/(2*pi)*D./[1;1;1];      % length scale of horizontal variation for phase frac [m]
Zf  = 1/(2*pi)*D./[1;1;1];      % length scale of vertical variation for phase frac [m]

Tp  = 1/(2*pi)*1e5*[1;1;1];     % period of variation for P [s]
P0  = [0;0;0];                  % mean value of pressure [Pa]
dP0 = 1e5*[1;1;1];             % scale of variation of pressure [Pa]
Xp  = 1/(2*pi)*D./[1;1;1];      % length scale of horizontal variation for pressure [m]
Zp  = 1/(2*pi)*D./[1;1;1];      % vertical decay length for pressure [m]

Tu  = 1/(2*pi)*1e5*[1;1;1];     % period of variation for u [s]
U0  = 1e-4*[-0.1;-1;1];      % mean value of u [m/s]
dU0 = 0.5*abs(U0);              % scale of variation of u [m/s]
Xu  = 1/(2*pi)*D./[1;1;1];      % length scale of horizontal variation for u [m]
Zu  = 1/(2*pi)*D./[1;1;1];      % length scale of vertical variation for u [m]

Tw  = 1/(2*pi)*1e5*[1;1;1];     % period of variation for w [s]
W0  = 1e-4*[-0.1;-1;1];      % mean value of w [m/s]
dW0 = 0.5*abs(W0);              % scale of variation of w [m/s]
Xw  = 1/(2*pi)*D./[1;1;1];      % length scale of horizontal variation for w [m]
Zw  = 1/(2*pi)*D./[1;1;1];      % length scale of vertical variation for w [m]
