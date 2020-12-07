function [u,w,p,f] = CalcFields (time, Din, Nin, Params, MfSol)

run(Params);
D = Din; N = Nin; h = D/N;

run(MfSol);

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% initialise coordinate arrays
x     = -D/2+h/2 : h : D/2-h/2;
[X,~] = meshgrid(x,x);
X     = repmat(X,1,1,NPHS);
Z     = permute(X,[3,2,1]);
X     = permute(X,[3,1,2]);

% initialise indexing for boundary condition stencils (periodic BC)
ic = [N,1:N,1]; im = [N,1:N]; ip = [1:N,1];

% initialize fields

% define x and z positions of field variables
XuGrid = cat(3,X-h/2,X(:,:,end)+h/2);   ZuGrid = cat(3,Z,Z(:,:,end));
XwGrid = cat(2,X,X(:,end,:));           ZwGrid = cat(2,Z-h/2,Z(:,end,:)+h/2);

% generate manufactured solutions
u        = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tu, Xu, Zu, U0, dU0);
w        = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tw, Xw, Zw, W0, dW0);
p        = MMSsource('Calc_p', time, X     , Z     , Tp, Xp, Zp, P0, dP0);
[f,dfdt] = MMSsource('Calc_f', time, X     , Z     , Tf, Xf, Zf, f0, df0);

end