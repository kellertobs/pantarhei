
% define x and z positions of field variables
XuGrid = cat(3,X-h/2,X(:,:,end)+h/2);   ZuGrid = cat(3,Z,Z(:,:,end));
XwGrid = cat(2,X,X(:,end,:));           ZwGrid = cat(2,Z-h/2,Z(:,end,:)+h/2);

% generate manufactured solutions
u = MMSsource('Calc_u', 0, XuGrid, ZuGrid, Tu, Xu, Zu, U0, dU0);
w = MMSsource('Calc_w', 0, XwGrid, ZwGrid, Tw, Xw, Zw, W0, dW0);
p = MMSsource('Calc_p', 0, X     , Z     , Tp, Xp, Zp, P0, dP0);
f = MMSsource('Calc_f', 0, X     , Z     , Tf, Xf, Zf, f0, df0);

ut0 = u; wt0 = w; pt0 = p; ft0 = f;

Gm  = Gmg.*ones(size(f));
rho = rho0.*ones(size(f));

% assign auxiliary fields
ui = u;  ustar = mean(u,1);  usegr = 0*u;  res_u = 0*u;  dtau_u = res_u;
wi = w;  wstar = mean(w,1);  wsegr = 0*w;  res_w = 0*w;  dtau_w = res_w;
pi = p;  pstar = mean(p,1);  pcmpt = 0*p;  res_p = 0*p;  dtau_p = res_p;
fo = f;  fi = f;                           res_f = 0*f;  dtau_f = res_f;

% check segregation compaction length is on order of domain
dsc = SegCompLength(f, eta0, d0, A, B, C);

% assign source functions for brevity
usrcFunc  = @(t) MMSsource('Calc_usrc', t,XuGrid(1,:,:),ZuGrid(1,:,:),...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm,...
    Tp,Xp,Zp,P0,dP0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0,dU0, Tw,Xw,Zw,W0,dW0,...
    thtlim, cfflim);

wsrcFunc  = @(t) MMSsource('Calc_wsrc', t,XwGrid(1,:,:),ZwGrid(1,:,:),...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm,...
    Tp,Xp,Zp,P0,dP0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0,dU0, Tw,Xw,Zw,W0,dW0,...
    thtlim, cfflim);

pfsrcFunc = @(t) MMSsource('Calc_pfsrc', t,X(1,:,:),Z(1,:,:),...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm,...
    Tp,Xp,Zp,P0,dP0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0,dU0, Tw,Xw,Zw,W0,dW0,...
    thtlim, cfflim);


% check source at t=0
if (1), time = 0; step = 0; mms_check_diff; end

% usrc_func  = @(t) Run_Calc_src('u',t,XuGrid(1,:,:),ZuGrid(1,:,:),...
%     d0,eta0,rho,grav,A,B,C,Gm,...
%     Tp,Xp,Zp,P0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0, Tw,Xw,Zw,W0);
% 
% wsrc_func  = @(t) Run_Calc_src('w',t,XwGrid(1,:,:),ZwGrid(1,:,:),...
%     d0,eta0,rho,grav,A,B,C,Gm,...
%     Tp,Xp,Zp,P0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0, Tw,Xw,Zw,W0);
% 
% pfsrc_func = @(t) Run_Calc_src('pf',t,X(1,:,:),Z(1,:,:),...
%     d0,eta0,rho,grav,A,B,C,Gm,...
%     Tp,Xp,Zp,P0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0, Tw,Xw,Zw,W0);