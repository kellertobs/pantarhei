% run the code for various times to check if residual - source function
% is below error
% YQW, 1 Dec 2020

clear all; restoredefaultpath
addpath('./SymsFuncs/');

MfSol_Params;
MfSol_BCPeriodic;

time = 1e4;

ttmp = 2*pi*linspace(0, Tf(1), 201);
% figure; plot(ttmp, cos(ttmp./Tf(1)), time, cos(time./Tf(1)), 'rx'); axis tight

%%

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

%% initialize fields

% define x and z positions of field variables
XuGrid = cat(3,X-h/2,X(:,:,end)+h/2);   ZuGrid = cat(3,Z,Z(:,:,end));
XwGrid = cat(2,X,X(:,end,:));           ZwGrid = cat(2,Z-h/2,Z(:,end,:)+h/2);

% generate manufactured solutions
u        = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tu, Xu, Zu, U0, dU0);
w        = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tw, Xw, Zw, W0, dW0);
p        = MMSsource('Calc_p', time, X     , Z     , Tp, Xp, Zp, P0, dP0);
[f,dfdt] = MMSsource('Calc_f', time, X     , Z     , Tf, Xf, Zf, f0, df0);

Gm  = Gmg.*ones(size(p));
rho = rho0.*ones(size(p));
rhomix = mean(mean(sum(f.*rho,1)));

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

%% check segregation compaction length is on order of domain

dsc = SegCompLength(f, eta0, d0, A, B, C);

%%

% initialise coefficient and auxiliary fields
run('../src/up2date');

%{
[min(dtau_u(:)), max(dtau_u(:))]
[min(dtau_w(:)), max(dtau_w(:))]
[min(dtau_p(:)), max(dtau_p(:))]
[min(dtau_f(:)), max(dtau_f(:))]
%}

%% calculate derivatives for Kv, Kf

dKvdx = diff(Kv(:,:,ic),1,3)./h;
dKvdz = diff(Kv(:,ic,:),1,2)./h;

dKfdx = diff(Kf(:,:,ic),1,3)./h;
dKfdz = diff(Kf(:,ic,:),1,2)./h;

domKfdx = diff(omfx(:,:,ic),1,3)./h;

%% get parameters from analytical solution

[f_x,dfdt_x,dfdx_x,dfdz_x] = MMSsource('Calc_f', time,XuGrid,ZuGrid,Tf,Xf,Zf,f0,df0);

[Kvs, Kfs, Cvs, Cfs] = MMSsource('CalcPermissions', f_x, eta0, d0, A, B, C, thtlim, cfflim);
[dKvs,dKfs] = MMSsource('GetPermissionDerivs', f_x, dfdt_x, dfdx_x, dfdz_x, d0, eta0, A, B, C);
[dKfs2] = MMSsource('CutPermissionDerivs', Kfs, dKfs, cfflim);
[dKvs2] = MMSsource('CutPermissionDerivs', Kvs, dKvs, cfflim);
[~  , omKf, omCv, omCf, domKf] = MMSsource('CalcWeights', Kvs, Kfs, Cvs, Cfs, dKfs2);

%%
% plot_field(dKfdx); plot_field(dKfs(:,:,:,2)); plot_field(dKfs2(:,:,:,2)); 
% plot_field(dKvdx); plot_field(dKvs(:,:,:,2)); plot_field(dKvs2(:,:,:,2)); 
% plot_field((dKvdx - dKvs2(:,:,:,2))./(dKvdx));
plot_field(domKfdx); plot_field(domKf(:,:,:,2)); 
plot_field((domKfdx-domKf(:,:,:,2)));


%%
% update reference and auxiliary fields
ustar     = sum(omvx.*u,1);
wstar     = sum(omvz.*w,1);
pstar     = sum(omfc.*p,1);
rhostar   = sum(omfc.*rho,1);
Gx_pstar  = sum(omfx.*diff(p(:,:,ic),1,3)./h,1);
Gz_pstar  = sum(omfz.*diff(p(:,ic,:),1,2)./h,1);
Div_v     = diff(u,1,3)./h + diff(w,1,2)./h;
pstar_Gfx = (pstar(:,:,im)+pstar(:,:,ip))./2.*diff(f(:,:,ic),1,3)./h;
pstar_Gfz = (pstar(:,im,:)+pstar(:,ip,:))./2.*diff(f(:,ic,:),1,2)./h;

% get momentum flux fields
qvxx  = - Kv .* (diff(u,1,3)./h - Div_v./3) + f.*p;
qvzz  = - Kv .* (diff(w,1,2)./h - Div_v./3) + f.*p;
qvxz  = -(Kv(:,im,im)+Kv(:,ip,im)+Kv(:,im,ip)+Kv(:,ip,ip))./4 ...
    .* (diff(u(:,ic,:),1,2)./h + diff(w(:,:,ic),1,3)./h)./2;

% get volume flux fields
qfx = - (Kf(:,:,im)+Kf(:,:,ip))./2 .* (diff(p(:,:,ic),1,3)./h - Gx_pstar) ...
    + ( f(:,:,im)+ f(:,:,ip))./2 .* u;
qfz = - (Kf(:,im,:)+Kf(:,ip,:))./2 .* (diff(p(:,ic,:),1,2)./h - Gz_pstar) ...
    + ( f(:,im,:)+ f(:,ip,:))./2 .* w;

% get momentum transfer fields
Gvx = - (Cv(:,:,im)+Cv(:,:,ip))./2 .* (u-ustar) + pstar_Gfx;
Gvz = - (Cv(:,im,:)+Cv(:,ip,:))./2 .* (w-wstar) + pstar_Gfz;

%get volume transfer field
Gf  = -             Cf             .* (p-pstar) + vstar_Gf + Gm./rhostar;

% get momentum source fields
Qvx = (f(:,:,im)+f(:,:,ip))./2.*((rho(:,:,im)+rho(:,:,ip))./2-rhomix).*grav(2);
Qvz = (f(:,im,:)+f(:,ip,:))./2.*((rho(:,im,:)+rho(:,ip,:))./2-rhomix).*grav(1);

% get physical time step
dt  = cfl/(max(abs([qfx(:);qfz(:)]))/(h/2) + max(abs(Gf(:)))./5e-3);  % [s]

% get residual fields
res_u =             + diff(qvxx(:,:,ic),1,3)./h + diff(qvxz,1,2)./h - Gvx - Qvx    ;
res_w =             + diff(qvzz(:,ic,:),1,2)./h + diff(qvxz,1,3)./h - Gvz - Qvz    ;
res_p =             + diff(qfx         ,1,3)./h + diff(qfz ,1,2)./h - Gf  - Gm./rho;
res_f =  dfdt                                                       + Gf           ;

%% symbolic source functions

[usrc]       =  usrcFunc(time);
[wsrc]       =  wsrcFunc(time);
[psrc, fsrc] = pfsrcFunc(time);

%%

udiff = res_u - usrc;
wdiff = res_w - wsrc;
pdiff = res_p - psrc;
fdiff = res_f - fsrc;

ur = norm(udiff(:).*dtau_u(:),2)./(norm(u(:),2)+1e-32);
wr = norm(wdiff(:).*dtau_w(:),2)./(norm(w(:),2)+1e-32);
pr = norm(pdiff(:).*dtau_p(:),2)./(norm(p(:),2)+1e-32);
fr = norm(fdiff(:).*dtau_f(:),2)./(norm(f(:),2)+1e-32);

res = ur + wr + pr + fr;

fprintf(1, '\n\n');
fprintf(1, 'mean dtau_u = %4.4e.\n', mean(dtau_u(:)));
fprintf(1, 'mean dtau_w = %4.4e.\n', mean(dtau_w(:)));
fprintf(1, 'mean dtau_p = %4.4e.\n', mean(dtau_p(:)));
fprintf(1, 'mean dtau_f = %4.4e.\n', mean(dtau_f(:)));
fprintf(1, '\n\n');

fprintf(1, 'mean res_u - usrc = %4.4e.\n', mean(udiff(:)));
fprintf(1, 'mean res_w - wsrc = %4.4e.\n', mean(wdiff(:)));
fprintf(1, 'mean res_p - psrc = %4.4e.\n', mean(pdiff(:)));
fprintf(1, 'mean res_f - fsrc = %4.4e.\n', mean(fdiff(:)));
fprintf(1, '\n\n');

fprintf(1, '(res_u x dtau_u)/u = %4.4e.\n', ur);
fprintf(1, '(res_w x dtau_w)/w = %4.4e.\n', wr);
fprintf(1, '(res_p x dtau_p)/p = %4.4e.\n', pr);
fprintf(1, '(res_f x dtau_f)/f = %4.4e.\n', fr);
fprintf(1, 'total residual     = %4.4e.\n', res);
fprintf(1, '\n\n');

% plot_field(udiff, ['udiff, N=' num2str(N)]);
% plot_field(wdiff, ['wdiff, N=' num2str(N)]);
% plot_field(pdiff, ['pdiff, N=' num2str(N)]);
% plot_field(fdiff, ['fdiff, N=' num2str(N)]);






