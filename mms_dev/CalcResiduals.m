function [res_u, res_w, res_p, res_f] = CalcResiduals (time, Din, Nin, Params, MfSol)

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

Gm  = Gmg.*ones(size(p));
rho = rho0.*ones(size(p));
rhomix = mean(mean(sum(f.*rho,1)));

% initialise coefficient and auxiliary fields
run('../src/up2date');

dKvdx = diff(Kv(:,:,ic),1,3)./h;
dKvdz = diff(Kv(:,ic,:),1,2)./h;

dKfdx = diff(Kf(:,:,ic),1,3)./h;
dKfdz = diff(Kf(:,ic,:),1,2)./h;

domKfdx = diff(omfx(:,:,ic),1,3)./h;

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
res_f =  (dfdt)                                                     + Gf           ;

end