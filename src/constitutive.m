% update reference and auxiliary fields
rhor   = sum(omfc.*rho,1);

Gx_pr  = sum(omfx.*diff(p(:,:,icx),1,3)./h,1);
Gz_pr  = sum(omfz.*diff(p(:,icz,:),1,2)./h,1);
pr_Gfx = (pr(:,:,imx)+pr(:,:,ipx))./2.*diff(f(:,:,icx),1,3)./h;
pr_Gfz = (pr(:,imz,:)+pr(:,ipz,:))./2.*diff(f(:,icz,:),1,2)./h;

% get volume flux fields
qfx = - Kfx .* (diff(p(:,:,icx),1,3)./h - Gx_pr) + fx .* u; 
qfz = - Kfz .* (diff(p(:,icz,:),1,2)./h - Gz_pr) + fz .* w;
if strcmp(BC{1},'closed'); qfz(:,[1,end],:) = 0; end
if strcmp(BC{2},'closed'); qfx(:,:,[1,end]) = 0; end

% divergence of volume fluxes (combine diffusive & advective parts)
Div_qf = diff(qfx,1,3)./h + diff(qfz,1,2)./h;

% get momentum flux fields
Div_v = diff(u,1,3)./h + diff(w,1,2)./h;
qvxx  = - Kv  .* (diff(u,1,3)./h - Div_v./3 - Pu + dmp.*sum(Div_qf)) + f.*p;
qvzz  = - Kv  .* (diff(w,1,2)./h - Div_v./3 + Pu + dmp.*sum(Div_qf)) + f.*p;
qvxz  = - Kvc .* (diff(u(:,icz,:),1,2)./h + diff(w(:,:,icx),1,3)./h + 2*Si)./2;

% divergence of momentum fluxes
Div_qvx = diff(qvxx(:,:,icx),1,3)./h + diff(qvxz,1,2)./h;
Div_qvz = diff(qvzz(:,icz,:),1,2)./h + diff(qvxz,1,3)./h;

% get parameterised mass transfer fields
Gm  = Gmg.*ones(size(f));

% get momentum transfer fields
Gvx = Cvx .* (u-ur) - pr_Gfx;
Gvz = Cvz .* (w-wr) - pr_Gfz;

% get volume transfer fields
Gf  = Cf  .* (p-pr) - vr_Gf - Gm./rhor;

% get momentum source fields
Qvx = -fx.*(rhox-rhomix).*grav(2);
Qvz = -fz.*(rhoz-rhomix).*grav(1);

% update physical time step [s]
dt   = min([ 2*dto;
             cfl/(max(abs([qfx(:);qfz(:)]))/(h/2) + max(abs(Gf(:)))./1e-2);
             cfl*0.5*h./max(abs([ushr(:);wshr(:)]) + TINY) ]);