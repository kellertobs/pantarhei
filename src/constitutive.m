
% update reference and auxiliary fields
ustar     = sum(omvx.*u,1);
wstar     = sum(omvz.*w,1);
pstar     = sum(omfc.*p,1);
rhostar   = sum(omfc.*rho,1);
Gx_pstar  = sum(ompx.*diff(p(:,:,icx),1,3)./h,1);
Gz_pstar  = sum(ompz.*diff(p(:,icz,:),1,2)./h,1);
pstar_Gfx = (pstar(:,:,imx)+pstar(:,:,ipx))./2.*diff(f(:,:,icx),1,3)./h;
pstar_Gfz = (pstar(:,imz,:)+pstar(:,ipz,:))./2.*diff(f(:,icz,:),1,2)./h;
Div_v     = diff(u,1,3)./h + diff(w,1,2)./h;


% volume flux fields
qfx = - Kfx.*(diff(p(:,:,icx),1,3)./h - Gx_pstar) + fx.*u; 
qfz = - Kfz.*(diff(p(:,icz,:),1,2)./h - Gz_pstar) + fz.*w;
if strcmp(BC{1},'closed'); qfz(:,[1,end],:) = 0; end
if strcmp(BC{2},'closed'); qfx(:,:,[1,end]) = 0; end

% divergence of volume fluxes (combine diffusive & advective parts)
Div_qf = diff(qfx,1,3)./h + diff(qfz,1,2)./h;

% parameterised mass transfer fields
Gm  = Gmg.*ones(size(f));

%get volume transfer fields
Gf  = Cf .*(p-pstar) - vstar_Gf - Gm./rhostar;


% momentum flux fields
qvxx  = - Kv .* (diff(u,1,3)./h - Div_v./3 + dmp*sum(Div_qf) - Pu) + f.*p;
qvzz  = - Kv .* (diff(w,1,2)./h - Div_v./3 + dmp*sum(Div_qf) + Pu) + f.*p;
qvxz  = - Kvn.* (diff(u(:,icz,:),1,2)./h + diff(w(:,:,icx),1,3)./h + 2*Si)./2;

% divergence of momentum fluxes
Div_qvx = diff(qvxx(:,:,icx),1,3)./h + diff(qvxz,1,2)./h;
Div_qvz = diff(qvzz(:,icz,:),1,2)./h + diff(qvxz,1,3)./h;

% momentum transfer fields
Gvx = Cvx.*(u-ustar) - pstar_Gfx;
Gvz = Cvz.*(w-wstar) - pstar_Gfz;

% momentum source fields
Qvx = -fx.*(rhox-rhomix).*grav(2);
Qvz = -fz.*(rhoz-rhomix).*grav(1);

