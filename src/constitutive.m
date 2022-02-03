% update reference and auxiliary fields
ustar     = sum(omvx.*u,1);
wstar     = sum(omvz.*w,1);
pstar     = sum(omfc.*p,1);
rhostar   = sum(omfc.*rho,1);
Gx_pstar  = sum(omfx.*diff(p(:,:,icx),1,3)./h,1);
Gz_pstar  = sum(omfz.*diff(p(:,icz,:),1,2)./h,1);
Div_v     = diff(u,1,3)./h + diff(w,1,2)./h;
pstar_Gfx = (pstar(:,:,imx)+pstar(:,:,ipx))./2.*diff(f(:,:,icx),1,3)./h;
pstar_Gfz = (pstar(:,imz,:)+pstar(:,ipz,:))./2.*diff(f(:,icz,:),1,2)./h;

% get momentum flux fields
qvxx  = - Kv .* (diff(u,1,3)./h - Div_v./3) + f.*p;
qvzz  = - Kv .* (diff(w,1,2)./h - Div_v./3) + f.*p;
qvxz  = -(Kv(:,imz,imx)+Kv(:,ipz,imx)+Kv(:,imz,ipx)+Kv(:,ipz,ipx))./4 ...
    .* (diff(u(:,icz,:),1,2)./h + diff(w(:,:,icx),1,3)./h)./2;
% if strcmp(BC,'closed'); qvxz(:,[1,end],:) = 0; end

% get volume flux fields
qfx = - (Kf(:,:,imx)+Kf(:,:,ipx))./2 .* (diff(p(:,:,icx),1,3)./h - Gx_pstar) ...
      + ( f(:,:,imx)+ f(:,:,ipx))./2 .* u;
qfz = - (Kf(:,imz,:)+Kf(:,ipz,:))./2 .* (diff(p(:,icz,:),1,2)./h - Gz_pstar) ...
      + ( f(:,imz,:)+ f(:,ipz,:))./2 .* w;
if strcmp(BC,'closed'); qfx(:,:,[1,end]) = 0; qfz(:,[1,end],:) = 0; end

% get parameterised mass transfer fields
Gm  = Gmg.*gsn;

% get momentum transfer fields
Gvx = - (Cv(:,:,imx)+Cv(:,:,ipx))./2 .* (u-ustar) + pstar_Gfx;
Gvz = - (Cv(:,imz,:)+Cv(:,ipz,:))./2 .* (w-wstar) + pstar_Gfz;

%get volume transfer fields
Gf  = -             Cf             .* (p-pstar) + vstar_Gf + Gm./rhostar;

% get momentum source fields
Qvx = (f(:,:,imx)+f(:,:,ipx))./2.*((rho(:,:,imx)+rho(:,:,ipx))./2-rhomix).*grav(2);
Qvz = (f(:,imz,:)+f(:,ipz,:))./2.*((rho(:,imz,:)+rho(:,ipz,:))./2-rhomix).*grav(1);