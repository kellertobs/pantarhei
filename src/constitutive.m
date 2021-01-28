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
if strcmp(BC,'closed'); qvxz(:,[1,end],:) = 0; end

% get volume flux fields
qfx = - (Kf(:,:,im)+Kf(:,:,ip))./2 .* (diff(p(:,:,ic),1,3)./h - Gx_pstar) ...
    + ( f(:,:,im)+ f(:,:,ip))./2 .* u;
qfz = - (Kf(:,im,:)+Kf(:,ip,:))./2 .* (diff(p(:,ic,:),1,2)./h - Gz_pstar) ...
    + ( f(:,im,:)+ f(:,ip,:))./2 .* w;
if strcmp(BC,'closed'); qfx(:,:,[1,end]) = 0; qfz(:,[1,end],:) = 0; end

% get parameterised mass transfer fields
Gm  = Gmg.*gsn;

% get momentum transfer fields
Gvx = - (Cv(:,:,im)+Cv(:,:,ip))./2 .* (u-ustar) + pstar_Gfx;
Gvz = - (Cv(:,im,:)+Cv(:,ip,:))./2 .* (w-wstar) + pstar_Gfz;

%get volume transfer fields
Gf  = -             Cf             .* (p-pstar) + vstar_Gf + Gm./rhostar;

% get momentum source fields
Qvx = (f(:,:,im)+f(:,:,ip))./2.*((rho(:,:,im)+rho(:,:,ip))./2-rhomix).*grav(2);
Qvz = (f(:,im,:)+f(:,ip,:))./2.*((rho(:,im,:)+rho(:,ip,:))./2-rhomix).*grav(1);