
% update reference fields
ustar     = sum(omvx.*u  ,1);
wstar     = sum(omvz.*w  ,1);
pstar     = sum(omfc.*p  ,1);


% volume flux fields
qfx = - Kfx.*(diff(p(:,:,NUM.icx),1,3)./NUM.h - Gx_pstar) + fx.*u; 
qfz = - Kfz.*(diff(p(:,NUM.icz,:),1,2)./NUM.h - Gz_pstar) + fz.*w;
if strcmp(NUM.BC{1},'closed'); qfz(:,[1,end],:) = 0; end
if strcmp(NUM.BC{2},'closed'); qfx(:,:,[1,end]) = 0; end
if strcmp(NUM.BC{1},'opentop'); qfz(:,end,:) = qfz(:,end-1,:); end
if strcmp(NUM.BC{1},'openbot'); qfz(:,1  ,:) = qfz(:,2    ,:); end

% divergence of volume fluxes
Div_qf = diff(qfx,1,3)./NUM.h + diff(qfz,1,2)./NUM.h;

% get volume transfer fields
Gf  = Cf.*(p-pstar) - vstar_Gf;


% momentum flux fields
Div_v = diff(u,1,3)./NUM.h + diff(w,1,2)./NUM.h;

qvxx  = - Kv .* (diff(u,1,3)./NUM.h - Div_v./3 + NUM.dmp*sum(Div_qf,1) - PHS.Pu) + f.*p;
qvzz  = - Kv .* (diff(w,1,2)./NUM.h - Div_v./3 + NUM.dmp*sum(Div_qf,1) + PHS.Pu) + f.*p;
qvxz  = - Kvn.* (diff(u(:,NUM.icz,:),1,2)./NUM.h + diff(w(:,:,NUM.icx),1,3)./NUM.h + 2*PHS.Si)./2;

% divergence of momentum fluxes
Div_qvx = diff(qvxx(:,:,NUM.icx),1,3)./NUM.h + diff(qvxz,1,2)./NUM.h;
Div_qvz = diff(qvzz(:,NUM.icz,:),1,2)./NUM.h + diff(qvxz,1,3)./NUM.h;

% momentum transfer fields
Gvx = Cvx.*(u-ustar) - pstar_Gfx;
Gvz = Cvz.*(w-wstar) - pstar_Gfz;




