
% get phase density array
rho = PHS.rho0.*ones(size(f));

% get permission weights
F  = permute(repmat(f,1,1,1,NUM.NPHS),[4,1,2,3]);
Sf = (F./PHS.B).^(1./PHS.C);  Sf = Sf./sum(Sf,2);
Xf = sum(PHS.A.*Sf,2).*F + (1-sum(PHS.A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(PHS.Mv.^Xf,2));  gmv = geomean(geomean(thtv,3),2);
thtf = squeeze(prod(PHS.Mf.^Xf,2));  gmf = geomean(geomean(thtf,3),2);
thtv = 1./(1./thtv + 1./(NUM.thtlim.^0.5*gmv)) + (gmv/NUM.thtlim.^0.5);
thtf = 1./(1./thtf + 1./(NUM.thtlim.^0.5*gmf)) + (gmf/NUM.thtlim.^0.5);

% get momentum and volume flux and transfer coefficients
Kv =    f .*PHS.eta0           .*thtv;
Kf =    f .*PHS.d0.^2./PHS.eta0.*thtf;
Cv = (1-f)./PHS.d0.^2.*Kv;
Cf = (1-f)./PHS.d0.^2.*Kf; 

Kv = Kv + max(geomean(geomean(Kv,3),2))./NUM.cfflim;
Kf = Kf + max(geomean(geomean(Kf,3),2))./NUM.cfflim;
Cv = 1./(1./Cv + 1./(min(geomean(geomean(Cv,3),2)).*NUM.cfflim));
Cf = 1./(1./Cf + 1./(min(geomean(geomean(Cf,3),2)).*NUM.cfflim));

% get cell face values to calculate dtau's
  fx = 1/2*(  f(:,:,NUM.imx)+  f(:,:,NUM.ipx));   fz = 1/2*(  f(:,NUM.imz,:)+  f(:,NUM.ipz,:));
rhox = 1/2*(rho(:,:,NUM.imx)+rho(:,:,NUM.ipx)); rhoz = 1/2*(rho(:,NUM.imz,:)+rho(:,NUM.ipz,:));
 Kvx = 1/2*( Kv(:,:,NUM.imx)+ Kv(:,:,NUM.ipx));  Kvz = 1/2*( Kv(:,NUM.imz,:)+ Kv(:,NUM.ipz,:));
 Kfx = 1/2*( Kf(:,:,NUM.imx)+ Kf(:,:,NUM.ipx));  Kfz = 1/2*( Kf(:,NUM.imz,:)+ Kf(:,NUM.ipz,:));
 Cvx = 1/2*( Cv(:,:,NUM.imx)+ Cv(:,:,NUM.ipx));  Cvz = 1/2*( Cv(:,NUM.imz,:)+ Cv(:,NUM.ipz,:));
 Cfx = 1/2*( Cf(:,:,NUM.imx)+ Cf(:,:,NUM.ipx));  Cfz = 1/2*( Cf(:,NUM.imz,:)+ Cf(:,NUM.ipz,:));

% cell corner
Kvn = 1/4*(Kv(:,NUM.imz,NUM.imx)+Kv(:,NUM.ipz,NUM.imx)+Kv(:,NUM.imz,NUM.ipx)+Kv(:,NUM.ipz,NUM.ipx)); 

% get coefficient-based weights
omvc = Cv./sum(Cv,1);  omvx = Cvx./sum(Cvx,1);  omvz = Cvz./sum(Cvz,1);
omfc = Cf./sum(Cf,1);  omfx = Cfx./sum(Cfx,1);  omfz = Cfz./sum(Cfz,1);
omdp = Kf./sum(Kf,1);  ompx = Kfx./sum(Kfx,1);  ompz = Kfz./sum(Kfz,1);

% update reference and auxiliary fields
ustar = sum(omvx.*u  ,1);  Du = u-ustar;
wstar = sum(omvz.*w  ,1);  Dw = w-wstar;
pstar = sum(omfc.*p  ,1);  Dp = p-pstar;

% advection on reference velocity
vstar_Gf  = advect(f, ustar, wstar, NUM.h, {NUM.advn, 'vdf'}, [2,3], NUM.BC);

% get iterative pseudo-time steps
dtau_u = 3./4./( (1-fx+NUM.dmp).*Kvx./(NUM.h/NUM.ndim)^2 + (1-omvx).*Cvx + (1-fx).*(fx.^2./Cfx)./(NUM.h/NUM.ndim)^2 ) ;
dtau_w = 3./4./( (1-fz+NUM.dmp).*Kvz./(NUM.h/NUM.ndim)^2 + (1-omvz).*Cvz + (1-fz).*(fz.^2./Cfz)./(NUM.h/NUM.ndim)^2 ) ;
dtau_p = 1./4./( (1-f         ).*Kf ./(NUM.h/NUM.ndim)^2 + (1-omfc).*Cf  + (1-f ).*(f .^2./Cv )./(NUM.h/NUM.ndim)^2 ) ;

