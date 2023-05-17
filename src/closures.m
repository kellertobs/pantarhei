
% advection on reference velocity
vstar_Gf  = advect(f, ustar, wstar, h, {advn, 'vdf'}, [2,3], BC);

% get permission weights
F  = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
Sf = (F./B).^(1./C);  Sf = Sf./sum(Sf,2);
Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));  gmv = geomean(geomean(thtv,3),2);
thtf = squeeze(prod(Mf.^Xf,2));  gmf = geomean(geomean(thtf,3),2);
thtv = 1./(1./thtv + 1./(thtlim.^0.5*gmv)) + (gmv/thtlim.^0.5);
thtf = 1./(1./thtf + 1./(thtlim.^0.5*gmf)) + (gmf/thtlim.^0.5);

% get momentum and volume flux and transfer coefficients
Kv =    f .*eta0       .*thtv;
Kf =    f .*d0.^2./eta0.*thtf;
Cv = (1-f)./d0.^2.*Kv;
Cf = (1-f)./d0.^2.*Kf; 

Kv = Kv + max(geomean(geomean(Kv,3),2))./cfflim;
Kf = Kf + max(geomean(geomean(Kf,3),2))./cfflim;
Cv = 1./(1./Cv + 1./(min(geomean(geomean(Cv,3),2)).*cfflim));
Cf = 1./(1./Cf + 1./(min(geomean(geomean(Cf,3),2)).*cfflim));

% get cell face values to calculate dtau's
  fx = 1/2*(  f(:,:,imx)+  f(:,:,ipx));   fz = 1/2*(  f(:,imz,:)+  f(:,ipz,:));
rhox = 1/2*(rho(:,:,imx)+rho(:,:,ipx)); rhoz = 1/2*(rho(:,imz,:)+rho(:,ipz,:));
 Kvx = 1/2*( Kv(:,:,imx)+ Kv(:,:,ipx));  Kvz = 1/2*( Kv(:,imz,:)+ Kv(:,ipz,:));
 Kfx = 1/2*( Kf(:,:,imx)+ Kf(:,:,ipx));  Kfz = 1/2*( Kf(:,imz,:)+ Kf(:,ipz,:));
 Cvx = 1/2*( Cv(:,:,imx)+ Cv(:,:,ipx));  Cvz = 1/2*( Cv(:,imz,:)+ Cv(:,ipz,:));
 Cfx = 1/2*( Cf(:,:,imx)+ Cf(:,:,ipx));  Cfz = 1/2*( Cf(:,imz,:)+ Cf(:,ipz,:));

% cell corner
Kvn = 1/4*(Kv(:,imz,imx)+Kv(:,ipz,imx)+Kv(:,imz,ipx)+Kv(:,ipz,ipx)); 

% get coefficient-based weights
omvc = Cv./sum(Cv,1);  omvx = Cvx./sum(Cvx,1);  omvz = Cvz./sum(Cvz,1);
omfc = Cf./sum(Cf,1);  omfx = Cfx./sum(Cfx,1);  omfz = Cfz./sum(Cfz,1);
omdp = Kf./sum(Kf,1);  ompx = Kfx./sum(Kfx,1);  ompz = Kfz./sum(Kfz,1);

% get iterative pseudo-time steps
dtau_u = 1./( (1-fx).*Kvx./(h/ndim)^2 + Cvx + (1-fx).*(fx.^2./Cfx)./(h/ndim)^2 ) ;
dtau_w = 1./( (1-fz).*Kvz./(h/ndim)^2 + Cvz + (1-fz).*(fz.^2./Cfz)./(h/ndim)^2 ) ;
dtau_p = 1./( (1-f ).*Kf ./(h/ndim)^2 + Cf  + (1-f ).*(f .^2./Cv )./(h/ndim)^2 ) ;
dtau_f = dt/2.*ones(size(f)); % [s]

dtau_ustar = 1./( sum(Kvx./(h/ndim)^2 + fx.^2./Cfx./(h/ndim)^2 , 1));
dtau_wstar = 1./( sum(Kvz./(h/ndim)^2 + fz.^2./Cfz./(h/ndim)^2 , 1));
dtau_pstar = 1./( sum(Kf ./(h/ndim)^2 + f .^2./Cv ./(h/ndim)^2 , 1));
