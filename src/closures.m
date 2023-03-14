
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
% dtau_u  = vfac./( (1+fx).*Kvx./Cvx./(h/2)^2 + (1-omvx) + (1+fx).*(fx.^2./Cfx./Cvx).*(1-omvx)./(h/2)^2 ) ;
% dtau_w  = vfac./( (1+fz).*Kvz./Cvz./(h/2)^2 + (1-omvz) + (1+fz).*(fz.^2./Cfz./Cvz).*(1-omvz)./(h/2)^2 ) ;
% dtau_p  = pfac./( (1+f ).*Kf ./Cf ./(h/2)^2 + (1-omfc) + (1+f ).*(f .^2./Cv ./Cf ).*(1-omfc)./(h/2)^2 ) ;
dtau_u  = 1./( (1+fx).*Kvx./(h/2)^2 + Cvx + (1+fx).*(fx.^2./Cfx)./(h/2)^2 ) ;
dtau_w  = 1./( (1+fz).*Kvz./(h/2)^2 + Cvz + (1+fz).*(fz.^2./Cfz)./(h/2)^2 ) ;
dtau_p  = 1./( (1+f ).*Kf ./(h/2)^2 + Cf  + (1+f ).*(f .^2./Cv )./(h/2)^2 ) ;

dtau_f = dt/2.*ones(size(f)); % [s]

% dtau_u  = 1./( (1+fx).*Kvx./(h/2)^2 + Cvx + (1+fx).*(fx.^2./Cfx)./(h/2)^2 ) ;
% dtau_w  = 1./( (1+fz).*Kvz./(h/2)^2 + Cvz + (1+fz).*(fz.^2./Cfz)./(h/2)^2 ) ;
% dtau_p  = 1./( (1+f ).*Kf ./(h/2)^2 + Cf  + (1+f ).*(f .^2./Cv )./(h/2)^2 ) ;

dtau_ustar = vfac./( sum((1+fx).*Kvx,1)./(h/2).^2 + sum((1+fx).*fx.^2./Cfx ,1)./(h/2)^2);
dtau_wstar = vfac./( sum((1+fz).*Kvz,1)./(h/2).^2 + sum((1+fz).*fz.^2./Cfz ,1)./(h/2)^2);
dtau_pstar = pfac./( sum((1+f ).*Kf ,1)./(h/2).^2 + sum((1+f ).*f .^2./Cv  ,1)./(h/2)^2);


% dtau_u  = dtau_f1./( 1./(1-fx).*d0.^2./(h/2)^2 + 1 + fx.^2./Cfx./Cvx./(h/2).^2);
% dtau_w  = dtau_f1./( 1./(1-fz).*d0.^2./(h/2)^2 + 1 + fz.^2./Cfz./Cvz./(h/2).^2);
% dtau_p  = dtau_f1./( 1./(1-f ).*d0.^2./(h/2)^2 + 1 + f .^2./Cv ./Cf ./(h/2).^2);
% 
% dtau_ustar  = dtau_f2./sum( 1./(1-fx).*d0.^2./(h/2)^2 + fx.^2./Cfx./Cvx./(h/2).^2, 1);
% dtau_wstar  = dtau_f2./sum( 1./(1-fz).*d0.^2./(h/2)^2 + fz.^2./Cfz./Cvz./(h/2).^2, 1);
% dtau_pstar  = dtau_f2./sum( 1./(1-f ).*d0.^2./(h/2)^2 + f .^2./Cv ./Cf ./(h/2).^2, 1);