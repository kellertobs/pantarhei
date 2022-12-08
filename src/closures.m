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

% get coefficient-based weights
omvc =  Cv./sum(Cv,1);
omvx = (Cv(:,:,imx)+Cv(:,:,ipx))./2./sum((Cv(:,:,imx)+Cv(:,:,ipx))./2,1);
omvz = (Cv(:,imz,:)+Cv(:,ipz,:))./2./sum((Cv(:,imz,:)+Cv(:,ipz,:))./2,1);
omfc =  Cf./sum(Cf,1);
omfx = (Kf(:,:,imx)+Kf(:,:,ipx))./2./sum((Kf(:,:,imx)+Kf(:,:,ipx))./2,1);
omfz = (Kf(:,imz,:)+Kf(:,ipz,:))./2./sum((Kf(:,imz,:)+Kf(:,ipz,:))./2,1);
omdp =  Kf./sum(Kf,1);

% get phase advection term
ustar    = sum(omvx.*u,1);
wstar    = sum(omvz.*w,1);
vstar_Gf = advect(f, ustar, wstar, h, {advn, 'vdf'}, [2,3], BC);

% get iterative pseudo-time steps
dtau_u = min(1.0./abs(-(Kv(:,:,imx)+Kv(:,:,ipx))./2./(h/2)^2 + (1-omvx).*(Cv(:,:,imx)+Cv(:,:,ipx))/2),[],1).*ones(size(u));  % [Pas/m2]
dtau_w = min(1.0./abs(-(Kv(:,imz,:)+Kv(:,ipz,:))./2./(h/2)^2 + (1-omvz).*(Cv(:,imz,:)+Cv(:,ipz,:))/2),[],1).*ones(size(w));  % [Pas/m2]
dtau_p = min(1.0./abs(- Kf                         ./(h/2)^2 + (1-omfc).* Cf                        ),[],1).*ones(size(p));  % [1/Pas]
dtau_f = dt/10.*ones(size(f));  % [s]

% dtau_u = 1.0./abs(-(Kv(:,:,imx)+Kv(:,:,ipx))./2./(h/2)^2 + (1-omvx).*(Cv(:,:,imx)+Cv(:,:,ipx))/2);  % [Pas/m2]
% dtau_w = 1.0./abs(-(Kv(:,imz,:)+Kv(:,ipz,:))./2./(h/2)^2 + (1-omvz).*(Cv(:,imz,:)+Cv(:,ipz,:))/2);  % [Pas/m2]
% dtau_p = 0.5./abs( -Kf.*(1-omdp)               ./(h/2)^2 + (1-omfc).* Cf                        );  % [1/Pas]