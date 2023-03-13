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

% get segregation-compaction lengths
for i=1:NPHS
    for k=1:NPHS
        if i~=k; delta(i,k,:,:) = f(i,:,:).*f(k,:,:)./sqrt(Cv(i,:,:).*Cf(k,:,:)); end
    end
end
fsusp = ones(size(f));%squeeze(d0./max(max(delta))).^2;

% get phase fractions and coefficients on faces
fx   = (f (:,:,imx)+f (:,:,ipx))./2;
fz   = (f (:,imz,:)+f (:,ipz,:))./2;
rhox = (rho(:,:,imx)+rho(:,:,ipx))./2;
rhoz = (rho(:,imz,:)+rho(:,ipz,:))./2;
Kvx  = (Kv(:,:,imx)+Kv(:,:,ipx))./2;
Kvz  = (Kv(:,imz,:)+Kv(:,ipz,:))./2;
Kvc  = (Kv(:,imz,imx)+Kv(:,ipz,imx)+Kv(:,imz,ipx)+Kv(:,ipz,ipx))./4;
Kfx  = (Kf(:,:,imx)+Kf(:,:,ipx))./2;
Kfz  = (Kf(:,imz,:)+Kf(:,ipz,:))./2;
Cvx  = (Cv(:,:,imx)+Cv(:,:,ipx))./2;
Cvz  = (Cv(:,imz,:)+Cv(:,ipz,:))./2;
Cfx  = (Cf(:,:,imx)+Cf(:,:,ipx))./2;
Cfz  = (Cf(:,imz,:)+Cf(:,ipz,:))./2;

% get coefficient-based weights
omvc =  Cv./sum(Cv,1);
omvx = (Cv(:,:,imx)+Cv(:,:,ipx))./2./sum((Cv(:,:,imx)+Cv(:,:,ipx))./2,1);
omvz = (Cv(:,imz,:)+Cv(:,ipz,:))./2./sum((Cv(:,imz,:)+Cv(:,ipz,:))./2,1);
omfc =  Cf./sum(Cf,1);
omfx = (Kf(:,:,imx)+Kf(:,:,ipx))./2./sum((Kf(:,:,imx)+Kf(:,:,ipx))./2,1);
omfz = (Kf(:,imz,:)+Kf(:,ipz,:))./2./sum((Kf(:,imz,:)+Kf(:,ipz,:))./2,1);
omdp =  Kf./sum(Kf,1);

% get phase advection term
vr_Gf  = advect(f, ur, wr, h, {advn, 'vdf'}, [2,3], BC);


% get iterative pseudo-time steps
dtau_u  = 1./((1+fx).*(Kvx)./(h/2)^2 + 1 + (1+fx).*(fx.^2./Cfx)./(h/2)^2);
dtau_w  = 1./((1+fz).*(Kvz)./(h/2)^2 + 1 + (1+fz).*(fz.^2./Cfz)./(h/2)^2);
dtau_p  = 1./((1+f ).*(Kf )./(h/2)^2 + 1 + (1+f ).*(f .^2./Cv )./(h/2)^2);

dtau_ur = vfact./sum((1+fx).*Kvx./(h/2)^2 + (1-omvx).*Cvx + (1+fx).*fx.^2./Cfx./(h/2)^2);
dtau_wr = vfact./sum((1+fz).*Kvz./(h/2)^2 + (1-omvz).*Cvz + (1+fz).*fz.^2./Cfz./(h/2)^2);
dtau_pr = pfact./sum((1+f ).*Kf ./(h/2)^2 + (1-omfc).*Cf  + (1+f ).*f .^2./Cv ./(h/2)^2);

dtau_f = dt/2.*ones(size(f)); % [s]
