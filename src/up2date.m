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

% apply cutoff to coefficients to safeguard numerical stability
Kv = Kv + max(Kv,[],1)./cfflim;
Kf = Kf + max(Kf,[],1)./cfflim;
Cv = 1./(1./Cv + 1./(min(Cv,[],1).*cfflim));
Cf = 1./(1./Cf + 1./(min(Cf,[],1).*cfflim));

% get coefficient-based weights
omvc =  Cv./sum(Cv,1);
omvx = (Cv(:,:,im)+Cv(:,:,ip))./2./sum((Cv(:,:,im)+Cv(:,:,ip))./2,1);
omvz = (Cv(:,im,:)+Cv(:,ip,:))./2./sum((Cv(:,im,:)+Cv(:,ip,:))./2,1);
omfc =  Cf./sum(Cf,1);
omfx = (Kf(:,:,im)+Kf(:,:,ip))./2./sum((Kf(:,:,im)+Kf(:,:,ip))./2,1);
omfz = (Kf(:,im,:)+Kf(:,ip,:))./2./sum((Kf(:,im,:)+Kf(:,ip,:))./2,1);

% get phase advection term
ustar     = sum(omvx.*u,1);
wstar     = sum(omvz.*w,1);
vstar_Gf  = diff((f(:,:,im)+f(:,:,ip))./2.*ustar,1,3)./h ...
          + diff((f(:,im,:)+f(:,ip,:))./2.*wstar,1,2)./h ...
          - f.*(diff(ustar,1,3)./h + diff(wstar,1,2)./h);

% get iterative pseudo-time steps
dtau_u = min(1.0./((Kv(:,:,im)+Kv(:,:,ip))./2./(h/2)^2 + (1-omvx).*(Cv(:,:,im)+Cv(:,:,ip))/2),[],1).*ones(size(u));  % [Pas/m2]
dtau_w = min(1.0./((Kv(:,im,:)+Kv(:,ip,:))./2./(h/2)^2 + (1-omvz).*(Cv(:,im,:)+Cv(:,ip,:))/2),[],1).*ones(size(w));  % [Pas/m2]
dtau_p = min(1.0./( Kf                       ./(h/2)^2 + (1-omfc).* Cf                      ),[],1).*ones(size(p));  % [1/Pas]
dtau_f = dt/10.*ones(size(f));  % [s]
