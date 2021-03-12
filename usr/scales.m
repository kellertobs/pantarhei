% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% get permission weights
F  = permute(repmat(f0,1,1,1,NPHS),[4,1,2,3]);
Sf = (F./B).^(1./C);  Sf = Sf./sum(Sf,2);
Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));  gmv = geomean(geomean(thtv,3),2);
thtf = squeeze(prod(Mf.^Xf,2));  gmf = geomean(geomean(thtf,3),2);
thtv = 1./(1./thtv + 1./(thtlim.^0.5*gmv)) + (gmv/thtlim.^0.5);
thtf = 1./(1./thtf + 1./(thtlim.^0.5*gmf)) + (gmf/thtlim.^0.5);

% get momentum and volume flux and transfer coefficients
Kv =    f0 .*eta0       .*thtv;
Kf =    f0 .*d0.^2./eta0.*thtf;
Cv = (1-f0)./d0.^2.*Kv;
Cf = (1-f0)./d0.^2.*Kf;

% apply cutoff to coefficients to safeguard numerical stability
Kv = Kv + max(Kv,[],1)./cfflim;
Kf = Kf + max(Kf,[],1)./cfflim;
Cv = 1./(1./Cv + 1./(min(Cv,[],1).*cfflim));
Cf = 1./(1./Cf + 1./(min(Cf,[],1).*cfflim));

% get segregtaion-compaction length scales
delta0 = f0.*f0.'./sqrt(Cv.*Cf.');
delta0 = delta0 - diag(diag(delta0));

% get Darcy speed scales
DeltaRho0 = abs(rho0 - rho0.');
w0        = DeltaRho0.*max(abs(grav)).*f0.^2./Cv;