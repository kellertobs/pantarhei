
function [delta0, w0] = scales (f0, grav, rho0, eta0, d0, A, B, C, thtlim, cfflim)

NPHS = length(f0);

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
Cv = f0.*(1-f0).*eta0.*thtv./d0.^2;
Cf = f0.*(1-f0).*      thtf./eta0 ;

% apply cutoff to coefficients to safeguard numerical stability
Cv = 1./(1./Cv + 1./(min(geomean(geomean(Cv,3),2)).*cfflim));
Cf = 1./(1./Cf + 1./(min(geomean(geomean(Cf,3),2)).*cfflim));

% get segregtaion-compaction length scales
% Cv of segregating phase, Cf of compacting phase
delta0 = f0.'.*f0./sqrt(Cv.'.*Cf);
delta0 = delta0 - diag(diag(delta0));

% get Darcy speed scales
% DeltaRho0 = abs(rho0 - rho0.');
DeltaRho0 = rho0 - sum(f0.*rho0,1);
w0        = DeltaRho0.*max(abs(grav)).*f0.^2./Cv;

end
