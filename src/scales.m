
function [PHS] = scales (PHS,NUM)

% get diffusivity contrasts
kv = PHS.eta0;              % momentum diffusivity
kf = PHS.d0.^2./PHS.eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% get permission weights
F  = permute(repmat(PHS.f0,1,1,1,NUM.NPHS),[4,1,2,3]);
Sf = (F./PHS.B).^(1./PHS.C);  Sf = Sf./sum(Sf,2);
Xf = sum(PHS.A.*Sf,2).*F + (1-sum(PHS.A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));  gmv = geomean(geomean(thtv,3),2);
thtf = squeeze(prod(Mf.^Xf,2));  gmf = geomean(geomean(thtf,3),2);
thtv = 1./(1./thtv + 1./(NUM.thtlim.^0.5*gmv)) + (gmv/NUM.thtlim.^0.5);
thtf = 1./(1./thtf + 1./(NUM.thtlim.^0.5*gmf)) + (gmf/NUM.thtlim.^0.5);

% get momentum and volume flux and transfer coefficients
Cv = PHS.f0.*(1-PHS.f0).*PHS.eta0.*thtv./PHS.d0.^2;
Cf = PHS.f0.*(1-PHS.f0).*          thtf./PHS.eta0 ;

% apply cutoff to coefficients to safeguard numerical stability
Cv = 1./(1./Cv + 1./(min(geomean(geomean(Cv,3),2)).*NUM.cfflim));
Cf = 1./(1./Cf + 1./(min(geomean(geomean(Cf,3),2)).*NUM.cfflim));

PHS.Cv0 = Cv;

% get segregtaion-compaction length scales
% Cv of segregating phase, Cf of compacting phase
delta0 = PHS.f0.*PHS.f0'./sqrt(Cv.*Cf');
delta0 = delta0 - diag(diag(delta0));

% extract dominant length scales
n=0;
for i=1:NUM.NPHS
    for j=i:NUM.NPHS
        if i~=j
        n=n+1;
        [d0(n),ind] = max([delta0(i,j),delta0(j,i)]);
        if ind==1; d0_segr(n) = i; else; d0_segr(n) = j; end
        if ind==1; d0_cmpt(n) = j; else; d0_cmpt(n) = i; end
        end
    end
end

[PHS.delta0,ind] = unique(d0);
PHS.delta0_segr = d0_segr(ind);
PHS.delta0_cmpt = d0_cmpt(ind);

% get Darcy speed scales
% DeltaRho0 = abs(rho0 - rho0.');
PHS.DeltaRho0 = PHS.rho0 - sum(PHS.f0.*PHS.rho0,1);
PHS.w0        = PHS.DeltaRho0.*PHS.grav(1).*PHS.f0.^2./Cv;
PHS.p0        = PHS.DeltaRho0.*PHS.grav(1).*delta0;

end
