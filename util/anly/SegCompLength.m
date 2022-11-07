function [dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength (f, eta0, d0, A, B, C, thtlim, cfflim)
% 
% [dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C, thtlim, cfflim)
% 
% Calculates the segregation-compaction length of a NPHS system given phase
% fractions, material properties, permission calibrations and limiters for
% numerical stability.
% 
% INPUTS
% f         phase fractions [NPHS x N]
% eta0      pure-phase viscosities [NPHS x 1]
% d0        pure-phase grain size [NPHS x 1]
% A, B, C   permission calibration fitting params [NPHS x NPHS]
% thtlim    numerical limiter on permission functions [1]
% cfflim    numerical limiter for coefficients [1]
% 
% OUTPUTS
% dsc       segregation-compaction length [NPHS x NPHS x N]
% Kv, Kf    momentum and volume flux coeffs [NPHS x N]
% Cv, Cf    momentum and volume transfer coeffs [NPHS x N]
% Xf        permission weights [NPHS x NPHS x N]
% 
% YQW, 12 Jan 2021


NPHS = size(f,1);

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% get permission weights
F  = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
Sf = (F./B).^(1./C);  Sf = Sf./sum(Sf,2);
Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));
thtf = squeeze(prod(Mf.^Xf,2));

if nargin>6
    gmv = geomean(geomean(thtv,3),2, 'omitnan');
    gmf = geomean(geomean(thtf,3),2, 'omitnan');
    thtv = 1./(1./thtv + 1./(thtlim.^0.5*gmv)) + (gmv/thtlim.^0.5);
    thtf = 1./(1./thtf + 1./(thtlim.^0.5*gmf)) + (gmf/thtlim.^0.5);
end

% get momentum and volume flux and transfer coefficients
Kv =    f .*eta0       .*thtv;
Kf =    f .*d0.^2./eta0.*thtf;
Cv = (1-f)./d0.^2.*Kv;
Cf = (1-f)./d0.^2.*Kf;

if nargin>7
    % apply cutoff to coefficients to safeguard numerical stability
    Kv = Kv + max(Kv,[],1,'omitnan')./cfflim;
    Kf = Kf + max(Kf,[],1,'omitnan')./cfflim;
    Cv = 1./(1./Cv + 1./(min(Cv,[],1).*cfflim));
    Cf = 1./(1./Cf + 1./(min(Cf,[],1).*cfflim));
end

dsc = zeros(NPHS,NPHS,size(f,2));

for k = 1:NPHS
    for i = 1:NPHS
        if i==k, continue; end
        dsc(k,i,:) = f(i,:).*f(k,:)./sqrt(Cv(i,:).*Cf(k,:));
    end
end

end