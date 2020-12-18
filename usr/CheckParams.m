
% use this script to check that parameters applied make sense to problem.


% check that initial dt is similar to timescale of variation in problem
save vars.mat
[dtOut] = CheckTimeStep('vars.mat')

% check the segregation compaction length
[dsc] = SegCompLength(f0, eta0, d0, A, B, C, thtlim, cfflim)


 
%%

function [dsc] = SegCompLength (f, eta0, d0, A, B, C, thtlim, cfflim)

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

dsc = f.'.*f./sqrt(Cv.'.*Cf);
dsc = dsc - diag(diag(dsc));

end

function [dtOut] = CheckTimeStep (MatFileName)
% calculate dt as expected by Courant Flow Condition

load(MatFileName);

dt     = 1;
tend   = 0.2;
% minits = 20;
% maxits = 20;

run('../src/pantarhei');

dtOut= cfl/(max(abs([qfx(:);qfz(:)]))/(h/2) + max(abs(Gf(:)))./5e-3);  % [s]

end
