
%% prepare folder and associated files

% prepare workspace
if svop && ~exist([outdir RunID],'dir'); mkdir([outdir RunID]); end

% make a log file
if svop
    logfile = [outdir,'/',RunID,'/',RunID,'.log'];
    if exist(logfile,'file'); delete(logfile); end
    diary(logfile)
end

load('ocean.mat','ocean');

% print header
fprintf(1,'\n\n');
fprintf(1,'**************************************************************\n');
fprintf(1,'**********  pantarhei | multi-phase flow simulator  **********\n');
fprintf(1,'**************************************************************\n\n');

% check if running verification using mms
if ~exist('mms','var'), mms = false; end
if (mms), addpath('../mms/mms_utils/'); fprintf(1, 'Running MMS...\n\n'); end

% print some information on the run
fprintf(1, 'NPHS\t= %d \n', NPHS);
fprintf(1, 'f0  \t= [  '); fprintf(1, '%.2f  ', f0); fprintf('] with ');
if dfg~=0, fprintf(1, 'gaussian '); end
if dfr~=0, fprintf(1, 'random '); end
if exist('fInit','var'), fprintf(1, 'fInit '); end
fprintf(1, '\n');
fprintf(1, 'N   \t= %.d \t\t D \t= %.0f x max(delta)\n', N, Dfac);
fprintf(1, 'alpha\t= %.2f \t\t beta \t= %.2f\n', alpha, beta);


%% calculate material properties,
% reset nondimensional variables to dimensional values

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% check problem scales (returns delta0, w0)
[delta0, w0] = scales(f0, grav, rho0, eta0, d0, A, B, C, thtlim, cfflim);

% reset domain depth to multiple of max segr-comp-length
D  = Dfac.*max(delta0(:));
h  = D(1)/N;

% reset shear rate to multiple of max segr velocity scale
[w0max,tmp] = max(abs(w0(:)));
[~   ,iphs] = ind2sub([NPHS,NPHS],tmp); % index of segregating phase
Pu = Pu * w0max / f0(iphs) / D(1);
Si = Si * w0max / f0(iphs) / D(1);

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));


%% initialise coordinate arrays and BCs

if length(D)==1, D = [D,D]; end
z     = -D(1)/2+h/2 : h : D(1)/2-h/2;
x     = -D(2)/2+h/2 : h : D(2)/2-h/2;
if all([dfg;dfr]==0), x = 0; end % use 1D domain if no phase frac variations
[Z,X] = meshgrid(z,x);
Z     = permute(repmat(Z,1,1,NPHS),[3,2,1]);
X     = permute(repmat(X,1,1,NPHS),[3,2,1]);
Nz    = length(z);
Nx    = length(x);

% initialise indexing for boundary condition stencils (order: {zBC, xBC})
if ~iscell(BC), BC = {BC, BC}; end

% first, z direction
if strcmp(BC{1},'periodic')
    icz = [Nz,1:Nz,1]; imz = [Nz,1:Nz]; ipz = [1:Nz,1];
elseif strcmp(BC{1},'open') || strcmp(BC{1},'closed')
    icz = [1,1:Nz,Nz]; imz = [1,1:Nz]; ipz = [1:Nz,Nz];
end

% now x direction
if strcmp(BC{2},'periodic')
    icx = [Nx,1:Nx,1]; imx = [Nx,1:Nx]; ipx = [1:Nx,1];
elseif strcmp(BC{2},'open') || strcmp(BC{2},'closed')
    icx = [1,1:Nx,Nx]; imx = [1,1:Nx]; ipx = [1:Nx,Nx];
end

%% intialise phase fraction

if exist('fInit','var')  
    f = fInit(X, Z, D);
else
    % initialise smoothed random and gaussian perturbation fields
    rng(15);
    rnd = randn(NPHS,Nz,Nx);
    for i = 1:smth
        rnd = rnd + diff(rnd(:,icz,:),2,2)./12 + diff(rnd(:,:,icx),2,3)./12;
    end
    rnd = rnd./max(abs(rnd(:)));
    gsn = exp(-X.^2./(D(2)/8).^2).*exp(-Z.^2./(D(1)/8).^2);
    f   = f0 + dfr.*rnd + dfg.*gsn;
end

f = max(1e-16,min(1-1e-16,f)); f = f./sum(f,1); 

%% initialise other fields

rho    = rho0.*ones(size(f));
rhomix = mean(mean(sum(f.*rho,1)));
if (mms), mms_init_phasefrac; end
if svop && restart==0, save([outdir RunID,'/',RunID,'_par.mat']); end
if (mms), mms_plot_truesol; end

fo = f;  fi = f;     res_f = 0*f;  dtau_f = res_f;
u = zeros(NPHS,Nz  ,Nx+1);  ui = u;  ustar = mean(u,1);  usegr = 0*u;  res_u = 0*u;  dtau_u = res_u;
w = zeros(NPHS,Nz+1,Nx  );  wi = w;  wstar = mean(w,1);  wsegr = 0*w;  res_w = 0*w;  dtau_w = res_w;
p = zeros(NPHS,Nz  ,Nx  );  pi = p;  pstar = mean(p,1);  pcmpt = 0*p;  res_p = 0*p;  dtau_p = res_p;

% shearing velocities on the staggered grid
ushr   = z'*Si - [x-0.5*h,x(end)+0.5*h] *Pu;    ushr = permute(ushr,[3,1,2]);
wshr   = x *Si + [z-0.5*h,z(end)+0.5*h]'*Pu;    wshr = permute(wshr,[3,1,2]);

% initialise auxiliary fields
qvxx   = zeros(NPHS,Nz,Nx  );  qvzz = zeros(NPHS,Nz,Nx  );  qvxz = zeros(NPHS,Nz+1,Nx+1);
qfx    = zeros(NPHS,Nz,Nx+1);  qfz  = zeros(NPHS,Nz+1,Nx);
Gvx    = zeros(NPHS,Nz,Nx+1);  Gvz  = zeros(NPHS,Nz+1,Nx);
Gf     = zeros(NPHS,Nz,Nx  );  Gfo  = Gf; 
fadv   = zeros(NPHS,Nz,Nx  );  fadvo= fadv;
delta  = zeros(NPHS,NPHS,Nz,Nx);
