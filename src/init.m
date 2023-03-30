
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

% stabilizing values
TINY = 1e-16;
HUGE = 1e+16;

% some plotting checks
figvis = 'off'; % toggle for figure visibility on/off
if ~exist('pltits','var'), pltits = false; end

%% calculate material properties
% reset nondimensional variables to dimensional values

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% check problem scales (returns delta0, w0)
[delta0, w0, p0, Cv0] = scales(f0, grav, rho0, eta0, d0, A, B, C, thtlim, cfflim);
[delta0_max, imax] = max(delta0(:));
[delmax_i,delmax_j] = ind2sub([NPHS, NPHS], imax);

% reset domain depth to multiple of max segr-comp-length
if length(Lfac)==1, Lfac = [Lfac, Lfac]; end
L  = Lfac.*max(delta0(:)); Lmax = max(L);
h  = L(1)/N;

% reset shear rate to multiple of max segr velocity scale
[w0max,tmp] = max(abs(w0(:)));
[~   ,iphs] = ind2sub([NPHS,NPHS],tmp); % index of segregating phase
Pu = Pu * w0max / f0(iphs) / L(1);
Si = Si * w0max / f0(iphs) / L(1);

% set appropriate initial time step size
dt = cfl.*h/2/max(w0(:));

%% initialise coordinate arrays and BCs

z     = (-L(1)/2+h/2) : h : (L(1)/2-h/2);
x     = (-L(2)/2+h/2) : h : (L(2)/2-h/2);
if all([dfg;dfr]==0), x = 0; end % use 1D domain if no phase frac variations

[Z,X] = meshgrid(z,x);
Z     = permute(repmat(Z,1,1,NPHS),[3,2,1]);
X     = permute(repmat(X,1,1,NPHS),[3,2,1]);
Nz    = length(z);
Nx    = length(x);
Nmax  = max(Nz, Nx);
ndim  = sum(double([Nx;Nz]>1));

zw    = [z-0.5*h,z(end)+0.5*h];  % z coordinates for z velocity (horiz    cell faces)
xu    = [x-0.5*h,x(end)+0.5*h];  % x coordinates for x velocity (vertical cell faces)

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

if ndim==1
    Re = 9*sqrt( 3)*pi/4; r = 0.25; 
elseif ndim==2
    Re = 3*sqrt(10)*pi/2; r = 0.5;
end
vfac = 0.9/Re/sqrt(ndim);
pfac = 0.9*Re/sqrt(ndim)/(r+2);

%% intialise phase fraction

if exist('fInit','var')  
    f = fInit(X, Z, L);
else
    % initialise smoothed random and gaussian perturbation fields
    rng(15);
    rnd = randn(NPHS,Nz,Nx);
    for i = 1:smth
        rnd = rnd + diff(rnd(:,icz,:),2,2)./12 + diff(rnd(:,:,icx),2,3)./12;
    end
    rnd = rnd./max(abs(rnd(:)));
    gsn = exp(-X.^2./(L(2)/12).^2).*exp(-Z.^2./(L(1)/12).^2);
    f   = f0 + dfr.*rnd + dfg.*gsn;
end

%% initialise other fields

if (mms), mms_init_phasefrac; end
if svop && restart==0, save([outdir RunID,'/',RunID,'_par.mat']); end
if (mms), mms_plot_truesol; end

f = max(flim,min(1-flim,f)); f = f./sum(f,1);   fo = f;  
u = zeros(NPHS,Nz  ,Nx+1);  ustar = mean(u,1);  
w = zeros(NPHS,Nz+1,Nx  );  wstar = mean(w,1);  
p = zeros(NPHS,Nz  ,Nx  );  pstar = mean(p,1);  

if exist('wInit', 'var'), [w,wstar] = wInit(w, w0); end
if exist('pInit', 'var'), [p,pstar] = pInit(p, p0); end

% initialise auxiliary fields
fi = f;                          res_f = 0*f;  dtau_f = 0*f;                   upd_f = 0*f;
ui = u;  usegr = 0*u;  Du = 0*u; res_u = 0*u;  dtau_u = 0*u; dtau_ustar = 0*u; upd_u = 0*u;
wi = w;  wsegr = 0*w;  Dw = 0*w; res_w = 0*w;  dtau_w = 0*w; dtau_wstar = 0*w; upd_w = 0*w;
pi = p;  pcmpt = 0*p;  Dp = 0*p; res_p = 0*p;  dtau_p = 0*p; dtau_pstar = 0*p; upd_p = 0*p;

qvxx   = zeros(NPHS,Nz,Nx  );  qvzz = zeros(NPHS,Nz  ,Nx);  qvxz = zeros(NPHS,Nz+1,Nx+1);
qfx    = zeros(NPHS,Nz,Nx+1);  qfz  = zeros(NPHS,Nz+1,Nx);
Gvx    = zeros(NPHS,Nz,Nx+1);  Gvz  = zeros(NPHS,Nz+1,Nx);
Gf     = zeros(NPHS,Nz,Nx  );  Gfo  = Gf; 
fadv   = zeros(NPHS,Nz,Nx  );  fadvo= fadv;
delta  = zeros(NPHS,NPHS,Nz,Nx);
delmax = zeros(NPHS,Nz,Nx  );

rho    = rho0.*ones(size(f));
rhomix = mean(mean(sum(f.*rho,1)));

% shearing velocities on the staggered grid
ushr   = z'*Si - xu *Pu;    ushr = permute(ushr,[3,1,2]);
wshr   = x *Si + zw'*Pu;    wshr = permute(wshr,[3,1,2]);

%% print some information on the run

if Nx==1, fprintf(1, 'Running 1D model'); 
else    , fprintf(1, 'Running 2D model'); end
fprintf(1, ' with RunID [ %s ] .\n\n', RunID);

fprintf(1, '    NPHS = %d \n', NPHS);
fprintf(1, '      f0 = [ '); fprintf(1, '%.2f ', f0); fprintf('] with ');
if dfg~=0, fprintf(1, 'gaussian '); end
if dfr~=0, fprintf(1, 'random '); end
if exist('fInit','var'), fprintf(1, 'fInit '); end
fprintf(1, '\n\n');
fprintf(1, '[Nz, Nx] = [  %d  %d ]  \n', Nz, Nx);
fprintf(1, '[Lz, Lx] = [   %.0f   %.0f ] x max(delta)\n\n', Lfac(1), Lfac(2).*(Nx>1));
fprintf(1, '   alpha = %.2f \n', alpha);
fprintf(1, '   beta  = %.2f \n', beta);
