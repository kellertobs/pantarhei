
%% prepare folder and associated files

% prepare workspace
if IO.svop && ~exist([IO.outdir IO.RunID],'dir'); mkdir([IO.outdir IO.RunID]); end

% make a log file
if IO.svop
    logfile = [IO.outdir,'/',IO.RunID,'/',IO.RunID,'.log'];
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
NUM.TINY = 1e-16;
NUM.HUGE = 1e+16;
NUM.pltits = IO.pltits;

% some plotting checks
if ~isfield(IO,'pltits'), IO.pltits = false; end

%% calculate material properties
% reset nondimensional variables to dimensional values

% get diffusivity contrasts
kv     = PHS.eta0;          % momentum diffusivity
kf     = PHS.d0.^2./PHS.eta0;   % volume diffusivity
PHS.Mv = kv.'./kv;      % momentum diffusivity ratios
PHS.Mf = kf.'./kf;      % volume diffusivity ratios

% check problem scales (returns delta0, w0)
[PHS] = scales(PHS,NUM);

% reset domain depth to multiple of max segr-comp-length
if length(NUM.Lfac)==1, NUM.Lfac = [NUM.Lfac, NUM.Lfac]; end
NUM.L  = NUM.Lfac.*max(PHS.delta0(:)); NUM.Lmax = max(NUM.L);
NUM.h  = NUM.L(1)/NUM.N;

% reset shear rate to multiple of max segr velocity scale
[w0max,tmp] = max(abs(PHS.w0(:)));
[~   ,iphs] = ind2sub([NUM.NPHS,NUM.NPHS],tmp); % index of segregating phase
PHS.Pu = PHS.Pu * w0max / PHS.f0(iphs) / NUM.L(1);
PHS.Si = PHS.Si * w0max / PHS.f0(iphs) / NUM.L(1);

% set appropriate initial time step size
dt = NUM.cfl.*NUM.h/2/max(PHS.w0(:));

%% initialise coordinate arrays and BCs

NUM.z     = (-NUM.L(1)/2+NUM.h/2) : NUM.h : (NUM.L(1)/2-NUM.h/2);
NUM.x     = (-NUM.L(2)/2+NUM.h/2) : NUM.h : (NUM.L(2)/2-NUM.h/2);
if all([PHS.dfg;PHS.dfr]==0), NUM.x = 0; end % use 1D domain if no phase frac variations

[Z,X]     = meshgrid(NUM.z,NUM.x);
NUM.Z     = permute(repmat(Z,1,1,NUM.NPHS),[3,2,1]);
NUM.X     = permute(repmat(X,1,1,NUM.NPHS),[3,2,1]);
NUM.Nz    = length(NUM.z);
NUM.Nx    = length(NUM.x);
NUM.Nmax  = max(NUM.Nz, NUM.Nx);
NUM.ndim  = sum(double([NUM.Nx;NUM.Nz]>1));

NUM.zw    = [NUM.z-0.5*NUM.h,NUM.z(end)+0.5*NUM.h];  % z coordinates for z velocity (horiz    cell faces)
NUM.xu    = [NUM.x-0.5*NUM.h,NUM.x(end)+0.5*NUM.h];  % x coordinates for x velocity (vertical cell faces)

if size(NUM.Z,3) > 1
    NUM.maxlvl = log2(min([size(NUM.Z,2),size(NUM.Z,3)]));
else
    NUM.maxlvl = log2(size(NUM.Z,2));
end

% initialise indexing for boundary condition stencils (order: {zBC, xBC})
if ~iscell(NUM.BC), NUM.BC = {NUM.BC, NUM.BC}; end

% create index arrays
get_indices;

% checks for direct solver
if NUM.ndim>1; NUM.direct = false; end
if NUM.direct; NUM.dmp = 0; end

%% intialise phase fraction

if exist('fInit','var')  
    f = fInit(NUM.X, NUM.Z, NUM.L);
else
    % initialise smoothed random and gaussian perturbation fields
    rng(15);
    rnd = randn(NUM.NPHS,NUM.Nz,NUM.Nx);
    for i = 1:PHS.smth
        rnd = rnd + diff(rnd(:,NUM.icz,:),2,2)./12 + diff(rnd(:,:,NUM.icx),2,3)./12;
    end
    rnd = rnd./max(abs(rnd(:)));
    gsn = exp(-NUM.X.^2./(NUM.L(2)/12).^2).*exp(-NUM.Z.^2./(NUM.L(1)/12).^2);
    f   = PHS.f0 + PHS.dfr.*rnd + PHS.dfg.*gsn;
end

%% initialise other fields

if (mms), mms_init_phasefrac; end
if IO.svop && IO.restart==0, save([IO.outdir IO.RunID,'/',IO.RunID,'_par.mat']); end
if (mms), mms_plot_truesol; end

f = max(NUM.flim,min(1-NUM.flim,f)); f = f./sum(f,1);  fo = f;  fin = f;
u = zeros(NUM.NPHS,NUM.Nz  ,NUM.Nx+1);  ustar = mean(u,1);  
w = zeros(NUM.NPHS,NUM.Nz+1,NUM.Nx  );  wstar = mean(w,1);  
p = zeros(NUM.NPHS,NUM.Nz  ,NUM.Nx  );  pstar = mean(p,1);  

if exist('wInit', 'var'), [w,wstar] = wInit(w, PHS.w0); end
if exist('pInit', 'var'), [p,pstar] = pInit(p, PHS.p0); end

% initialise auxiliary fields
fi = f;                          res_f = 0*f;  dtau_f = 0*f;                   upd_f = 0*f;
ui = u;  usegr = 0*u;  Du = 0*u; res_u = 0*u;  dtau_u = 0*u; dtau_ustar = 0*u; upd_u = 0*u;
wi = w;  wsegr = 0*w;  Dw = 0*w; res_w = 0*w;  dtau_w = 0*w; dtau_wstar = 0*w; upd_w = 0*w;
pi = p;  pcmpt = 0*p;  Dp = 0*p; res_p = 0*p;  dtau_p = 0*p; dtau_pstar = 0*p; upd_p = 0*p;

qvxx   = zeros(NUM.NPHS,NUM.Nz,NUM.Nx  );  qvzz = zeros(NUM.NPHS,NUM.Nz  ,NUM.Nx);  qvxz = zeros(NUM.NPHS,NUM.Nz+1,NUM.Nx+1);
qfx    = zeros(NUM.NPHS,NUM.Nz,NUM.Nx+1);  qfz  = zeros(NUM.NPHS,NUM.Nz+1,NUM.Nx);
Gvx    = zeros(NUM.NPHS,NUM.Nz,NUM.Nx+1);  Gvz  = zeros(NUM.NPHS,NUM.Nz+1,NUM.Nx);
Gf     = zeros(NUM.NPHS,NUM.Nz,NUM.Nx  );  Gfo  = Gf; 
dfdt   = zeros(NUM.NPHS,NUM.Nz,NUM.Nx  );  dfdto = dfdt;
delta  = zeros(NUM.NPHS,NUM.NPHS,NUM.Nz,NUM.Nx);
delmax = zeros(NUM.NPHS,NUM.Nz,NUM.Nx  );

rho    = PHS.rho0.*ones(size(f));
rhomix = mean(mean(sum(f.*rho,1)));

% shearing velocities on the staggered grid
ushr   = NUM.z'*PHS.Si - NUM.xu *PHS.Pu;    ushr = permute(ushr,[3,1,2]);
wshr   = NUM.x *PHS.Si + NUM.zw'*PHS.Pu;    wshr = permute(wshr,[3,1,2]);

%% print some information on the run

if NUM.Nx==1, fprintf(1, 'Running 1D model'); 
else    , fprintf(1, 'Running 2D model'); end
fprintf(1, ' with RunID [ %s ] .\n\n', IO.RunID);

fprintf(1, '    NPHS = %d \n', NUM.NPHS);
fprintf(1, '      f0 = [ '); fprintf(1, '%.4f ', PHS.f0); fprintf(']\n');
fprintf(1, '  delta0 = [ '); fprintf(1, '%.4f ', PHS.delta0); fprintf(']\n');
fprintf(1, '    segr = [ '); fprintf(1, '    %d  ', PHS.delta0_segr); fprintf(']\n');
fprintf(1, '    cmpt = [ '); fprintf(1, '    %d  ', PHS.delta0_cmpt); fprintf(']');

if PHS.dfg~=0, fprintf(1, 'gaussian '); end
if PHS.dfr~=0, fprintf(1, 'random '); end
if exist('fInit','var'), fprintf(1, 'fInit '); end
fprintf(1, '\n\n');
fprintf(1, '[Nz, Nx] = [  %d  %d ]  \n', NUM.Nz, NUM.Nx);
fprintf(1, '[Lz, Lx] = [   %.0f   %.0f ] x max(delta)\n\n', NUM.Lfac(1), NUM.Lfac(2).*(NUM.Nx>1));
