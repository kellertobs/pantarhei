
% print header
fprintf(1,'\n\n****************************************************\n');
fprintf(1,'*****  pantarhei | multi-phase flow simulator  *****\n');
fprintf(1,'****************************************************\n\n');

% prepare workspace
if svop && ~exist([outdir RunID],'dir'); mkdir([outdir RunID]); end
if svop && restart==0, save([outdir RunID,'/',RunID,'_par.mat']); end

% make a log file
if svop
    logfile = [outdir,'/',RunID,'/',RunID,'.log'];
    if exist(logfile,'file'); delete(logfile); end
    diary(logfile)
end

load('ocean.mat','ocean');

% check if running verification using mms
if ~exist('mms','var'), mms = false; end
if (mms), addpath('../mms/mms_utils/'); fprintf(1, 'Running MMS...\n'); end

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% initialise coordinate arrays
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

% intialise solution and residual fields. First phase fraction.
if exist('fInit','var')  
    f = fInit(X, Z);
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

f = max(1e-16,min(1-1e-16,f));  f = f./sum(f,1);  fo = f;  fi = f;     res_f = 0*f;  dtau_f = res_f;
u = zeros(NPHS,Nz  ,Nx+1);  ui = u;  ustar = mean(u,1);  usegr = 0*u;  res_u = 0*u;  dtau_u = res_u;
w = zeros(NPHS,Nz+1,Nx  );  wi = w;  wstar = mean(w,1);  wsegr = 0*w;  res_w = 0*w;  dtau_w = res_w;
p = zeros(NPHS,Nz  ,Nx  );  pi = p;  pstar = mean(p,1);  pcmpt = 0*p;  res_p = 0*p;  dtau_p = res_p;
if (mms), mms_init_phasefrac; end

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
rho    = rho0.*ones(size(f));
rhomix = mean(mean(sum(f.*rho,1)));

% initialise time stepping loop
time = 0;
step = 0;

% restart run?
if restart>0
    % load previous file to get ui, wi, pi, fi
    % NB: these are not exactly the same as (ui, wi, pi, fi) if nop~=1, 
    % but closer than zeros initialisation matrices
    fname = [outdir RunID,'/',RunID,'_',num2str(restart-1),'.mat'];
    load(fname); fprintf(1,'Loaded %s to get [u,w,p,f] of previous time step.\n', fname);
    ui = u; wi = w; pi = p;  fi = f;
    
    % load current file 
    fname = [outdir RunID,'/',RunID,'_',num2str(restart),'.mat'];
    load(fname); fprintf(1,'Loaded %s.\n', fname);
    rsstep = restart*nop;  % current step
    step   = rsstep + 1;
end

% initialise coefficient closures and constitutive relations
closures; constitutive;

while time <= tend && step <= NtMax  % keep stepping until final run time reached
    
    tic;
    
    % print time step diagnostic
    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    
    % store phase fractions of previous step
    fo = f;  dto = dt; Gfo = Gf; fadvo = fadv;
    
    % initialise non-linear iteration loop
    res  = 1e3;
    res0 = res;
    it   = 0;
    
    if step == 0 % require more stringent convergence criterion for first step
        conv_crit = @(resc,res0c,itc) (resc >= atol && itc <= 10*maxits);
    else
        conv_crit = @(resc,res0c,itc) (resc >= atol && resc/res0c >= rtol && itc <= maxits || itc < minits);
    end
    
    % make iterative convergence plot
    if (nop>0) && ~mod(step,abs(nop)), f10 = figure(10);  f10.Visible = 'off'; end
    
    while  conv_crit(res,res0,it) > 0 % keep stepping until convergence criterion reached
        
        % store solution of previous two iterations
        uii = ui;  ui = u;  wii = wi;  wi = w;  pii = pi;  pi = p;  fii = fi;  fi = f;
        
        it = it+1;  % update iteration count
        
        % update coefficient closures (only every 'nupd' iterations)
        if ~mod(it,nupd); closures; end
        
        % update constitutive relations
        constitutive;
        
        % update phase advection term when shear velocities are included
        fadv = advect(f, ushr, wshr, h, {advn, 'vdf'}, [2,3], BC);
           
        % update physical time step [s]
        dt   = min([ 2*dto; 
                     cfl/(max(abs([qfx(:);qfz(:)]+advscl))/(h/2) + max(abs(Gf(:)))./1e-2); 
                     cfl*0.5*h./max(abs([ushr(:);wshr(:)]) + 1e-16) ]);  

        % update residual fields
        res_u =             + diff(qvxx(:,:,icx),1,3)./h + diff(qvxz,1,2)./h + Gvx + Qvx    ;
        res_w =             + diff(qvzz(:,icz,:),1,2)./h + diff(qvxz,1,3)./h + Gvz + Qvz    ;
        res_p =             + dqf + Gf  + Gm./rho;
        res_f = (f-fo)./dt                                                   -(Gf + Gfo)./2 + (fadv+fadvo)./2;
        
        % call manufactured solution (if benchmarking)
        if (mms); mms_calc_source; end
        
        % do not update phase fractions during initial step
        if step==0; res_f = 0.*res_f; end
        
        % prepare solution updates
        if any( strcmp(BC,'closed') )
            if strcmp(BC{1},'closed'), res_w(:,[1,end],:) = 0; end
            if strcmp(BC{2},'closed'), res_u(:,:,[1,end]) = 0; end
            upd_u = res_u.*dtau_u;
            upd_w = res_w.*dtau_w;
            upd_p = res_p.*dtau_p;
            upd_f = res_f.*dtau_f;
        else
            upd_u = res_u.*dtau_u - mean(res_u(:).*dtau_u(:));
            upd_w = res_w.*dtau_w - mean(res_w(:).*dtau_w(:));
            upd_p = res_p.*dtau_p - mean(res_p(:).*dtau_p(:));
            upd_f = res_f.*dtau_f - mean(res_f(:).*dtau_f(:));
        end
        
        % update velocity-pressure solution
        u = ui - alpha.*upd_u + beta.*(ui-uii);
        w = wi - alpha.*upd_w + beta.*(wi-wii);
        p = pi - alpha.*upd_p + beta.*(pi-pii);
        
        % update phase fractions (only every 'nupd' iterations)
        if ~mod(it,nupd); f = fi - alpha.*upd_f; f = max(flim,min(1-flim,f)); end
        
        % print iteration diagnostics
        if ~mod(it,nupd) || it==1
            % get residual norm
            resflds = [ norm(res_u(:).*dtau_u(:),2)./(norm(u(:),2)+1e-32); 
                        norm(res_w(:).*dtau_w(:),2)./(norm(w(:),2)+1e-32);
                        norm(res_p(:).*dtau_p(:),2)./(norm(p(:),2)+1e-32);
                        norm(res_f(:).*dtau_f(:),2)./(norm(f(:),2)+1e-32)];
            res = sum(resflds);
            if res>=10*res0 && it>maxits/4 || isnan(res); error('!!! solver diverged, try again !!!'); end
            if max(abs(u(:)))>1e2 || max(abs(w(:)))>1e2, error('!!! solution is blowing up, try again !!!'); end
            if it==1; res0 = res; end
            fprintf(1,'    ---  it = %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,res,res/res0);
            if (nop>0) && ~mod(step,abs(nop)), semilogy(it,res,'r.','MarkerSize',10); axis tight; hold on; end
        end
        
    end  % iteration loop
    
    
    % update closures, constitutives, then plot and store model output
    if (nop~=0) && ~mod(step,abs(nop)); closures; constitutive; output; end
    
    % update time and step count
    time = time+dt;
    step = step+1;
    
    
    fprintf(1,'    solver time: %4.2f min \n',toc/60);
    
end  % time loop

if (mms), mms_results; end