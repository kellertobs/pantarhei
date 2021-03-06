
% print header
fprintf(1,'\n\n****************************************************\n');
fprintf(1,'*****  pantarhei | multi-phase flow simulator  *****\n');
fprintf(1,'****************************************************\n\n');

% prepare workspace
if svop && ~exist(['../out/',RunID],'dir'); mkdir(['../out/',RunID]); end
if svop, save(['../out/',RunID,'/',RunID,'_par.mat']); end

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
x     = -D/2+h/2:h:D/2-h/2;
[X,~] = meshgrid(x,x);
X     = repmat(X,1,1,NPHS);
Z     = permute(X,[3,2,1]);
X     = permute(X,[3,1,2]);

% initialise indexing for boundary condition stencils
if     strcmp(BC,'periodic'); ic = [N,1:N,1]; im = [N,1:N]; ip = [1:N,1];
elseif strcmp(BC,'open') || strcmp(BC,'closed'); ic = [1,1:N,N]; im = [1,1:N]; ip = [1:N,N]; end

% initialise smoothed random and gaussian perturbation fields
rng(15);
rnd = randn(NPHS,N,N);
for i = 1:smth
    rnd = rnd + diff(rnd(:,ic,:),2,2)./8 + diff(rnd(:,:,ic),2,3)./8;
end
rnd = rnd./max(abs(rnd(:)));
gsn = exp(-X.^2./(D/5).^2).*exp(-Z.^2./(D/5).^2);

% intialise solution, auxiliary, and residual fields
u      = zeros(NPHS,N  ,N+1);  ui = u;  ustar = mean(u,1);  usegr = 0*u;  res_u = 0*u;  dtau_u = res_u;
w      = zeros(NPHS,N+1,N  );  wi = w;  wstar = mean(w,1);  wsegr = 0*w;  res_w = 0*w;  dtau_w = res_w;
p      = zeros(NPHS,N  ,N  );  pi = p;  pstar = mean(p,1);  pcmpt = 0*p;  res_p = 0*p;  dtau_p = res_p;

if exist('fInit','var'),    addpath('./inits/'); f = fInit(X, Z);
else,                       f = f0 + dfr.*rnd + dfg.*gsn;
end
f = max(1e-16,min(1-1e-16,f));  f = f./sum(f,1);  fo = f;  fi = f;  res_f = 0*f;  dtau_f = res_f;
if (mms), mms_init_phasefrac; end

qvxx   = zeros(NPHS,N,N  );  qvzz = zeros(NPHS,N,N  );  qvxz = zeros(NPHS,N+1,N+1);
qfx    = zeros(NPHS,N,N+1);  qfz  = zeros(NPHS,N+1,N);
Gvx    = zeros(NPHS,N,N+1);  Gvz  = zeros(NPHS,N+1,N);
Gf     = zeros(NPHS,N,N  );  Gfo  = Gf;
delta  = zeros(NPHS,NPHS,N,N);
rho    = rho0.*ones(size(f));
rhomix = mean(mean(sum(f.*rho,1)));

% initialise coefficient closures and constitutive relations
closures; constitutive;

% initialise time stepping loop
time = 0;
step = 0;
while time <= tend && step <= NtMax  % keep stepping until final run time reached

    tic;
    
    % print time step diagnostics
    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    
    % store phase fractions of previous step
    fo = f;  dto = dt; Gfo = Gf;
    
    % initialise non-linear iteration loop
    res  = 1e3;
    res0 = res;
    it   = 0;
    while res >= atol && res/res0 >= rtol && it <= maxits || it < minits  % keep stepping until residual norm below tolerance
                
        % store solution of previous two iterations
        uii = ui;  ui = u;  wii = wi;  wi = w;  pii = pi;  pi = p;  fii = fi;  fi = f;

        it = it+1;  % update iteration count

        % update coefficient closures (only every 'nupd' iterations)
        if ~mod(it,nupd); closures; end
        
        % update constitutive relations
        constitutive;
        
        % update physical time step
        dt   = min(2*dto,cfl/(max(abs([qfx(:);qfz(:)]))/(h/2) + max(abs(Gf(:)))./1e-2));  % [s]
        
        % update residual fields
        res_u =             + diff(qvxx(:,:,ic),1,3)./h + diff(qvxz,1,2)./h - Gvx - Qvx    ;
        res_w =             + diff(qvzz(:,ic,:),1,2)./h + diff(qvxz,1,3)./h - Gvz - Qvz    ;
        res_p =             + diff(qfx         ,1,3)./h + diff(qfz ,1,2)./h - Gf  - Gm./rho;
        res_f = (f-fo)./dt                                                  +(Gf + Gfo)./2 ;
        
        % call manufactured solution (if benchmarking)
        if (mms); mms_calc_source; end
        
        % do not update phase fractions during initial step
        if step==0; res_f = 0.*res_f; end
        
        % prepare solution updates
        if strcmp(BC,'closed')
            res_u(:,:,[1,end]) = 0; 
            res_w(:,[1,end],:) = 0;
            upd_u = res_u.*dtau_u;
            upd_w = res_w.*dtau_w;
            upd_p = res_p.*dtau_p;
            upd_f = res_f.*dtau_f;
        else
            upd_u = res_u.*dtau_u-mean(res_u(:).*dtau_u(:));
            upd_w = res_w.*dtau_w-mean(res_w(:).*dtau_w(:));
            upd_p = res_p.*dtau_p-mean(res_p(:).*dtau_p(:));
            upd_f = res_f.*dtau_f-mean(res_f(:).*dtau_f(:));
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
            res = norm(res_u(:).*dtau_u(:),2)./(norm(u(:),2)+1e-32) ...
                + norm(res_w(:).*dtau_w(:),2)./(norm(w(:),2)+1e-32) ...
                + norm(res_p(:).*dtau_p(:),2)./(norm(p(:),2)+1e-32) ...
                + norm(res_f(:).*dtau_f(:),2)./(norm(f(:),2)+1e-32);
            if res>=10*res0 && it>maxits/4 || isnan(res); error('!!! solver diverged, try again !!!'); end
            if max(abs(u(:)))>1e2 || max(abs(w(:)))>1e2, error('!!! solution is blowing up, try again !!!'); end
            if it==1; res0 = res; end
            fprintf(1,'    ---  it = %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,res,res/res0);
            if (nop); f10 = figure(10); if it==1; clf; end; semilogy(it,res,'r.','MarkerSize',10); axis tight; box on; hold on; drawnow; end
        end
        
    end  % iteration loop
        
    % update closures, constitutives, then plot and store model output
    if (nop) && ~mod(step,nop); closures; constitutive; output; end
    
    % update time and step count
    time = time+dt;
    step = step+1;
    
    fprintf(1,'    solver time: %4.2f min \n',toc/60);
    
end  % time loop

if (mms), mms_results; end