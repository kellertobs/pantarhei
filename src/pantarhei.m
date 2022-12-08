
% initialise simulation information and matrices
init;

% initialise time stepping loop
time = 0;
step = 0;
it   = 0;

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
    if (nop>0) && ~mod(step,abs(nop)), f10 = figure(10);  f10.Visible = figvis; end
    
    while  conv_crit(res,res0,it) > 0 % keep stepping until convergence criterion reached
        
        % store solution of previous two iterations
        uii = ui;  ui = u;  wii = wi;  wi = w;  pii = pi;  pi = p;  fii = fi;  fi = f;
        
        it = it+1;  % update iteration count
        
        % update coefficient closures and phase advection on shear velocity
        % (only every 'nupd' iterations)
        if ~mod(it,nupd)
            closures; 
            fadv = advect(f, ushr, wshr, h, {advn, 'vdf'}, [2,3], BC);
        end
        
        % update constitutive relations
        constitutive;        
           
        % update physical time step [s]
        dt   = min([ 2*dto; 
                     cfl/(max(abs([qfx(:);qfz(:)]))/(h/2) + max(abs(Gf(:)))./1e-2); 
                     cfl*0.5*h./max(abs([ushr(:);wshr(:)]) + TINY) ]);  

        % update residual fields
        res_u =            +    Div_qvx      + Gvx + Qvx     ;
        res_w =            +    Div_qvz      + Gvz + Qvz     ;
        res_p =            +    Div_qf       + Gf  + Gm./rho ;
        res_f = (f-fo)./dt + (fadv+fadvo)./2 -(Gf + Gfo)./2 ;
        
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
            resflds = [ norm(res_u(:).*dtau_u(:),2)./(norm(u(:),2)+TINY) ; 
                        norm(res_w(:).*dtau_w(:),2)./(norm(w(:),2)+TINY) ;
                        norm(res_p(:).*dtau_p(:),2)./(norm(p(:),2)+TINY) ;
                        norm(res_f(:).*dtau_f(:),2)./(norm(f(:),2)+TINY)];
            res = sum(resflds);
            if (res>=10*res0 && it>maxits/4) || isnan(res) || (step>0 && res>1e-2); error('!!! solver diverged, try again !!!'); end
            if max(abs(u(:)))>1e2 || max(abs(w(:)))>1e2, error('!!! solution is blowing up, try again !!!'); end
            if it==1; res0 = res; end
            fprintf(1,'    ---  it = %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,res,res/res0);
            if (nop>0) && ~mod(step,abs(nop)), semilogy(it,res,'r.','MarkerSize',10); axis tight; hold on; drawnow; end
        end
        
    end  % iteration loop
    
    
    % update closures, constitutives, then plot and store model output
    if (nop~=0) && ~mod(step,abs(nop))
        closures; constitutive; output; 
        if (mms), mms_results; end
    end
    
    % update time and step count
    time = time+dt;
    step = step+1;
    
    
    fprintf(1,'    solver time: %4.2f min \n',toc/60);
    
end  % time loop

