
% initialise simulation information and matrices
init;

% initialise time stepping loop
time = 0;
step = 0;
it   = 0;
dpi = 0.*p;
dwi = 0.*w;
dui = 0.*u;

history;

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
        
        it = it+1;  % update iteration count

        % store solution of previous two iterations
        dui = u-ui;  dwi = w-wi;  dpi = p-pi;
        uii = ui; ui = u;  wii = wi; wi = w;  pii = pi; pi = p;  fi = f;

        % update coefficient closures and phase advection on shear velocity
        % (only every 'nupd' iterations)
        if ~mod(it,nupd)
            closures; 
            fadv = advect(f, ushr, wshr, h, {advn, 'vdf'}, [2,3], BC);
        end
        
        % update constitutive relations
        constitutive;

        % update residual fields
        res_u =            +  Div_qvx        +  Gvx + Qvx     ;
        res_w =            +  Div_qvz        +  Gvz + Qvz     ;
        res_p =            +  Div_qf         +  Gf  + Gm./rho ;
        res_f = (f-fo)./dt + (fadv+fadvo)./2 - (Gf + Gfo)./2  ;

        % call manufactured solution (if benchmarking)
        if (mms); mms_calc_source; end
        
        % do not update phase fractions during initial step
        if step==0; res_f = 0.*res_f; end
        
        % prepare solution updates
        if any( strcmp(BC,'closed') )
            if strcmp(BC{1},'closed'), res_w (:,[1,end],:) = 0; end
            if strcmp(BC{2},'closed'), res_u (:,:,[1,end]) = 0; end
            upd_u = res_u.*dtau_u + sum(res_u).*dtau_ur;
            upd_w = res_w.*dtau_w + sum(res_w).*dtau_wr;
            upd_p = res_p.*dtau_p + sum(res_p).*dtau_pr;
            upd_f = res_f.*dtau_f;
            upd_f = upd_f - mean(upd_f(:));
        else
            upd_u = res_u.*dtau_u + sum(res_u).*dtau_ur;
            upd_w = res_w.*dtau_w + sum(res_w).*dtau_wr;
            upd_p = res_p.*dtau_p + sum(res_p).*dtau_pr;
            upd_f = res_f.*dtau_f;
            upd_u = upd_u - mean(upd_u(:));
            upd_w = upd_w - mean(upd_w(:));
            upd_p = upd_p - mean(upd_p(:));
            upd_f = upd_f - mean(upd_f(:));
        end

        u = (2*ui - 1/2*uii - alpha.*upd_u + beta.*dui)/(3/2) - diff(upd_u(:,:,[1 1:end end]),2,3)./8; upd_u(:,:,[1,end]) = 0;
        w = (2*wi - 1/2*wii - alpha.*upd_w + beta.*dwi)/(3/2) - diff(upd_w(:,[1 1:end end],:),2,2)./8; upd_w(:,[1,end],:) = 0;
        p = (2*pi - 1/2*pii - alpha.*upd_p + beta.*dpi)/(3/2) - diff(upd_p(:,[1 1:end end],:),2,2)./8 - diff(upd_p(:,:,[1 1:end end]),2,3)./8;

        ur = sum(omvx.*u,1);
        wr = sum(omvz.*w,1);
        pr = sum(omfc.*p,1);

        Du = u-ur;
        Dw = w-wr;
        Dp = p-pr;

        % update phase fractions (only every 'nupd' iterations)
        if ~mod(it,nupd); f = fi - alpha.*upd_f; f = max(flim,min(1-flim,f)); end
        
        % print iteration diagnostics
        if ~mod(it,nupd)
            % get residual norm
            resflds = [ norm(res_u.*dtau_u,'fro')./(norm(u,'fro')+TINY) ; 
                        norm(res_w.*dtau_w,'fro')./(norm(w,'fro')+TINY) ;
                        norm(res_p.*dtau_p,'fro')./(norm(p,'fro')+TINY) ;
                        norm(sum(res_u).*dtau_ur,'fro')./(norm(ur,'fro')+TINY) ; 
                        norm(sum(res_w).*dtau_wr,'fro')./(norm(wr,'fro')+TINY) ;
                        norm(sum(res_p).*dtau_pr,'fro')./(norm(pr,'fro')+TINY)];
            res = sum(resflds(:));
            if (res>=100*res0 && it>maxits/4) || isnan(res) || (step>0 && res>1); error('!!! solver diverged, try again !!!'); end
            if max(abs(u(:)))>1e2 || max(abs(w(:)))>1e2, error('!!! solution is blowing up, try again !!!'); end
            if it==nupd || res>res0; res0 = res; end
            fprintf(1,'    ---  it = %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,res,res/res0);
        end
        if (nop>0) && ~mod(step,abs(nop)) && ~mod(it,5*nupd)
            figure(f10);
            semilogy(it,resflds(2,:),'r.','MarkerSize',10); axis tight; hold on;
            semilogy(it,resflds(3,:),'b.','MarkerSize',10);
            semilogy(it,resflds(5,:),'m.','MarkerSize',10);
            semilogy(it,resflds(6,:),'c.','MarkerSize',10);
            semilogy(it,res,'k.','MarkerSize',10);
            drawnow;
            figure(20);
            subplot(1,2,1)
            plot(res_w.*dtau_w./norm(w(:)).*sqrt(length(w(:))),1:Nz+1); axis tight; hold on;
            plot(sum(res_w).*dtau_wr./norm(wr(:)).*sqrt(length(wr(:))),1:Nz+1); hold off;
            subplot(1,2,2)
            plot(res_p.*dtau_p./norm(p(:)).*sqrt(length(p(:))),1:Nz+0); axis tight; hold on;
            plot(sum(res_p).*dtau_pr./norm(pr(:)).*sqrt(length(pr(:))),1:Nz+0); hold off; drawnow;
        end
    end  % iteration loop
    
    % populate history matrices
    history;

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

