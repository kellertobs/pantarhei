
% initialise simulation information and matrices
init;

% initialise time stepping loop
time = 0;
step = 0;
it   = 0;

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
    itvec = nan(ceil(maxits/nupd), 8); % collect iteration information

    if step == 0 % require more stringent convergence criterion for first step
        conv_crit = @(resc,res0c,itc) (resc >= atol && itc <= 10*maxits || itc < minits);
    else
        conv_crit = @(resc,res0c,itc) (resc >= atol && resc/res0c >= rtol && itc <= maxits || itc < minits);
    end
    
    % make iterative convergence plot
    if pltits, f10 = figure(10); clf; set(f10,'Position',[600,550,1200,500]); rc = lines(4); end
    
    while  conv_crit(res,res0,it) > 0 % keep stepping until convergence criterion reached
        
        % store solution of previous two iterations
        uii = ui;  ui = u;  upd_ui = upd_u;
        wii = wi;  wi = w;  upd_wi = upd_w;
        pii = pi;  pi = p;  upd_pi = upd_p;
        fii = fi;  fi = f;  upd_fi = upd_f;
        
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
        res_f = (f-fo)./dt + (fadv+fadvo)./2 -(Gf + Gfo)./2  ;
        
        % call manufactured solution (if benchmarking)
        if (mms); mms_calc_source; end
        
        % do not update phase fractions during initial step
        if step==0; res_f = 0.*res_f; end

        if strcmp(BC{1},'closed'), res_w(:,[1,end],:) = 0; end
        if strcmp(BC{2},'closed'), res_u(:,:,[1,end]) = 0; end

        upd_u = res_u.*dtau_u +  sum(res_u,1).*dtau_ustar;
        upd_w = res_w.*dtau_w +  sum(res_w,1).*dtau_wstar;
        upd_p = res_p.*dtau_p +  sum(res_p,1).*dtau_pstar; 
        upd_f = res_f.*dtau_f - mean(res_f(:).*dtau_f(:));
    
        if any(strcmp(BC, 'periodic'))
            upd_u = upd_u - mean(upd_u(:));
            upd_w = upd_w - mean(upd_w(:));
            upd_p = upd_p - mean(upd_p(:));          
        end

        u = ui - alpha.*upd_u + beta.*(ui-uii);
        w = wi - alpha.*upd_w + beta.*(wi-wii);
        p = pi - alpha.*upd_p + beta.*(pi-pii);

        ustar     = sum(omvx.*u,1);
        wstar     = sum(omvz.*w,1);
        pstar     = sum(omfc.*p,1);

        % update phase fractions (only every 'nupd' iterations)
        if ~mod(it,nupd); f = fi - alpha.*upd_f; f = max(flim,min(1-flim,f)); end
     
        % print iteration diagnostics
        if ~mod(it,nupd) || it==1
            % get residual norm
            resflds = [ (norm(    res_u   ,'fro'))./(norm(u    ./dtau_u    ,'fro')+TINY);
                        (norm(    res_w   ,'fro'))./(norm(w    ./dtau_w    ,'fro')+TINY);
                        (norm(    res_p   ,'fro'))./(norm(p    ./dtau_p    ,'fro')+TINY);
                        (norm(    res_f   ,'fro'))./(norm(f    ./dtau_f    ,'fro')+TINY); 
                        (norm(sum(res_u,1),'fro'))./(norm(ustar./dtau_ustar,'fro')+TINY);
                        (norm(sum(res_w,1),'fro'))./(norm(wstar./dtau_wstar,'fro')+TINY);
                        (norm(sum(res_p,1),'fro'))./(norm(pstar./dtau_pstar,'fro')+TINY)];
            
            res = sum(resflds);
            
            itvec(round(it/nupd)+1,:) = [it, resflds(:)'];

            fprintf(1,'    ---  it = %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,res,res/res0);

            if (step>0 && res>=100*res0 && it>maxits/4) || isnan(res) || (step>0 && res>1); error('!!! solver diverged, try again !!!'); end
            if max(abs(u(:)))>1e2 || max(abs(w(:)))>1e2, error('!!! solution is blowing up, try again !!!'); end
            if it==1; res0 = res; end

            if pltits && it>1
                % plot the residuals as function of the iterations
                figure(f10);
                subplot(1,4,[1,2]);
                for jvar=1:4, semilogy(it,resflds(jvar),'.','MarkerSize',10,'Color',rc(jvar  ,:),'linewidth',1); hold on; end
                for jvar=5:7, semilogy(it,resflds(jvar),'*','MarkerSize', 6,'Color',rc(jvar-4,:),'linewidth',0.2); hold on; end
                semilogy(it, res,'ko','MarkerSize',8,'linewidth',1); axis tight; drawnow;

                if ~mod(it,nupd*2)
                    subplot(143); plot(res_p(:,:,1).*dtau_p(:,:,1), z , sum(res_p(:,:,1),1).*dtau_pstar(:,:,1), z , 'k-'); grid on; axis tight; title('upd p')
                    subplot(144); plot(res_w(:,:,1).*dtau_w(:,:,1), zw, sum(res_w(:,:,1),1).*dtau_wstar(:,:,1), zw, 'k-'); grid on; axis tight; title('upd w')

                    if ~mod(it,nupd*20)
                        subplot(143); xlimits1 = xlim; subplot(144); xlimits2 = xlim;
                    elseif it>20*nupd
                        subplot(143); xlim(xlimits1); subplot(144); xlim(xlimits2);
                    end
                    drawnow;
                end
            end
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

