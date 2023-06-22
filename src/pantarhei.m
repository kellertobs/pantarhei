
% initialise simulation information and matrices
init;

% initialise time stepping loop
time = 0;
step = 0;
it   = 0;

history;

% restart run?
if IO.restart>0
    % load previous file to get ui, wi, pi, fi
    % NB: these are not exactly the same as (ui, wi, pi, fi) if nop~=1, 
    % but closer than zeros initialisation matrices
    fname = [IO.outdir IO.RunID,'/',IO.RunID,'_',num2str(IO.restart-1),'.mat'];
    load(fname); fprintf(1,'Loaded %s to get [u,w,p,f] of previous time step.\n', fname);
    ui = u; wi = w; pi = p;  fi = f;
    
    % load current file 
    fname = [IO.outdir IO.RunID,'/',IO.RunID,'_',num2str(restart),'.mat'];
    load(fname); fprintf(1,'Loaded %s.\n', fname);
    rsstep = IO.restart*nop;  % current step
    step   = rsstep + 1;
end

% initialise coefficient closures and constitutive relations
closures; constitutive;

while time <= NUM.tend && step <= NUM.NtMax  % keep stepping until final run time reached
    
    tic;
    
    % print time step diagnostic
    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    
    if     strcmp(NUM.tint,'be1im') || step==1 % first step / 1st-order backward-Euler implicit scheme
        a1 = 1; a2 = 1; a3 = 0;
        b1 = 1; b2 = 0; b3 = 0;
    elseif strcmp(NUM.tint,'bd2im') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 1;   b2 =  0;  b3 = 0;
    elseif strcmp(NUM.tint,'cn2si')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        a1 = 1;   a2 = 1;   a3 = 0;
        b1 = 1/2; b2 = 1/2; b3 = 0;
    elseif strcmp(NUM.tint,'bd2si')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 3/4; b2 = 2/4; b3 = -1/4;
    end

    % store phase fractions of previous step
    foo = fo; fo = f;  
    dfdtoo = dfdto; dfdto = dfdt;  
    dto = dt; 
    
    % initialise non-linear iteration loop
    res  = 1e3;
    res0 = 1;
    it   = 0;
    itvec = nan(ceil(NUM.maxits), 4); % collect iteration information

    if step == 0 % require more stringent convergence criterion for first step
        conv_crit = @(resc,res0c,itc) (resc >= NUM.atol && itc <= 10*NUM.maxits || itc < NUM.minits);
    else
        conv_crit = @(resc,res0c,itc) (resc >= NUM.atol && resc/res0c >= NUM.rtol && itc <= NUM.maxits || itc < NUM.minits);
    end
    
    % make iterative convergence plot
    if IO.pltits, f10 = figure(10); clf; set(f10,'Position',[600,550,1200,500]); rc = lines(4); end
    
    while  conv_crit(res,res0,it) > 0 % keep stepping until convergence criterion reached
        
        if step>0
            % update physical time step [s]
            dt   = min([ 2*dto;
                         NUM.cfl/(max(abs([qfx(:);qfz(:)]))/(NUM.h/2) + max(abs(Gf(:)))./1e-2);
                         NUM.cfl*(NUM.h/2)./max(abs([ushr(:);wshr(:)]) + NUM.TINY) ]);

            fadv = -advect(f, ushr, wshr, NUM.h, {NUM.advn, 'vdf'}, [2,3], NUM.BC);
            dfdt = fadv + Gf;
            f    = (a2*fo+a3*foo + (b1*dfdt + b2*dfdto + b3*dfdtoo)*dt)/a1;
            f    = max(NUM.flim,min(1-NUM.flim, f ));
        end

        % update source fields
        rhomix = mean(mean(sum(f.*rho,1)));
        Qvx    = fx.*(rhox-rhomix).*PHS.grav(2);
        Qvz    = fz.*(rhoz-rhomix).*PHS.grav(1);
        Qf     = zeros(size(f));

        % update closures and residuals
        closures; residuals;

        if NUM.direct
            [upd_u, upd_w, upd_p] = direct_solve(u, w, p, f, Qvz, Qf, PHS, NUM);
        else
            [upd_u, upd_w, upd_p] = multigrid_solve(res_u, res_w, res_p, f, PHS, NUM);
        end

        % % apply correction
        u = u - upd_u;
        w = w - upd_w;
        p = p - upd_p;

        % update closures and residuals
        closures; residuals;

        it = it+1;  % update iteration count

        % get residual norms
        resflds = [ (norm(upd_u,'fro'))./(norm(u,'fro')+NUM.TINY);
                    (norm(upd_w,'fro'))./(norm(w,'fro')+NUM.TINY);
                    (norm(upd_p,'fro'))./(norm(p,'fro')+NUM.TINY)  ];

        res = sum(resflds)/(NUM.ndim+1);
        if it==1 || res>res0; res0 = res; end

        itvec(it+1,:) = [it, resflds(:)'];

        fprintf(1,'    ---  it = %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,res,res/res0);

        if (step>0 && res>=100*res0 && it>maxits/4) || isnan(res) || (step>0 && res>1); error('!!! solver diverged, try again !!!'); end
        if max(abs(u(:)))>1e2 || max(abs(w(:)))>1e2, error('!!! solution is blowing up, try again !!!'); end

        if IO.pltits && ~mod(it,1)
            % plot the residuals as function of the iterations
            figure(f10);
            subplot(1,4,[1,2]);
            for jvar=1:3, semilogy(it,resflds(jvar),'.','MarkerSize',10,'Color',rc(jvar  ,:),'linewidth',1); hold on; end
            semilogy(it, res,'ko','MarkerSize',8,'linewidth',1); axis tight;

            subplot(143); plot(upd_p(:,:,1)./norm(p,'fro')*sqrt(length(p(:))), NUM.z ); grid on; axis tight; title('upd p')
            subplot(144); plot(upd_w(:,:,1)./norm(w,'fro')*sqrt(length(w(:))), NUM.zw); grid on; axis tight; title('upd w')
            drawnow;
        end
     end  % iteration loop
    
    % populate history matrices
    history;

    % update closures, constitutives, then plot and store model output
    if (IO.nop~=0) && ~mod(step,abs(IO.nop))
        closures; constitutive; output; 
        if (mms), mms_results; end
    end

    % update time and step count
    if step>0; time = time+dt; end
    step = step+1;
    
    
    fprintf(1,'    solver time: %4.2f min \n',toc/60);
    
end  % time loop

