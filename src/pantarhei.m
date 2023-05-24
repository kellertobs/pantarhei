
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
    
    % store phase fractions of previous step
    fo = f;  dto = dt; dfdto = dfdt;
    
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
        
        closures; constitutive;
        if step>0
            % update physical time step [s]
            dt   = min([ 2*dto;
                         NUM.cfl/(max(abs([qfx(:);qfz(:)]))/(NUM.h/2) + max(abs(Gf(:)))./1e-2);
                         NUM.cfl*(NUM.h/2)./max(abs([ushr(:);wshr(:)]) + NUM.TINY) ]);

            fadv = -advect(f, ushr, wshr, NUM.h, {NUM.advn, 'vdf'}, [2,3], NUM.BC);
            dfdt = fadv + Gf;
            f    = fo + (dfdt+dfdto)/2*dt;
        end

        % source fields
        rhomix = mean(mean(sum(f.*rho,1)));
        Qvx    = fx.*(rhox-rhomix).*PHS.grav(2);
        Qvz    = fz.*(rhoz-rhomix).*PHS.grav(1);
        Qf     = zeros(size(f));

        ui = u; wi = w; pi = p;
        [u, w, p] = multigrid_solve(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM);
        upd_u = u-ui; upd_w = w-wi;  upd_p = p-pi;

        closures; constitutive;
        res_u = Div_qvx + Gvx - Qvx;
        res_w = Div_qvz + Gvz - Qvz;
        res_p = Div_qf  + Gf  - Qf ;

        it = it+1;  % update iteration count

        % get residual norm
        resflds = [ (norm(upd_u,'fro'))./(norm(u,'fro')+NUM.TINY);
                    (norm(upd_w,'fro'))./(norm(w,'fro')+NUM.TINY);
                    (norm(upd_p,'fro'))./(norm(p,'fro')+NUM.TINY)  ];

        res = sum(resflds);
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

            subplot(143); plot(upd_p(:,:,1)./norm(p,'fro')*sqrt(length(p(:))), NUM.z , sum(-res_p,1)./norm(p./dtau_p,'fro')*(length(p(:))), NUM.z , 'k-'); grid on; axis tight; title('upd p')
            subplot(144); plot(upd_w(:,:,1)./norm(w,'fro')*sqrt(length(w(:))), NUM.zw, sum(-res_w,1)./norm(w./dtau_w,'fro')*(length(w(:))), NUM.zw, 'k-'); grid on; axis tight; title('upd w')
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

