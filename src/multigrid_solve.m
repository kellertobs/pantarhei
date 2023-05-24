function [u, w, p] = multigrid_solve(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM)
% Multigrid parameters
if size(f,3) > 1
    max_levels = log2(min([size(f,2),size(f,3)]));
else
    max_levels = log2(size(f,2));
end
level = max_levels;
TINY  = 1e-16;

% Coarsest grid solve
if  level == NUM.minlvl
    % solve coarsest level
    [u, w, p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, true);
    level = level + 1;
else
    % Pre-smoothing
    [u, w, p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, false);

    closures;
    constitutive;

    % update residual fields
    res_u = Div_qvx + Gvx - Qvx;
    res_w = Div_qvz + Gvz - Qvz;
    res_p = Div_qf  + Gf  - Qf ;

    if strcmp(NUM.BC{1},'closed'), res_w(:,[1,end],:) = 0; end
    if strcmp(NUM.BC{2},'closed'), res_u(:,:,[1,end]) = 0; end

    % Restriction (coarsening)
    res_p_coarse = restrict_p(res_p);
    res_w_coarse = restrict_w(res_w);
    res_u_coarse = restrict_u(res_u);
    f_coarse     = restrict_p(f    );

    NUM.Nz = size(res_p_coarse,2); 
    if strcmp(NUM.BC{1},'periodic')
        NUM.icz = [NUM.Nz,1:NUM.Nz,1]; NUM.imz = [NUM.Nz,1:NUM.Nz]; NUM.ipz = [1:NUM.Nz,1];
    elseif strcmp(NUM.BC{1},'open') || strcmp(NUM.BC{1},'closed')
        NUM.icz = [1,1:NUM.Nz,NUM.Nz]; NUM.imz = [1,1:NUM.Nz]; NUM.ipz = [1:NUM.Nz,NUM.Nz];
    end

    NUM.Nx = size(res_p_coarse,3);
    % now x direction
    if strcmp(NUM.BC{2},'periodic')
        NUM.icx = [NUM.Nx,1:NUM.Nx,1]; NUM.imx = [NUM.Nx,1:NUM.Nx]; NUM.ipx = [1:NUM.Nx,1];
    elseif strcmp(NUM.BC{2},'open') || strcmp(NUM.BC{2},'closed')
        NUM.icx = [1,1:NUM.Nx,NUM.Nx]; NUM.imx = [1,1:NUM.Nx]; NUM.ipx = [1:NUM.Nx,NUM.Nx];
    end

    NUM.h     = NUM.h*2;
    NUM.alpha = NUM.alpha*0.8;
    NUM.smth  = NUM.smth*1.2;

    % Recursively solve on coarser grid
    [upd_u_coarse, upd_w_coarse, upd_p_coarse] = multigrid_solve(0.*res_u_coarse, 0.*res_w_coarse, 0.*res_p_coarse, res_u_coarse, res_w_coarse, res_p_coarse, f_coarse, PHS, NUM);

    % Prolongation (interpolation)
    upd_p = prolong_p(upd_p_coarse, p, NUM.icz, NUM.icx);
    upd_w = prolong_w(upd_w_coarse, w, NUM.icz, NUM.icx);
    upd_u = prolong_u(upd_u_coarse, u, NUM.icz, NUM.icx);

    upd_p = upd_p - mean(upd_p(:));
    if any(strcmp(NUM.BC, 'periodic'))
        upd_u = upd_u - mean(upd_u(:));
        upd_w = upd_w - mean(upd_w(:));
    end

    % apply correction
    u = u - upd_u;
    w = w - upd_w;
    p = p - upd_p;

    NUM.Nz = size(p,2);
    if strcmp(NUM.BC{1},'periodic')
        NUM.icz = [NUM.Nz,1:NUM.Nz,1]; NUM.imz = [NUM.Nz,1:NUM.Nz]; NUM.ipz = [1:NUM.Nz,1];
    elseif strcmp(NUM.BC{1},'open') || strcmp(NUM.BC{1},'closed')
        NUM.icz = [1,1:NUM.Nz,NUM.Nz]; NUM.imz = [1,1:NUM.Nz]; NUM.ipz = [1:NUM.Nz,NUM.Nz];
    end

    NUM.Nx = size(p,3);
    % now x direction
    if strcmp(NUM.BC{2},'periodic')
        NUM.icx = [NUM.Nx,1:NUM.Nx,1]; NUM.imx = [NUM.Nx,1:NUM.Nx]; NUM.ipx = [1:NUM.Nx,1];
    elseif strcmp(NUM.BC{2},'open') || strcmp(NUM.BC{2},'closed')
        NUM.icx = [1,1:NUM.Nx,NUM.Nx]; NUM.imx = [1,1:NUM.Nx]; NUM.ipx = [1:NUM.Nx,NUM.Nx];
    end

    NUM.h     = NUM.h/2;
    NUM.alpha = NUM.alpha/0.8;
    NUM.smth  = NUM.smth/1.2;

    % Post-smoothing
    [u, w, p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, false);
end
end

function [u, w, p, upd_u, upd_w, upd_p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, minlvl)

if minlvl
    tol    = NUM.minlvl_tol;
    maxits = NUM.minlvl_smth;
else
    tol    = 1e-16;
    maxits = NUM.smth;
end

alpha  = NUM.alpha;
beta   = NUM.beta;

ui = u; 
wi = w; 
pi = p;

% update coefficient closures and phase advection on shear velocity
closures;

resnorm = 1; it = 1;
while resnorm >= tol && it < maxits

    % store solution of previous two iterations
    uii = ui;  ui = u;
    wii = wi;  wi = w;
    pii = pi;  pi = p;

    % update constitutive relations
    constitutive;

    % update residual fields
    res_u = Div_qvx + Gvx - Qvx;
    res_w = Div_qvz + Gvz - Qvz;
    res_p = Div_qf  + Gf  - Qf ;

    if strcmp(NUM.BC{1},'closed'), res_w(:,[1,end],:) = 0; end
    if strcmp(NUM.BC{2},'closed'), res_u(:,:,[1,end]) = 0; end

    upd_u = res_u.*dtau_u;
    upd_w = res_w.*dtau_w;
    upd_p = res_p.*dtau_p;

    upd_p = upd_p - mean(upd_p(:));
    if any(strcmp(NUM.BC, 'periodic'))
        upd_u = upd_u - mean(upd_u(:));
        upd_w = upd_w - mean(upd_w(:));
    end

    u = u - alpha.*upd_u + beta.*(u-uii);
    w = w - alpha.*upd_w + beta.*(w-wii);
    p = p - alpha.*upd_p + beta.*(p-pii);

    if ~mod(it-1,10)
        resnorm = norm(upd_w,'fro')./norm(w,'fro') + norm(upd_p,'fro')./(norm(p,'fro')+1e-16);
    end

    it = it+1;
end
end

function coarse = restrict_p(fine)
    % Restriction (coarsening) operation on p, f fields
    if size(fine,3) > 1
        coarse = (fine(:, 1:2:end-1, 1:2:end  ) + fine(:, 2:2:end, 1:2:end) ...
               +  fine(:, 1:2:end  , 1:2:end-1) + fine(:, 1:2:end, 2:2:end))/4;
    else
        coarse = (fine(:, 1:2:end-1, 1) + fine(:, 2:2:end, 1))/2;
    end
end

function coarse = restrict_w(fine)
    % Restriction (coarsening) operation on w fields
    if size(fine,3) > 1
        coarse = (fine(:, 1:2:end, 1:2:end-1) + fine(:, 1:2:end, 2:2:end))/2;
    else
        coarse = fine(:, 1:2:end, 1);
    end
end

function coarse = restrict_u(fine)
    % Restriction (coarsening) operation on u fields
    if size(fine,3) > 2
        coarse = (fine(:, 1:2:end-1, 1:2:end) + fine(:, 2:2:end, 1:2:end))/2;
    else
        coarse = (fine(:, 1:2:end-1, :) + fine(:, 2:2:end, :))/2;
    end
end

function fine = prolong_p(coarse, fine, icz, icx)
    % Prolongation (interpolation) operation on p, f fields
    fine(:, 1:2:end-1, 1) = (coarse(:, icz(1:end-2), 1)*1/4 + coarse(:, icz(2:end-1), 1)*3/4);
    fine(:, 2:2:end  , 1) = (coarse(:, icz(2:end-1), 1)*3/4 + coarse(:, icz(3:end-0), 1)*1/4);
end

function fine = prolong_w(coarse, fine, icz, icx)
    % Prolongation (interpolation) operation on w fields
    fine(:, 1:2:end  , 1:2:end) = coarse;
    fine(:, 2:2:end-1, 1:2:end) = (coarse(:, 1:end-1, :) + coarse(:, 2:end, :))/2;
end

function fine = prolong_u(coarse, fine, icz, icx)
    % Prolongation (interpolation) operation on p, f fields
    fine(:, 1:2:end-1, :) = (coarse(:, icz(1:end-2), :)*1/4 + coarse(:, icz(2:end-1), :)*3/4);
    fine(:, 2:2:end  , :) = (coarse(:, icz(2:end-1), :)*3/4 + coarse(:, icz(3:end-0), :)*1/4);
end