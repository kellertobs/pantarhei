function [u, w, p] = multigrid_solve(Qvx, Qvz, Qf, f, PHS, NUM)

% get multigrid level
if size(f,3) > 1
    level = log2(min([size(f,2),size(f,3)]));
else
    level = log2(size(f,2));
end
TINY  = 1e-16;

u = 0.*Qvx;
w = 0.*Qvz;
p = 0.*Qf;

if  level == NUM.minlvl

    % solve coarsest level
    [u, w, p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, true);
    level = level + 1;

else

    % pre-smoothing
    [u, w, p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, false);

    % update closures and residuals
    closures; residuals;

    % restrict residuals to coarse grid (coarsening)
    res_p_coarse = restrict_p(res_p);
    res_w_coarse = restrict_w(res_w);
    res_u_coarse = restrict_u(res_u);
    f_coarse     = restrict_p(f    );

    NUM.Nz = NUM.Nz/2;
    NUM.Nx = NUM.Nx/(2*(NUM.Nx>=2)+1*(NUM.Nx<2));
    NUM.h  = NUM.h*2;
    NUM.alpha = NUM.alpha*0.95;
    NUM.beta  = NUM.beta*0.90;

    get_indices;

    % recursively solve update on coarser grid
    [upd_u_coarse, upd_w_coarse, upd_p_coarse] = multigrid_solve(res_u_coarse, res_w_coarse, res_p_coarse, f_coarse, PHS, NUM);

    % prolongate update to fine grid (interpolation)
    upd_p = prolong_p(upd_p_coarse, 0.*p, NUM);
    upd_w = prolong_w(upd_w_coarse, 0.*w, NUM);
    upd_u = prolong_u(upd_u_coarse, 0.*u, NUM);

    % remove mean from update to avoid solution drift
    upd_p = upd_p - mean(upd_p(:));
    if any(strcmp(NUM.BC, 'periodic'))
        upd_u = upd_u - mean(upd_u(:));
        upd_w = upd_w - mean(upd_w(:));
    end

    % apply correction
    u = u - NUM.gamma*upd_u;
    w = w - NUM.gamma*upd_w;
    p = p - NUM.gamma*upd_p;

    NUM.Nz = NUM.Nz*2;
    NUM.Nx = NUM.Nx*(2*(NUM.Nx>=2)+1*(NUM.Nx<2));
    NUM.h  = NUM.h/2;
    NUM.alpha = NUM.alpha/0.95;
    NUM.beta  = NUM.beta/0.90;

    get_indices;

    % post-smoothing
    [u, w, p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, false);

    if NUM.pltits
        figure(200);
        subplot(1,2,1);
        plot(p.'.*(-1).^level,1:2^(NUM.maxlvl-level):2^NUM.maxlvl); axis ij tight;
        subplot(1,2,2);
        plot(w.'.*(-1).^level,1:2^(NUM.maxlvl-level):2^NUM.maxlvl+1); axis ij tight;
        drawnow; pause(0.2)
    end
end
end

function coarse = restrict_p(fine)
    % Restriction (coarsening) operation on p, f fields
    if size(fine,3) > 1
        coarse = (fine(:, 1:2:end-1, 1:2:end-1) + fine(:, 2:2:end, 1:2:end-1) ...
               +  fine(:, 1:2:end-1, 2:2:end  ) + fine(:, 2:2:end, 2:2:end))/4;
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

function fine = prolong_p(coarse, fine, NUM)
    % Prolongation (interpolation) operation on p, f fields
    if NUM.Nx==1
        tmp = zeros(size(coarse,1),size(coarse,2)+2,1);
        tmp(:,2:end-1,1) = coarse;
        if strcmp(NUM.BC{1}, 'periodic')
            tmp(:,1  ,:) = tmp(:,end-1,:);
            tmp(:,end,:) = tmp(:,2    ,:);
        elseif strcmp(NUM.BC{1}, 'closed')
            tmp(:,1  ,:) = 2*tmp(:,2    ,:) - tmp(:,3    ,:);
            tmp(:,end,:) = 2*tmp(:,end-1,:) - tmp(:,end-2,:);
        else
            tmp(:,1  ,:) = tmp(:,2    ,:);
            tmp(:,end,:) = tmp(:,end-1,:);
        end
        fine(:, 1:2:end-1, 1) = (tmp(:, 1:end-2, 1)*1/4 + tmp(:, 2:end-1, 1)*3/4);
        fine(:, 2:2:end  , 1) = (tmp(:, 2:end-1, 1)*3/4 + tmp(:, 3:end-0, 1)*1/4);
    else
        tmp = zeros(size(coarse,1),size(coarse,2)+2,size(coarse,3)+2);
        tmp(:,2:end-1,2:end-1) = coarse;
        if strcmp(NUM.BC{1}, 'periodic')
            tmp(:,1  ,:) = tmp(:,end-1,:);
            tmp(:,end,:) = tmp(:,2    ,:);
        elseif strcmp(NUM.BC{1}, 'closed')
            tmp(:,1  ,:) = 2*tmp(:,2    ,:) - tmp(:,3    ,:);
            tmp(:,end,:) = 2*tmp(:,end-1,:) - tmp(:,end-2,:);
        else
            tmp(:,1  ,:) = tmp(:,2    ,:);
            tmp(:,end,:) = tmp(:,end-1,:);
        end
        if strcmp(NUM.BC{2}, 'periodic')
            tmp(:,:,1  ) = tmp(:,:,end-1);
            tmp(:,:,end) = tmp(:,:,2    );
        elseif strcmp(NUM.BC{2}, 'closed')
            tmp(:,:,1  ) = 2*tmp(:,:,2    ) - tmp(:,:,3    );
            tmp(:,:,end) = 2*tmp(:,:,end-1) - tmp(:,:,end-2);
        else
            tmp(:,:,1  ) = tmp(:,:,2    );
            tmp(:,:,end) = tmp(:,:,end-1);
        end
        fine(:, 1:2:end-1, 1:2:end-1) = (tmp(:, 1:end-2, 1:end-2)*1/16 + tmp(:, 2:end-1, 1:end-2)*3/16  ...
                                      +  tmp(:, 1:end-2, 2:end-1)*3/16 + tmp(:, 2:end-1, 2:end-1)*9/16);

        fine(:, 2:2:end  , 1:2:end-1) = (tmp(:, 2:end-1, 1:end-2)*3/16 + tmp(:, 3:end-0, 1:end-2)*1/16  ...
                                      +  tmp(:, 2:end-1, 2:end-1)*9/16 + tmp(:, 3:end-0, 2:end-1)*3/16);

        fine(:, 1:2:end-1, 2:2:end  ) = (tmp(:, 1:end-2, 2:end-1)*3/16 + tmp(:, 2:end-1, 2:end-1)*9/16  ...
                                      +  tmp(:, 1:end-2, 3:end-0)*1/16 + tmp(:, 2:end-1, 3:end-0)*3/16);

        fine(:, 2:2:end  , 2:2:end  ) = (tmp(:, 2:end-1, 2:end-1)*9/16 + tmp(:, 3:end-0, 2:end-1)*3/16  ...
                                      +  tmp(:, 2:end-1, 3:end-0)*3/16 + tmp(:, 3:end-0, 3:end-0)*1/16);
    end
end

function fine = prolong_w(coarse, fine, NUM)
    % Prolongation (interpolation) operation on w fields
    if NUM.Nx==1
        fine(:, 1:2:end  , :) = coarse;
        fine(:, 2:2:end-1, :) = coarse(:, 1:end-1, :)*1/2 + coarse(:, 2:end, :)*1/2;
    else
        tmp = zeros(size(coarse,1),size(coarse,2),size(coarse,3)+2);
        tmp(:,:,2:end-1) = coarse;
        if strcmp(NUM.BC{2}, 'periodic')
            tmp(:,:,1  ) = tmp(:,:,end-1);
            tmp(:,:,end) = tmp(:,:,2    );
        elseif strcmp(NUM.BC{2}, 'closed')
            tmp(:,:,1  ) = 2*tmp(:,:,2    ) - tmp(:,:,3    );
            tmp(:,:,end) = 2*tmp(:,:,end-1) - tmp(:,:,end-2);
        else
            tmp(:,:,1  ) = tmp(:,:,2    );
            tmp(:,:,end) = tmp(:,:,end-1);
        end
        fine(:, 1:2:end  , 1:2:end-1) = (tmp(:, :, 1:end-2)*1/4 + tmp(:, :, 2:end-1)*3/4);

        fine(:, 2:2:end-1, 1:2:end-1) = (tmp(:, 1:end-1, 1:end-2)*1/8 + tmp(:, 2:end, 1:end-2)*1/8  ...
                                      +  tmp(:, 1:end-1, 2:end-1)*3/8 + tmp(:, 2:end, 2:end-1)*3/8);

        fine(:, 1:2:end  , 2:2:end  ) = (tmp(:, :, 2:end-1)*3/4 + tmp(:, :, 3:end)*1/4);

        fine(:, 2:2:end-1, 2:2:end  ) = (tmp(:, 1:end-1, 2:end-1)*3/8 + tmp(:, 2:end, 2:end-1)*3/8  ...
                                      +  tmp(:, 1:end-1, 3:end  )*1/8 + tmp(:, 2:end, 3:end  )*1/8);
    end
end

function fine = prolong_u(coarse, fine, NUM)
    % Prolongation (interpolation) operation on p, f fields
    if NUM.Nx==1
        fine = 0.*fine;
    else
        tmp = zeros(size(coarse,1),size(coarse,2)+2,size(coarse,3));
        tmp(:,2:end-1,:) = coarse;
        if strcmp(NUM.BC{1}, 'periodic')
            tmp(:,1  ,:) = tmp(:,end-1,:);
            tmp(:,end,:) = tmp(:,2    ,:);
        elseif strcmp(NUM.BC{1}, 'closed')
            tmp(:,1  ,:) = 2*tmp(:,2    ,:) - tmp(:,3    ,:);
            tmp(:,end,:) = 2*tmp(:,end-1,:) - tmp(:,end-2,:);
        else
            tmp(:,1  ,:) = tmp(:,2    ,:);
            tmp(:,end,:) = tmp(:,end-1,:);
        end
        fine(:, 1:2:end-1, 1:2:end  ) = (tmp(:, 1:end-2, :)*1/4 + tmp(:, 2:end-1, :)*3/4);

        fine(:, 1:2:end-1, 2:2:end-1) = (tmp(:, 1:end-2, 1:end-1)*1/8 + tmp(:, 1:end-2, 2:end)*1/8  ...
                                      +  tmp(:, 2:end-1, 1:end-1)*3/8 + tmp(:, 2:end-1, 2:end)*3/8);

        fine(:, 2:2:end  , 1:2:end  ) = (tmp(:, 2:end-1, :)*3/4 + tmp(:, 3:end, :)*1/4);

        fine(:, 2:2:end  , 2:2:end-1) = (tmp(:, 2:end-1, 1:end-1)*3/8 + tmp(:, 2:end-1, 2:end)*3/8  ...
                                      +  tmp(:, 3:end  , 1:end-1)*1/8 + tmp(:, 3:end  , 2:end)*1/8);
    end
end