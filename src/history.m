
% collect history of useful values

% if we are in the zeroth step, initialize the variables to allocate memory
if step==0
    hist = struct('t'  , nan(NUM.NtMax,1), ...
                  'u'  , nan(NUM.NtMax,NUM.NPHS,3), ...
                  'w'  , nan(NUM.NtMax,NUM.NPHS,3), ...
                  'p'  , nan(NUM.NtMax,NUM.NPHS,3), ...
                  'f'  , nan(NUM.NtMax,NUM.NPHS,3), ...
                  'res', nan(NUM.NtMax,NUM.NPHS,4), ...
                  'restot',nan(NUM.NtMax,1)    );
end

%%

hist.t(step+1) = time;

% solution fields
hist.u(step+1,:,1) =  min(u,[],[2,3]);
hist.u(step+1,:,2) = mean(u,   [2,3]);
hist.u(step+1,:,3) =  max(u,[],[2,3]);

hist.w(step+1,:,1) =  min(w,[],[2,3]);
hist.w(step+1,:,2) = mean(w,   [2,3]);
hist.w(step+1,:,3) =  max(w,[],[2,3]);

hist.p(step+1,:,1) =  min(p,[],[2,3]);
hist.p(step+1,:,2) = mean(p,   [2,3]);
hist.p(step+1,:,3) =  max(p,[],[2,3]);

hist.f(step+1,:,1) =  min(f,[],[2,3]);
hist.f(step+1,:,2) = mean(f,   [2,3]);
hist.f(step+1,:,3) =  max(f,[],[2,3]);


% residuals
hist.res(step+1,:,1) = vecnorm(vecnorm(upd_u,2,3),2,2)./(vecnorm(vecnorm(u,2,3),2,2)+NUM.TINY);
hist.res(step+1,:,2) = vecnorm(vecnorm(upd_w,2,3),2,2)./(vecnorm(vecnorm(w,2,3),2,2)+NUM.TINY);
hist.res(step+1,:,3) = vecnorm(vecnorm(upd_p,2,3),2,2)./(vecnorm(vecnorm(p,2,3),2,2)+NUM.TINY);
hist.res(step+1,:,4) = vecnorm(vecnorm(upd_f,2,3),2,2)./(vecnorm(vecnorm(f,2,3),2,2)+NUM.TINY);

if step>0, hist.restot(step+1) = res; end
