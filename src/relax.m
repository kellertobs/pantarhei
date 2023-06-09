function [u, w, p, upd_u, upd_w, upd_p] = relax(u, w, p, Qvx, Qvz, Qf, f, PHS, NUM, minlvl)

if minlvl
    atol   = NUM.TINY;
    rtol   = NUM.TINY;
    maxits = NUM.minlvl_maxits;
else
    atol   = NUM.atol;
    rtol   = NUM.rtol;
    maxits = NUM.smth;
end

alpha  = NUM.alpha;
beta   = NUM.beta;

ui = u; 
wi = w; 
pi = p;

% update closures
closures;

resnorm = 1; resnorm0 = 0.1; it = 1;
while resnorm/resnorm0 >= rtol && resnorm >= atol && it < maxits

    % store solution of previous two iterations
    uii = ui;  ui = u;
    wii = wi;  wi = w;
    pii = pi;  pi = p;

    if ~mod(it,50); closures; end

    % update residuals
    residuals;

    % get solution updates
    upd_u = res_u.*dtau_u;
    upd_w = res_w.*dtau_w;
    upd_p = res_p.*dtau_p;

    upd_p = upd_p - mean(upd_p(:));
    if any(strcmp(NUM.BC, 'periodic'))
        upd_u = upd_u - mean(upd_u(:));
        upd_w = upd_w - mean(upd_w(:));
    end

    % apply updates
    u = ui - alpha.*upd_u + beta.*(ui-uii);
    w = wi - alpha.*upd_w + beta.*(wi-wii);
    p = pi - alpha.*upd_p + beta.*(pi-pii);

    % get residual norm
    if ~mod(it-1,10)
        resnorm = norm(upd_w,'fro')./(norm(w,'fro')+1e-16) ...
                + norm(upd_p,'fro')./(norm(p,'fro')+1e-16);
    end
    if it==1; resnorm0 = resnorm; end

    it = it+1;
end
end