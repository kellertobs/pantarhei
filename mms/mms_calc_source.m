
% calculate mms source term
% may want to put some limits for numerical stability.
usrc        =  usrcFunc(time);
wsrc        =  wsrcFunc(time);
[psrc,fsrc] = pfsrcFunc(time);

% psrc (abs(psrc )<1e-14) = 1e-8;
% res_p(abs(res_p)<1e-14) = 1e-8;
% 
% fsrc (abs(fsrc )<1e-14) = 1e-8;
% res_f(abs(res_f)<1e-14) = 1e-8;

res_u = res_u - usrc;
res_w = res_w - wsrc;
res_p = res_p - psrc;
res_f = res_f - fsrc;

