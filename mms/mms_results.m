addpath('../usr/');

figure; plot(tev, cos(tev./Tf(1)), time, cos(time./Tf(1)), 'rx');

%% print individual residuals

fprintf(1, '\n\n');
fprintf(1, '(res_u x dtau_u)/u = %4.4e.\n', norm(res_u(:).*dtau_u(:),2)./(norm(u(:),2)+1e-32));
fprintf(1, '(res_w x dtau_w)/w = %4.4e.\n', norm(res_w(:).*dtau_w(:),2)./(norm(w(:),2)+1e-32));
fprintf(1, '(res_p x dtau_p)/p = %4.4e.\n', norm(res_p(:).*dtau_p(:),2)./(norm(p(:),2)+1e-32));
fprintf(1, '(res_f x dtau_f)/f = %4.4e.\n', norm(res_f(:).*dtau_f(:),2)./(norm(f(:),2)+1e-32));
fprintf(1, '\n\n');

%% generate manufactured solutions

uTrue = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tu, Xu, Zu, U0, dU0);
wTrue = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tw, Xw, Zw, W0, dW0);
pTrue = MMSsource('Calc_p', time, X     , Z     , Tp, Xp, Zp, P0, dP0);
fTrue = MMSsource('Calc_f', time, X     , Z     , Tf, Xf, Zf, f0, df0);

fprintf(1, '\n\n');
fprintf(1,'norm u error = %.4f.\n', norm(uTrue(:)-u(:),2)./norm(uTrue(:),2));
fprintf(1,'norm w error = %.4f.\n', norm(wTrue(:)-w(:),2)./norm(wTrue(:),2));
fprintf(1,'norm p error = %.4f.\n', norm(pTrue(:)-p(:),2)./norm(pTrue(:),2));
fprintf(1,'norm f error = %.4f.\n', norm(fTrue(:)-f(:),2)./norm(fTrue(:),2));
fprintf(1, '\n\n');

%%

%{
plot_field(uTrue, 'uTrue'); plot_field(u, 'uNum');
plot_field(wTrue, 'wTrue'); plot_field(w, 'wNum');
plot_field(pTrue, 'pTrue'); plot_field(p, 'pNum');
plot_field(fTrue, 'fTrue'); plot_field(f, 'fNum');

plot_field((uTrue-u)./uTrue, 'u normalized error'); 
plot_field((wTrue-w)./wTrue, 'w normalized error'); 
plot_field((pTrue-p)./pTrue, 'p normalized error'); 
plot_field((fTrue-f)./fTrue, 'f normalized error'); 
%}



