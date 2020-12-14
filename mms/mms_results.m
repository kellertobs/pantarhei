
%% print individual residuals

% fprintf(1, '\n\n');
% fprintf(1, '(res_f x dtau_f)/f = %4.4e.\n', norm(res_f(:).*dtau_f(:),2)./(norm(f(:),2)+1e-32));
% fprintf(1, '(res_p x dtau_p)/p = %4.4e.\n', norm(res_p(:).*dtau_p(:),2)./(norm(p(:),2)+1e-32));
% fprintf(1, '(res_u x dtau_u)/u = %4.4e.\n', norm(res_u(:).*dtau_u(:),2)./(norm(u(:),2)+1e-32));
% fprintf(1, '(res_w x dtau_w)/w = %4.4e.\n', norm(res_w(:).*dtau_w(:),2)./(norm(w(:),2)+1e-32));
% fprintf(1, '\n\n');

% generate manufactured solutions
fTrue = MMSsource('Calc_f', time, X     , Z     , Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
pTrue = MMSsource('Calc_p', time, X     , Z     , Tmf(:,2), Xmf(:,2), Zmf(:,2), Amf(:,2), dmf(:,2));
uTrue = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tmf(:,3), Xmf(:,3), Zmf(:,3), Amf(:,3), dmf(:,3));
wTrue = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tmf(:,4), Xmf(:,4), Zmf(:,4), Amf(:,4), dmf(:,4));

fNormErr = 100*1./N.^2.*norm(reshape((f-fTrue)./dmf(:,1),[],1));
pNormErr = 100*1./N.^2.*norm(reshape((p-pTrue)./dmf(:,2),[],1));
uNormErr = 100*1./N.^2.*norm(reshape((u-uTrue)./dmf(:,3),[],1));
wNormErr = 100*1./N.^2.*norm(reshape((w-wTrue)./dmf(:,4),[],1));

fMaxErr = 100*norm(reshape((f-fTrue)./dmf(:,1),[],1), inf);
pMaxErr = 100*norm(reshape((p-pTrue)./dmf(:,2),[],1), inf);
uMaxErr = 100*norm(reshape((u-uTrue)./dmf(:,3),[],1), inf);
wMaxErr = 100*norm(reshape((w-wTrue)./dmf(:,4),[],1), inf);

fprintf(1, '\n\n');
fprintf(1,'2-norm [f,p,u,w] percent error: \n%.6e, %.6f, %.6f, %.6f.\n\n', ...
    fNormErr, pNormErr, uNormErr, wNormErr); 
fprintf(1,'Max [f,p,u,w] percent error: \n%.6e, %.6f, %.6f, %.6f.', ...
    fMaxErr, pMaxErr, uMaxErr, wMaxErr);
fprintf(1, '\n\n');

%%

%{
plot_field(fTrue, 'fTrue'); plot_field(f, 'fNum');

plot_field(pTrue, 'pTrue'); plot_field(p, 'pNum');
plot_field((pTrue-p)./dmf(:,2), ' p error');

plot_field(uTrue, 'uTrue'); plot_field(u, 'uNum'); 
plot_field((u-uTrue)./dmf(:,3), 'u error');

plot_field(wTrue, 'wTrue'); plot_field(w, 'wNum');
plot_field((w-wTrue)./dmf(:,4), 'w error');


%}



