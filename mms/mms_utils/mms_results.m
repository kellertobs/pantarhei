
%% print individual residuals

% generate manufactured solutions
fTrue = MMSsource('Calc_f', time, X     , Z     , Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
pTrue = MMSsource('Calc_p', time, X     , Z     , Tmf(:,2), Xmf(:,2), Zmf(:,2), Amf(:,2), dmf(:,2));
uTrue = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tmf(:,3), Xmf(:,3), Zmf(:,3), Amf(:,3), dmf(:,3));
wTrue = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tmf(:,4), Xmf(:,4), Zmf(:,4), Amf(:,4), dmf(:,4));

% 2-norm error
fNormErr = 100*norm(f(:)-fTrue(:),2)./norm(fTrue(:),2);
pNormErr = 100*norm(p(:)-pTrue(:),2)./norm(pTrue(:),2);
uNormErr = 100*norm(u(:)-uTrue(:),2)./norm(uTrue(:),2);
wNormErr = 100*norm(w(:)-wTrue(:),2)./norm(wTrue(:),2);

% maximum error
fMaxErr = 100*norm(f(:)-fTrue(:),inf)./norm(fTrue(:),inf);
pMaxErr = 100*norm(p(:)-pTrue(:),inf)./norm(pTrue(:),inf);
uMaxErr = 100*norm(u(:)-uTrue(:),inf)./norm(uTrue(:),inf);
wMaxErr = 100*norm(w(:)-wTrue(:),inf)./norm(wTrue(:),inf);

fprintf(1, '\n\n');
fprintf(1,'2-norm [f,p,u,w] percent error: \n%.6e, %.6f, %.6f, %.6f.\n\n', ...
    fNormErr, pNormErr, uNormErr, wNormErr); 
fprintf(1,'Max [f,p,u,w] percent error: \n%.6e, %.6f, %.6f, %.6f.', ...
    fMaxErr, pMaxErr, uMaxErr, wMaxErr);
fprintf(1, '\n\n');

