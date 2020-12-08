
% calculate mms source term
% may want to put some limits for numerical stability.
usrc        = MMSsource('Calc_usrc' ,time,XuGrid(1,:,:),ZuGrid(1,:,:),d0,eta0,rho0,rhomix,grav,A,B,C,Gm, Amf,dmf,Tmf,Xmf,Zmf, thtlim, cfflim);
wsrc        = MMSsource('Calc_wsrc' ,time,XwGrid(1,:,:),ZwGrid(1,:,:),d0,eta0,rho0,rhomix,grav,A,B,C,Gm, Amf,dmf,Tmf,Xmf,Zmf, thtlim, cfflim);
[psrc,fsrc] = MMSsource('Calc_pfsrc',time,     X(1,:,:),     Z(1,:,:),d0,eta0,rho0,rhomix,grav,A,B,C,Gm, Amf,dmf,Tmf,Xmf,Zmf, thtlim, cfflim);


res_u = res_u - usrc;
res_w = res_w - wsrc;
res_p = res_p - psrc;
res_f = res_f - fsrc;

