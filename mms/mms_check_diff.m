
% initialise coefficient and auxiliary fields
run('../src/up2date');

% update reference and auxiliary fields
ustar     = sum(omvx.*u,1);
wstar     = sum(omvz.*w,1);
pstar     = sum(omfc.*p,1);
rhostar   = sum(omfc.*rho,1);
Gx_pstar  = sum(omfx.*diff(p(:,:,ic),1,3)./h,1);
Gz_pstar  = sum(omfz.*diff(p(:,ic,:),1,2)./h,1);
Div_v     = diff(u,1,3)./h + diff(w,1,2)./h;
pstar_Gfx = (pstar(:,:,im)+pstar(:,:,ip))./2.*diff(f(:,:,ic),1,3)./h;
pstar_Gfz = (pstar(:,im,:)+pstar(:,ip,:))./2.*diff(f(:,ic,:),1,2)./h;

% get momentum flux fields
qvxx  = - Kv .* (diff(u,1,3)./h - Div_v./3) + f.*p;
qvzz  = - Kv .* (diff(w,1,2)./h - Div_v./3) + f.*p;
qvxz  = -(Kv(:,im,im)+Kv(:,ip,im)+Kv(:,im,ip)+Kv(:,ip,ip))./4 ...
    .* (diff(u(:,ic,:),1,2)./h + diff(w(:,:,ic),1,3)./h)./2;
if strcmp(BC,'closed'); qvxz(:,[1,end],:) = 0; end

% get volume flux fields
qfx = - (Kf(:,:,im)+Kf(:,:,ip))./2 .* (diff(p(:,:,ic),1,3)./h - Gx_pstar) ...
    + ( f(:,:,im)+ f(:,:,ip))./2 .* u;
qfz = - (Kf(:,im,:)+Kf(:,ip,:))./2 .* (diff(p(:,ic,:),1,2)./h - Gz_pstar) ...
    + ( f(:,im,:)+ f(:,ip,:))./2 .* w;
if strcmp(BC,'closed'); qfx(:,:,[1,end]) = 0; qfz(:,[1,end],:) = 0; end

% get momentum transfer fields
Gvx = - (Cv(:,:,im)+Cv(:,:,ip))./2 .* (u-ustar) + pstar_Gfx;
Gvz = - (Cv(:,im,:)+Cv(:,ip,:))./2 .* (w-wstar) + pstar_Gfz;

%get volume transfer field
Gf  = -             Cf             .* (p-pstar) + vstar_Gf + Gm./rhostar;

% get momentum source fields
Qvx = (f(:,:,im)+f(:,:,ip))./2.*((rho(:,:,im)+rho(:,:,ip))./2-rhomix).*grav(2);
Qvz = (f(:,im,:)+f(:,ip,:))./2.*((rho(:,im,:)+rho(:,ip,:))./2-rhomix).*grav(1);

% get physical time step
dttmp = cfl/(max(abs([qfx(:);qfz(:)]))/(h/2) + max(abs(Gf(:)))./5e-3);  % [s]

% get residual fields
res_u =             + diff(qvxx(:,:,ic),1,3)./h + diff(qvxz,1,2)./h - Gvx - Qvx    ;
res_w =             + diff(qvzz(:,ic,:),1,2)./h + diff(qvxz,1,3)./h - Gvz - Qvz    ;
res_p =             + diff(qfx         ,1,3)./h + diff(qfz ,1,2)./h - Gf  - Gm./rho;
res_f =  (f-fo)./dt                                                 + Gf           ;

% calculate mms source term
% may want to put some limits for numerical stability.
usrc        =  usrcFunc(time);
wsrc        =  wsrcFunc(time);
[psrc,fsrc] = pfsrcFunc(time);

res_u = res_u - usrc;
res_w = res_w - wsrc;
res_p = res_p - psrc;
res_f = res_f - fsrc;

%%
% check residual norm
res_u0 = norm(res_u(:).*dtau_u(:),2)./(norm(u(:),2)+1e-32);
res_w0 = norm(res_w(:).*dtau_w(:),2)./(norm(w(:),2)+1e-32);
res_p0 = norm(res_p(:).*dtau_p(:),2)./(norm(p(:),2)+1e-32);
res_f0 = norm(res_f(:).*dtau_f(:),2)./(norm(f(:),2)+1e-32);
res_0  = res_u0 + res_w0 + res_p0 + res_f0;

fprintf(1, '\n\n');
fprintf(1, 'initial dt = %.4f.\n', dttmp);
fprintf(1, '\n\n');

fprintf(1, '(res_u x dtau_u)/u = %4.4e.\n', res_u0);
fprintf(1, '(res_w x dtau_w)/w = %4.4e.\n', res_w0);
fprintf(1, '(res_p x dtau_p)/p = %4.4e.\n', res_p0);
fprintf(1, '(res_f x dtau_f)/f = %4.4e.\n', res_f0);
fprintf(1, 'total residual     = %4.4e.\n', res_0);
fprintf(1, '\n\n');

fprintf(1, 'mean dtau_u = %4.4e.\n', mean(dtau_u(:)));
fprintf(1, 'mean dtau_w = %4.4e.\n', mean(dtau_w(:)));
fprintf(1, 'mean dtau_p = %4.4e.\n', mean(dtau_p(:)));
fprintf(1, 'mean dtau_f = %4.4e.\n', mean(dtau_f(:)));
fprintf(1, '\n\n');

fprintf(1, 'mean res_u = %4.4e.\n', mean(res_u(:).*dtau_u(:)));
fprintf(1, 'mean res_w = %4.4e.\n', mean(res_w(:).*dtau_w(:)));
fprintf(1, 'mean res_p = %4.4e.\n', mean(res_p(:).*dtau_p(:)));
fprintf(1, 'mean res_f = %4.4e.\n', mean(res_f(:).*dtau_f(:)));
fprintf(1, '\n\n');


