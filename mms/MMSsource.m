function varargout = MMSsource (varargin)
% calculate MMS source
[varargout{1:nargout}] = feval(varargin{:});
end

function [usrc, dqvx, qvxx, qvxz, dqvxxdx, dqvxzdz, Gvx, Qvx] = Calc_usrc (t, x, z, ...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm,...
    Tp,Xp,Zp,P0,dP0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0,dU0, Tw,Xw,Zw,W0,dW0,...
    thtlim, cfflim)

[f, dfdt, dfdx, dfdz                      ] = Calc_f(t,x,z,Tf,Xf,Zf,f0,df0);
[p, ~   , dpdx, ~                         ] = Calc_p(t,x,z,Tp,Xp,Zp,P0,dP0);
[u, ~   , dudx, dudz, dudx2, dudz2, ~     ] = Calc_u(t,x,z,Tu,Xu,Zu,U0,dU0);
[~, ~   , dwdx, dwdz, ~    , ~    , dwdxdz] = Calc_w(t,x,z,Tw,Xw,Zw,W0,dW0);

rho    = rho0.*ones(size(p)); 

[Kv , Kf, Cv  , Cf  ] = CalcPermissions(f, eta0, d0, A, B, C, thtlim, cfflim);
[dKv                ] = GetPermissionDerivs(f, dfdt, dfdx, dfdz, d0, eta0, A, B, C);
[dKv                ] = CutPermissionDerivs(Kv, dKv, cfflim);
[~  , ~ , omCv, omCf] = CalcWeights(Kv, Kf, Cv, Cf);

% calculate reference fields
ustar = sum(omCv.*u, 1);
pstar = sum(omCf.*p, 1);

% derivatives of momentum flux
div_v    = dudx  + dwdz;
d_divv_x = dudx2 + dwdxdz;

qvxx = -Kv.*(dudx - 1./3*div_v) + f.*p;
qvxz = -0.5*Kv.*(dudz + dwdx);

% momentum flux
dqvxxdx  = - Kv.*          (dudx2 - 1/3*d_divv_x) ...
           - dKv(:,:,:,2).*(dudx  - 1/3*div_v) + f.*dpdx + p.*dfdx;
dqvxzdz  = -0.5*Kv.*(dwdxdz + dudz2) - 0.5*dKv(:,:,:,3).*(dwdx + dudz);
dqvx     = dqvxxdx + dqvxzdz;

% momentum transfer
Gvx   = - Cv.*(u-ustar) + pstar.*dfdx;

% momentum source
Qvx = f.*(rho-rhomix).*grav(2); 

% sum it all
usrc = + dqvx - Gvx - Qvx;
end

function [wsrc] = Calc_wsrc (t, x, z, ...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm,...
    Tp,Xp,Zp,P0,dP0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0,dU0, Tw,Xw,Zw,W0,dW0,...
    thtlim, cfflim)

[f, dfdt, dfdx, dfdz                      ] = Calc_f(t,x,z,Tf,Xf,Zf,f0,df0);
[p, ~   , ~   , dpdz                      ] = Calc_p(t,x,z,Tp,Xp,Zp,P0,dP0);
[~, ~   , dudx, dudz, ~    , ~    , dudxdz] = Calc_u(t,x,z,Tu,Xu,Zu,U0,dU0);
[w, ~   , dwdx, dwdz, dwdx2, dwdz2, ~     ] = Calc_w(t,x,z,Tw,Xw,Zw,W0,dW0);

rho    = rho0.*ones(size(p)); 

[Kv , Kf, Cv  , Cf  ] = CalcPermissions(f, eta0, d0, A, B, C, thtlim, cfflim);
[dKv                ] = GetPermissionDerivs(f, dfdt, dfdx, dfdz, d0, eta0, A, B, C);
[dKv                ] = CutPermissionDerivs(Kv, dKv, cfflim);
[~  , ~ , omCv, omCf] = CalcWeights(Kv, Kf, Cv, Cf);

% update reference fields
wstar = sum(omCv.*w, 1);
pstar = sum(omCf.*p, 1);

% derivatives of momentum flux
divv     = dudx   + dwdz;
d_divv_z = dudxdz + dwdz2;
dqvxzdx = -0.5*Kv.*(dwdx2 + dudxdz) - 0.5*dKv(:,:,:,2).*(dwdx + dudz);
dqvzzdz = -    Kv         .*(dwdz2 - 1/3.*d_divv_z) ...
          -   dKv(:,:,:,3).*(dwdz  - 1/3.*divv) + f.*dpdz + p.*dfdz;
dqvz = dqvxzdx + dqvzzdz;

% momentum transfer
Gvz = - Cv.*(w-wstar) + pstar.*dfdz;

% momentum source
Qvz = f.*(rho-rhomix).*grav(1); 

wsrc = + dqvz - Gvz - Qvz;

end

function [psrc, fsrc, dqfxdx, dqfzdz, Gf] = Calc_pfsrc (t, x, z, ...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm,...
    Tp,Xp,Zp,P0,dP0, Tf,Xf,Zf,f0,df0, Tu,Xu,Zu,U0,dU0, Tw,Xw,Zw,W0,dW0,...
    thtlim, cfflim)

[f, dfdt, dfdx, dfdz                 ] = Calc_f(t,x,z,Tf,Xf,Zf,f0,df0);
[p, ~   , dpdx, dpdz, dpdx2, dpdz2   ] = Calc_p(t,x,z,Tp,Xp,Zp,P0,dP0);
[u, ~   , dudx, ~   , ~    , ~    , ~] = Calc_u(t,x,z,Tu,Xu,Zu,U0,dU0);
[w, ~   , ~   , dwdz, ~    , ~    , ~] = Calc_w(t,x,z,Tw,Xw,Zw,W0,dW0);

rho = rho0.*ones(size(p)); 

[Kv , Kf  , Cv  , Cf         ] = CalcPermissions(f, eta0, d0, A, B, C, thtlim, cfflim);
[~  , dKf                    ] = GetPermissionDerivs(f, dfdt, dfdx, dfdz, d0, eta0, A, B, C);
[dKf                         ] = CutPermissionDerivs(Kf, dKf, cfflim);
[~  , omKf, omCv, omCf, domKf] = CalcWeights(Kv, Kf, Cv, Cf, dKf);

% reference fields
ustar     = sum(omCv.*u   , 1);
wstar     = sum(omCv.*w   , 1);
pstar     = sum(omCf.*p   , 1);
rhostar   = sum(omCf.*rho , 1);
dpdxstar  = sum(omKf.*dpdx, 1);
dpdzstar  = sum(omKf.*dpdz, 1);
dpdx2star = sum(omKf.*dpdx2 + domKf(:,:,:,2).*dpdx, 1);
dpdz2star = sum(omKf.*dpdz2 + domKf(:,:,:,3).*dpdz, 1);

% qfx = -Kf.*(dpdx - dpdxstar) + f.*u;
% qfz = -Kf.*(dpdz - dpdzstar) + f.*w;

dqfxdx = -Kf.*(dpdx2 - dpdx2star) - dKf(:,:,:,2).*(dpdx - dpdxstar) + f.*dudx + u.*dfdx;
dqfzdz = -Kf.*(dpdz2 - dpdz2star) - dKf(:,:,:,3).*(dpdz - dpdzstar) + f.*dwdz + w.*dfdz;

Gf = - Cf.*(p-pstar) + ustar.*dfdx + wstar.*dfdz + Gm./rhostar;

psrc =       dqfxdx + dqfzdz - Gf - Gm./rho;
fsrc = dfdt                  + Gf;

end

function [omKv, omKf, omCv, omCf, domKf] = CalcWeights (Kv, Kf, Cv, Cf, dKf)

omKv = Kv./sum(Kv,1);
omKf = Kf./sum(Kf,1);
omCv = Cv./sum(Cv,1);
omCf = Cf./sum(Cf,1);

if nargout>4
    Kfsum = sum(Kf,1);
    domKf = dKf./Kfsum - Kf.*(sum(dKf,1))./(Kfsum.^2);
end

end

function [Kv, Kf, Cv, Cf] = CalcPermissions (f, eta0, d0, A, B, C, thtlim, cfflim)

NPHS = length(eta0);

% pure-phase diffusion parameters
kv = eta0;
kf = d0.^2./eta0;

% get diffusivity contrasts
Mv = kv.'./kv;          % momentum diffusivity ratios
Mf = kf.'./kf;          % volume diffusivity ratios

% permission weights
F  = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
Sf = (F./B).^(1./C);
Sf = Sf./sum(Sf,2);
Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;

% permissions
thtv = squeeze(prod(Mv.^Xf,2));  gmv = geomean(geomean(thtv,3),2);
thtf = squeeze(prod(Mf.^Xf,2));  gmf = geomean(geomean(thtf,3),2);
thtv = 1./(1./thtv + 1./(thtlim.^0.5*gmv)) + (gmv/thtlim.^0.5);
thtf = 1./(1./thtf + 1./(thtlim.^0.5*gmf)) + (gmf/thtlim.^0.5);

% get momentum and volume flux and transfer coefficients
Kv =    f .*eta0       .*thtv;
Kf =    f .*d0.^2./eta0.*thtf;
Cv = (1-f)./d0.^2.*Kv;
Cf = (1-f)./d0.^2.*Kf;

% apply cutoff to coefficients to safeguard numerical stability
Kv = Kv + max(Kv,[],1)./cfflim;
Kf = Kf + max(Kf,[],1)./cfflim;
Cv = 1./(1./Cv + 1./(min(Cv,[],1).*cfflim));
Cf = 1./(1./Cf + 1./(min(Cf,[],1).*cfflim));
end

function [dKout] = CutPermissionDerivs (K, dKin, cfflim)

[~, ind] = max(K,[],1,'linear');

dKout = dKin;
for i = 1:3
    tmp = dKin(:,:,:,i);
    dKout(:,:,:,i) = tmp + tmp(ind)./cfflim;
end

end

function [dKv, dKf] = GetPermissionDerivs (f, dfdt, dfdx, dfdz, d0, eta0, A, B, C)

f1    = f(1,:,:);   f2    = f(2,:,:);   f3    = f(3,:,:);
d01   = d0(1);      d02   = d0(2);      d03   = d0(3);
eta01 = eta0(1);    eta02 = eta0(2);    eta03 = eta0(3);

% fitting parameters for phase permissions---------------------------------
A1_1 = A(1,1);      A1_2 = A(1,2);      A1_3 = A(1,3);
A2_1 = A(2,1);      A2_2 = A(2,2);      A2_3 = A(2,3);
A3_1 = A(3,1);      A3_2 = A(3,2);      A3_3 = A(3,3);

B1_1 = B(1,1);      B1_2 = B(1,2);      B1_3 = B(1,3);
B2_1 = B(2,1);      B2_2 = B(2,2);      B2_3 = B(2,3);
B3_1 = B(3,1);      B3_2 = B(3,2);      B3_3 = B(3,3);

C1_1 = C(1,1);      C1_2 = C(1,2);      C1_3 = C(1,3);
C2_1 = C(2,1);      C2_2 = C(2,2);      C2_3 = C(2,3);
C3_1 = C(3,1);      C3_2 = C(3,2);      C3_3 = C(3,3);

% get derivatives for flux and transfer coeffs-----------------------------
if nargout == 1
    
    [dKvdf1,dKvdf2,dKvdf3] = coeff_derivs(...
        f1,f2,f3,d01,d02,d03,eta01,eta02,eta03,...
        A1_1,A1_2,A1_3,A2_1,A2_2,A2_3,A3_1,A3_2,A3_3,...
        B1_1,B1_2,B1_3,B2_1,B2_2,B2_3,B3_1,B3_2,B3_3,...
        C1_1,C1_2,C1_3,C2_1,C2_2,C2_3,C3_1,C3_2,C3_3);
    
    % change to derivs wrt t, x, z
    dKvdt = dKvdf1.*dfdt(1,:,:) + dKvdf2.*dfdt(2,:,:) + dKvdf3.*dfdt(3,:,:);
    dKvdx = dKvdf1.*dfdx(1,:,:) + dKvdf2.*dfdx(2,:,:) + dKvdf3.*dfdx(3,:,:);
    dKvdz = dKvdf1.*dfdz(1,:,:) + dKvdf2.*dfdz(2,:,:) + dKvdf3.*dfdz(3,:,:);
    dKv   = cat(4, dKvdt, dKvdx, dKvdz);
else
    
    [dKvdf1,dKvdf2,dKvdf3,dKfdf1,dKfdf2,dKfdf3] = coeff_derivs(...
        f1,f2,f3,d01,d02,d03,eta01,eta02,eta03,...
        A1_1,A1_2,A1_3,A2_1,A2_2,A2_3,A3_1,A3_2,A3_3,...
        B1_1,B1_2,B1_3,B2_1,B2_2,B2_3,B3_1,B3_2,B3_3,...
        C1_1,C1_2,C1_3,C2_1,C2_2,C2_3,C3_1,C3_2,C3_3);
    
    % change to derivs wrt t, x, z
    dKvdt = dKvdf1.*dfdt(1,:,:) + dKvdf2.*dfdt(2,:,:) + dKvdf3.*dfdt(3,:,:);
    dKvdx = dKvdf1.*dfdx(1,:,:) + dKvdf2.*dfdx(2,:,:) + dKvdf3.*dfdx(3,:,:);
    dKvdz = dKvdf1.*dfdz(1,:,:) + dKvdf2.*dfdz(2,:,:) + dKvdf3.*dfdz(3,:,:);
    dKv   = cat(4, dKvdt, dKvdx, dKvdz);
    
    dKfdt = dKfdf1.*dfdt(1,:,:) + dKfdf2.*dfdt(2,:,:) + dKfdf3.*dfdt(3,:,:);
    dKfdx = dKfdf1.*dfdx(1,:,:) + dKfdf2.*dfdx(2,:,:) + dKfdf3.*dfdx(3,:,:);
    dKfdz = dKfdf1.*dfdz(1,:,:) + dKfdf2.*dfdz(2,:,:) + dKfdf3.*dfdz(3,:,:);
    dKf   = cat(4, dKfdt, dKfdx, dKfdz);
end

end

function [X,Z] = MakeGrid (D,h, NPHS)
% initialise coordinate arrays
x     = -D/2+h/2 : h : D/2-h/2;
[X,~] = meshgrid(x,x);
X     = repmat(X,1,1,NPHS);
Z     = permute(X,[3,2,1]);
X     = permute(X,[3,1,2]);

end

function [f, dfdt, dfdx, dfdz] = Calc_f (t, x, z, T, X, Z, f0, df0)
% calculate phase fractions (cosines)

tfrac    = t./T;
xfrac    = x./X;
zfrac    = z./Z;
f        = f0 + df0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
f(2,:,:) = 1 - f(1,:,:) - f(3,:,:);

if nargout>1
    dfdt = -1./T.*df0.*sin(tfrac).*cos(xfrac).*cos(zfrac);
    dfdx = -1./X.*df0.*cos(tfrac).*sin(xfrac).*cos(zfrac);
    dfdz = -1./Z.*df0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
    
    dfdt(2,:,:) = - dfdt(1,:,:) - dfdt(3,:,:);
    dfdx(2,:,:) = - dfdx(1,:,:) - dfdx(3,:,:);
    dfdz(2,:,:) = - dfdz(1,:,:) - dfdz(3,:,:);
end
end

function [p, dpdt, dpdx, dpdz, dpdx2, dpdz2] = Calc_p (t, x, z, T, X, Z, P0, dP0)
% pressure (cosines)

tfrac = t./T;
xfrac = x./X;
zfrac = z./Z;
p     = P0 + dP0.*cos(tfrac).*cos(xfrac).*cos(zfrac);

if nargout>1
    dpdt = -1./T.*dP0.*sin(tfrac).*cos(xfrac).*cos(zfrac);
    dpdx = -1./X.*dP0.*cos(tfrac).*sin(xfrac).*cos(zfrac);
    dpdz = -1./Z.*dP0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
end

if nargout>4
    dpdx2  = -1./(X.*X).*dP0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
    dpdz2  = -1./(Z.*Z).*dP0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
end

end

function [u, dudt, dudx, dudz, dudx2, dudz2, dudxdz] = Calc_u (t, x, z, T, X, Z, U0, dU0)
% horizontal velocity (cosines)

tfrac = t./T;
xfrac = x./X;
zfrac = z./Z;
u     = U0 + dU0.*cos(tfrac).*cos(xfrac).*cos(zfrac);

if nargout>1
    dudt = -1./T.*dU0.*sin(tfrac).*cos(xfrac).*cos(zfrac);
    dudx = -1./X.*dU0.*cos(tfrac).*sin(xfrac).*cos(zfrac);
    dudz = -1./Z.*dU0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
end

if nargout>4
    dudx2  = -1./(X.*X).*dU0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
    dudz2  = -1./(Z.*Z).*dU0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
    dudxdz =  1./(X.*Z).*dU0.*cos(tfrac).*sin(xfrac).*sin(zfrac);
end

end

function [w, dwdt, dwdx, dwdz, dwdx2, dwdz2, dwdxdz] = Calc_w (t, x, z, T, X, Z, W0, dW0)
% vertical velocity (cosines)

tfrac = t./T;
xfrac = x./X;
zfrac = z./Z;
w     = W0 + dW0.*cos(tfrac).*cos(xfrac).*cos(zfrac);

if nargout>1
    dwdt = -1./T.*dW0.*sin(tfrac).*cos(xfrac).*cos(zfrac);
    dwdx = -1./X.*dW0.*cos(tfrac).*sin(xfrac).*cos(zfrac);
    dwdz = -1./Z.*dW0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
end

if nargout>4
    dwdx2  = -1./(X.*X).*dW0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
    dwdz2  = -1./(Z.*Z).*dW0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
    dwdxdz =  1./(X.*Z).*dW0.*cos(tfrac).*sin(xfrac).*sin(zfrac);
end

end