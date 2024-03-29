function varargout = MMSsource (varargin)
% calculate MMS source
[varargout{1:nargout}] = feval(varargin{:});
end




% -------------------------------------------------------------------------
% functions to calculate sources
% (need separate functions because u, w, p/f are on different grids

function [usrc, dqvx, Gvx, Qvx, dqvxxdx, dqvxzdz] = Calc_usrc (t, x, z, ...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm, Amf,dmf,Tmf,Xmf,Zmf, thtlim, cfflim)
% on the vertical cell faces where u is defined

[f, dfdt, dfdx, dfdz                      ] = Calc_f(t,x,z,Tmf(:,1),Xmf(:,1),Zmf(:,1),Amf(:,1),dmf(:,1));
[p, ~   , dpdx, ~                         ] = Calc_p(t,x,z,Tmf(:,2),Xmf(:,2),Zmf(:,2),Amf(:,2),dmf(:,2));
[u, ~   , dudx, dudz, dudx2, dudz2, ~     ] = Calc_u(t,x,z,Tmf(:,3),Xmf(:,3),Zmf(:,3),Amf(:,3),dmf(:,3));
[~, ~   , dwdx, dwdz, ~    , ~    , dwdxdz] = Calc_w(t,x,z,Tmf(:,4),Xmf(:,4),Zmf(:,4),Amf(:,4),dmf(:,4));

rho    = rho0.*ones(size(p));

[Kv , Kf, Cv  , Cf  ] = closures_mms(f, eta0, d0, A, B, C, thtlim, cfflim);
[dKv                ] = closurederivs_mms(f, dfdt, dfdx, dfdz, d0, eta0, A, B, C);
[~  , ~ , omCv, omCf] = weights_mms(Kv, Kf, Cv, Cf);

% calculate reference fields
ustar = sum(omCv.*u, 1);
pstar = sum(omCf.*p, 1);

% derivatives of momentum flux
div_v    = dudx  + dwdz;
d_divv_x = dudx2 + dwdxdz;

% momentum flux
dqvxxdx  = - Kv.*(dudx2 - 1/3*d_divv_x)  - dKv(:,:,:,2).*(dudx  - 1/3*div_v) + f.*dpdx + p.*dfdx;
dqvxzdz  = -0.5*Kv.*(dwdxdz + dudz2) - 0.5*dKv(:,:,:,3).*(dwdx + dudz);
dqvx     = dqvxxdx + dqvxzdz;

% momentum transfer
Gvx  = Cv.*(u-ustar) - pstar.*dfdx;

% momentum source
Qvx = -f.*(rho-rhomix).*grav(2);

% sum it all
usrc = + dqvx + Gvx + Qvx;
end


function [wsrc, dqvz, Gvz, Qvz, gvz1, gvz2] = Calc_wsrc (t, x, z, ...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm, Amf,dmf,Tmf,Xmf,Zmf, thtlim, cfflim)
% on the horizontal cell faces where w is defined

[f, dfdt, dfdx, dfdz                      ] = Calc_f(t,x,z,Tmf(:,1),Xmf(:,1),Zmf(:,1),Amf(:,1),dmf(:,1));
[p, ~   , ~   , dpdz                      ] = Calc_p(t,x,z,Tmf(:,2),Xmf(:,2),Zmf(:,2),Amf(:,2),dmf(:,2));
[~, ~   , dudx, dudz, ~    , ~    , dudxdz] = Calc_u(t,x,z,Tmf(:,3),Xmf(:,3),Zmf(:,3),Amf(:,3),dmf(:,3));
[w, ~   , dwdx, dwdz, dwdx2, dwdz2, ~     ] = Calc_w(t,x,z,Tmf(:,4),Xmf(:,4),Zmf(:,4),Amf(:,4),dmf(:,4));

rho    = rho0.*ones(size(p));

[Kv , Kf, Cv  , Cf  ] = closures_mms(f, eta0, d0, A, B, C, thtlim, cfflim);
[dKv                ] = closurederivs_mms(f, dfdt, dfdx, dfdz, d0, eta0, A, B, C);
[~  , ~ , omCv, omCf] = weights_mms(Kv, Kf, Cv, Cf);

% update reference fields
wstar = sum(omCv.*w, 1);
pstar = sum(omCf.*p, 1);

% derivatives of momentum flux
div_v    = dudx   + dwdz;
d_divv_z = dudxdz + dwdz2;

dqvxzdx = -0.5*Kv.*(dwdx2 + dudxdz)    - 0.5*dKv(:,:,:,2).*(dwdx + dudz);
dqvzzdz = -    Kv.*(dwdz2 - 1/3.*d_divv_z) - dKv(:,:,:,3).*(dwdz - 1/3.*div_v) + f.*dpdz + p.*dfdz;
dqvz = dqvxzdx + dqvzzdz;

% momentum transfer
gvz1 =  Cv.*(w-wstar);
gvz2 = pstar.*dfdz;
Gvz = Cv.*(w-wstar) - pstar.*dfdz;

% momentum source
Qvz = -f.*(rho-rhomix).*grav(1);

wsrc = + dqvz + Gvz + Qvz;

end


function [psrc, fsrc, adv] = Calc_pfsrc (t, x, z, ...
    d0,eta0,rho0,rhomix,grav,A,B,C,Gm, Amf,dmf,Tmf,Xmf,Zmf, thtlim, cfflim)
% on the cell centers where p and f are defined

[f, dfdt, dfdx, dfdz                 ] = Calc_f(t,x,z,Tmf(:,1),Xmf(:,1),Zmf(:,1),Amf(:,1),dmf(:,1));
[p, ~   , dpdx, dpdz, dpdx2, dpdz2   ] = Calc_p(t,x,z,Tmf(:,2),Xmf(:,2),Zmf(:,2),Amf(:,2),dmf(:,2));
[u, ~   , dudx, ~   , ~    , ~    , ~] = Calc_u(t,x,z,Tmf(:,3),Xmf(:,3),Zmf(:,3),Amf(:,3),dmf(:,3));
[w, ~   , ~   , dwdz, ~    , ~    , ~] = Calc_w(t,x,z,Tmf(:,4),Xmf(:,4),Zmf(:,4),Amf(:,4),dmf(:,4));

rho = rho0.*ones(size(p));

[Kv , Kf  , Cv  , Cf         ] = closures_mms(f, eta0, d0, A, B, C, thtlim, cfflim);
[~  , dKf                    ] = closurederivs_mms(f, dfdt, dfdx, dfdz, d0, eta0, A, B, C);
[~  , omKf, omCv, omCf, domKf] = weights_mms(Kv, Kf, Cv, Cf, dKf);

% reference fields
ustar     = sum(omCv.*u   , 1);
wstar     = sum(omCv.*w   , 1);
pstar     = sum(omCf.*p   , 1);
rhostar   = sum(omCf.*rho , 1);
dpdxstar  = sum(omKf.*dpdx, 1);
dpdzstar  = sum(omKf.*dpdz, 1);
dpdx2star = sum(omKf.*dpdx2 + domKf(:,:,:,2).*dpdx, 1);
dpdz2star = sum(omKf.*dpdz2 + domKf(:,:,:,3).*dpdz, 1);

dqfxdx = -Kf.*(dpdx2 - dpdx2star) - dKf(:,:,:,2).*(dpdx - dpdxstar) + f.*dudx + u.*dfdx;
dqfzdz = -Kf.*(dpdz2 - dpdz2star) - dKf(:,:,:,3).*(dpdz - dpdzstar) + f.*dwdz + w.*dfdz;
dqf    = dqfxdx + dqfzdz;

adv = ustar.*dfdx + wstar.*dfdz;
Gf = Cf.*(p-pstar) - ustar.*dfdx - wstar.*dfdz - Gm./rhostar;

psrc =       dqf + Gf + Gm./rho;
fsrc = dfdt      - Gf;

end







% -------------------------------------------------------------------------
% functions to calculate permissions and weights

function [omKv, omKf, omCv, omCf, domKf] = weights_mms (Kv, Kf, Cv, Cf, dKf)

omKv = Kv./sum(Kv,1);
omKf = Kf./sum(Kf,1);
omCv = Cv./sum(Cv,1);
omCf = Cf./sum(Cf,1);

if nargout>4
    Kfsum = sum(Kf,1);
    domKf = dKf./Kfsum - Kf.*(sum(dKf,1))./(Kfsum.^2);
end

end


function [Kv, Kf, Cv, Cf, Kvcut, Kfcut] = closures_mms (f, eta0, d0, A, B, C, thtlim, cfflim)

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

% calculate cutoffs
Kvcut = max(geomean(Kv,[2,3]))./cfflim;
Kfcut = max(geomean(Kf,[2,3]))./cfflim;

% apply cutoff to coefficients to safeguard numerical stability
Kv = Kv + Kvcut;
Kf = Kf + Kfcut;
Cv = 1./(1./Cv + 1./(min(geomean(Cv,[2,3])).*cfflim));
Cf = 1./(1./Cf + 1./(min(geomean(Cf,[2,3])).*cfflim));
end






% -------------------------------------------------------------------------
% functions to calculate permissions derivatives

function [dKv, dKf] = closurederivs_mms (f, dfdt, dfdx, dfdz, d0, eta0, A, B, C)

NPHS = size(f,1);

f1    = f(1,:,:);   f2    = f(2,:,:);
d01   = d0(1);      d02   = d0(2);
eta01 = eta0(1);    eta02 = eta0(2); 

% fitting parameters for phase permissions
A1_1 = A(1,1);      A1_2 = A(1,2);      
A2_1 = A(2,1);      A2_2 = A(2,2);      

B1_1 = B(1,1);      B1_2 = B(1,2);
B2_1 = B(2,1);      B2_2 = B(2,2);

C1_1 = C(1,1);      C1_2 = C(1,2);   
C2_1 = C(2,1);      C2_2 = C(2,2);     

if NPHS==2
    
    if nargout == 1
        
        [dKvdf1,dKvdf2] = coeff_derivs_twophase(...
            f1,f2,d01,d02,eta01,eta02,...
            A1_1,A1_2,A2_1,A2_2,B1_1,B1_2,B2_1,B2_2,C1_1,C1_2,C2_1,C2_2);
        
        % change to derivs wrt t, x, z
        dKvdt = dKvdf1.*dfdt(1,:,:) + dKvdf2.*dfdt(2,:,:);
        dKvdx = dKvdf1.*dfdx(1,:,:) + dKvdf2.*dfdx(2,:,:);
        dKvdz = dKvdf1.*dfdz(1,:,:) + dKvdf2.*dfdz(2,:,:);
        dKv   = cat(4, dKvdt, dKvdx, dKvdz);
    else
        
        [dKvdf1,dKvdf2,dKfdf1,dKfdf2] = coeff_derivs_twophase(...
            f1,f2,d01,d02,eta01,eta02,...
            A1_1,A1_2,A2_1,A2_2,B1_1,B1_2,B2_1,B2_2,C1_1,C1_2,C2_1,C2_2);
        
        % change to derivs wrt t, x, z
        dKvdt = dKvdf1.*dfdt(1,:,:) + dKvdf2.*dfdt(2,:,:);
        dKvdx = dKvdf1.*dfdx(1,:,:) + dKvdf2.*dfdx(2,:,:);
        dKvdz = dKvdf1.*dfdz(1,:,:) + dKvdf2.*dfdz(2,:,:);
        dKv   = cat(4, dKvdt, dKvdx, dKvdz);
        
        dKfdt = dKfdf1.*dfdt(1,:,:) + dKfdf2.*dfdt(2,:,:);
        dKfdx = dKfdf1.*dfdx(1,:,:) + dKfdf2.*dfdx(2,:,:);
        dKfdz = dKfdf1.*dfdz(1,:,:) + dKfdf2.*dfdz(2,:,:);
        dKf   = cat(4, dKfdt, dKfdx, dKfdz);
    end
    
elseif NPHS==3
    
    f3    = f(3,:,:);   d03   = d0(3);      eta03 = eta0(3);
    
    A1_3  = A(1,3);     A2_3  = A(2,3);  
    A3_1  = A(3,1);     A3_2  = A(3,2);     A3_3  = A(3,3);

    B1_3  = B(1,3);     B2_3  = B(2,3);
    B3_1  = B(3,1);     B3_2  = B(3,2);     B3_3  = B(3,3);
    
    C1_3  = C(1,3);     C2_3  = C(2,3);
    C3_1  = C(3,1);     C3_2  = C(3,2);     C3_3  = C(3,3);
    
    % get derivatives for flux and transfer coeffs
    if nargout == 1
        
        [dKvdf1,dKvdf2,dKvdf3] = coeff_derivs_threephase(...
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
        
        [dKvdf1,dKvdf2,dKvdf3,dKfdf1,dKfdf2,dKfdf3] = coeff_derivs_threephase(...
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

end











% -------------------------------------------------------------------------
% mms source function initialisations

function [f, dfdt, dfdx, dfdz] = Calc_f (t, x, z, T, X, Z, f0, df0)
% calculate phase fractions (cosines)

tfrac    = t./T;
xfrac    = x./X;
zfrac    = z./Z;
f        = f0 + df0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
f(end,:,:) = 1 - sum(f(1:end-1,:,:),1);

if nargout>1
    dfdt = -1./T.*df0.*sin(tfrac).*cos(xfrac).*cos(zfrac);
    dfdx = -1./X.*df0.*cos(tfrac).*sin(xfrac).*cos(zfrac);
    dfdz = -1./Z.*df0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
    
    dfdt(end,:,:) = - sum(dfdt(1:end-1,:,:),1);
    dfdx(end,:,:) = - sum(dfdx(1:end-1,:,:),1);
    dfdz(end,:,:) = - sum(dfdz(1:end-1,:,:),1);
end
end


function [p, dpdt, dpdx, dpdz, dpdx2, dpdz2] = Calc_p (t, x, z, T, X, Z, P0, dP0)
% pressure (cosines)

tfrac = t./T;
xfrac = x./X;
zfrac = z./Z;
p     = P0 + dP0.*cos(tfrac).*cos(xfrac).*sin(zfrac);

if nargout>1
    dpdt = -1./T.*dP0.*sin(tfrac).*cos(xfrac).*sin(zfrac);
    dpdx = -1./X.*dP0.*cos(tfrac).*sin(xfrac).*sin(zfrac);
    dpdz =  1./Z.*dP0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
end

if nargout>4
    dpdx2  = -1./(X.*X).*dP0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
    dpdz2  = -1./(Z.*Z).*dP0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
end

end


function [u, dudt, dudx, dudz, dudx2, dudz2, dudxdz] = Calc_u (t, x, z, T, X, Z, U0, dU0)
% horizontal velocity (cosines)

tfrac = t./T;
xfrac = x./X;
zfrac = z./Z;
u     = U0 + dU0.*cos(tfrac).*sin(xfrac).*sin(zfrac);

if nargout>1
    dudt = -1./T.*dU0.*sin(tfrac).*sin(xfrac).*sin(zfrac);
    dudx =  1./X.*dU0.*cos(tfrac).*cos(xfrac).*sin(zfrac);
    dudz =  1./Z.*dU0.*cos(tfrac).*sin(xfrac).*cos(zfrac);
end

if nargout>4
    dudx2  = -1./(X.*X).*dU0.*cos(tfrac).*sin(xfrac).*sin(zfrac);
    dudz2  = -1./(Z.*Z).*dU0.*cos(tfrac).*sin(xfrac).*sin(zfrac);
    dudxdz =  1./(X.*Z).*dU0.*cos(tfrac).*cos(xfrac).*cos(zfrac);
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