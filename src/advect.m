function [adv] = advect (f, u, w, h, scheme, dim, BC)
%
% [adv] = advect (f, u, w, h, scheme, dim, BC)
%
% example:
% adv = advect(f, u, w, h, {'quick','vdf'}, [2,3], {'periodic','periodic'});
% 
% advects a quantity f on a velocity field represented by (u,w),
% for 2D problems
% 
% INPUTS
% f         quantity you want to advect [Nz x Nx, or Nphs x Nz x Nx]
% u         STAGGERED horiz velocity field [Nz x Nx+1, or Nphs x Nz x Nx+1]
% u         STAGGERED vert  velocity field [Nz+1 x Nx, or Nphs x Nz+1 x Nx]
% h         grid spacing
% scheme    info about advection scheme [2-element cell]
%               1st element is the advection scheme name
%               2nd element is '' (returns div(v x f)) or 'vdf' (returns v grad(f))
% dim       which dimensions correspond to z, x [2-element array]
% BC        boundary conditions on z,x [2-element cell]
% 
% OUTPUT
% adv       div(v x f) or v grad(f) depending on scheme{2}, same size as f
% 
% YQW, 11 Oct 2022
% 



% collect information on the dimensions and BCs corresponding to (z,x)
zdim = dim(1); xdim = dim(2);
zBC  = BC{1};  xBC  = BC{2};

fcc = f;
fdv = f.*(diff(u,1,xdim)./h + diff(w,1,zdim)./h);   % f x div(v)


switch scheme{1}
    case 'centr'        
        % centered differences, prone to num dispersion
        [um , up ] = facevels(u, xdim);
        [wm , wp ] = facevels(w, zdim);
        
        [fmx, fpx] = makestencil(f, xdim, xBC);
        [fmz, fpz] = makestencil(f, zdim, zBC);
        
        adv = ( up.*(fcc+fpx)./2 - um.*(fcc+fmx)./2 )./h + ...
              ( wp.*(fcc+fpz)./2 - wm.*(fcc+fmz)./2 )./h;
      
    case 'upwd1'
        % upwind differences, prone to num diffusion
        % NB upwind designed for v x grad(f), so we add f x div(v) here
        uc = centervels(u, xdim);
        wc = centervels(w, zdim);
        
        [fmx, fpx] = makestencil(f, xdim, xBC);
        [fmz, fpz] = makestencil(f, zdim, zBC);
        
        adv = max(uc,0).*(fcc-fmx)./h + min(uc,0).*(-fcc+fpx)./h + ...
              max(wc,0).*(fcc-fmz)./h + min(wc,0).*(-fcc+fpz)./h + fdv;
          
    case 'quick'
        % quick scheme == 3rd order upwind
        uc = centervels(u, xdim);
        wc = centervels(w, zdim);
        
        [fmx, fpx, fmmx, fppx] = makestencil(f, xdim, xBC);
        [fmz, fpz, fmmz, fppz] = makestencil(f, zdim, zBC);

        adv = max(uc,0).*( 2*fpx + 3*fcc - 6*fmx + fmmx)./6./h + ...
              min(uc,0).*(-2*fmx - 3*fcc + 6*fpx - fppx)./6./h + ...
              max(wc,0).*( 2*fpz + 3*fcc - 6*fmz + fmmz)./6./h + ...
              min(wc,0).*(-2*fmz - 3*fcc + 6*fpz - fppz)./6./h + fdv;
        
          
          
          
    % the schemes below here definitely work for periodic BCs, 
    % still not sure about closed BCs
    case 'fromm'
        % Fromm scheme that reduces oscillations
        [um, up] = facevels(u, xdim);
        [wm, wp] = facevels(w, zdim);
        
        [fmx, fpx, fmmx, fppx] = makestencil(f, xdim, xBC);
        [fmz, fpz, fmmz, fppz] = makestencil(f, zdim, zBC);
        
        adv = +     up .*( -(fppx-fpx)./8 + (fpx+fcc)./2 + (fcc-fmx )./8 )./h ...
              - abs(up).*( -(fppx-fpx)./8 + (fpx-fcc)./4 - (fcc-fmx )./8 )./h ...
              -     um .*( -(fpx -fcc)./8 + (fcc+fmx)./2 + (fmx-fmmx)./8 )./h ...
              + abs(um).*( -(fpx -fcc)./8 + (fcc-fmx)./4 - (fmx-fmmx)./8 )./h ...
              +     wp .*( -(fppz-fpz)./8 + (fpz+fcc)./2 + (fcc-fmz )./8 )./h ...
              - abs(wp).*( -(fppz-fpz)./8 + (fpz-fcc)./4 - (fcc-fmz )./8 )./h ...
              -     wm .*( -(fpz -fcc)./8 + (fcc+fmz)./2 + (fmz-fmmz)./8 )./h ...
              + abs(wm).*( -(fpz -fcc)./8 + (fcc-fmz)./4 - (fmz-fmmz)./8 )./h;
          
    case 'weno5'
        % 5th order WENO (slow) from Jiang & Shu 1996, J Comp Physics
        % TODO: figure out a more elegant way to make use of staggered grid
        [um, up] = facevels(u, xdim);
        [wm, wp] = facevels(w, zdim);
        
        xdflux = weno5(f, up, xdim, xBC);
        zdflux = weno5(f, wp, zdim, zBC);
        
        adv    = 1/h*(xdflux + zdflux);
        
    case 'tvdim'
        % total variation diminishing approach from Sramek et al. 2010, GJI
        % TODO: figure out a more elegant way to make use of staggered grid
        [um, up] = facevels(u, xdim);
        [wm, wp] = facevels(w, zdim);
        
        xdflux = tvd(f, up, xdim, xBC);
        zdflux = tvd(f, wp, zdim, zBC);
        
        adv    = 1/h*(xdflux + zdflux);
        
end

if strcmp(scheme{2}, 'vdf')
    % v x grad(f) = div (v x f) - f x div(v)
    adv = adv - fdv;    
end

end



%%  utility functions used for many schemes

function [vm, vp] = facevels (v, dim)
% return cell face velocities upwind (vm) and downwind (vp) of cell center
if     dim==1, vm = v(1:end-1,:,:);    vp = v(2:end,:,:);
elseif dim==2, vm = v(:,1:end-1,:);    vp = v(:,2:end,:);
elseif dim==3, vm = v(:,:,1:end-1);    vp = v(:,:,2:end);
end
end

function [vc] = centervels (v, dim)
% return cell-centered velocities by taking averages
if     dim==1, vc = 0.5*( v(1:end-1,:,:) + v(2:end,:,:) );
elseif dim==2, vc = 0.5*( v(:,1:end-1,:) + v(:,2:end,:) );
elseif dim==3, vc = 0.5*( v(:,:,1:end-1) + v(:,:,2:end) );
end
end

function [flux] = shiftflux (flux, shift, dim, BC)
% function needed to deal with circshift when applied to closed BCs
% fluxes = 0 at closed boundaries 

flux = circshift(flux, shift);

if strcmp(BC, 'closed')
    % which direction is the shift? + is right, - is left
    dir = sign(sum(shift)); 
    if dir>0
        if     dim==1, flux(1,:,:) = 0;
        elseif dim==2, flux(:,1,:) = 0;
        elseif dim==3, flux(:,:,1) = 0;
        end
    else
        if     dim==1, flux(end,:,:) = 0;
        elseif dim==2, flux(:,end,:) = 0;
        elseif dim==3, flux(:,:,end) = 0;
        end
    end
end

end

function [fm, fp, fmm, fpp, fppp] = makestencil (f, dim, BC)
% 
% makes stencil for calculate differences
% use circshift which is faster than slicing
% 
% fmm  (i-2)
% fm   (i-1)
% fcc  ( i )
% fp   (i+1)
% fpp  (i+2)

shift = circshift([1, 0, 0], [0, dim-1]);
sten5 = (nargout>3) ;

% 3 point stencil
fm = circshift(f,  shift);
fp = circshift(f, -shift);

if (sten5) 
    % extras for 5 point stencil
    fmm  = circshift(f,  2*shift);
    fpp  = circshift(f, -2*shift);
    fppp = circshift(f, -3*shift);
end

% if the boundary is closed (could use more elegance) but i wanted to be
% able to handle any specified dimension
if strcmp(BC,'closed')
    if dim==1
        fm( 1 ,:,:) = f( 1 ,:,:);
        fp(end,:,:) = f(end,:,:);
        
        if (sten5)
            fmm (    1:2  ,:,:) = f([  1,1  ],:,:);
            fpp (end-1:end,:,:) = f([end,end],:,:);
            fppp(end-2:end,:,:) = f(end-2:end,:,:);
        end
    elseif dim==2
        fm(:, 1 ,:) = f(:, 1 ,:);
        fp(:,end,:) = f(:,end,:);
        
        if (sten5)
            fmm (:,    1:2  ,:) = f(:,[  1,1  ],:);
            fpp (:,end-1:end,:) = f(:,[end,end],:);
            fppp(:,end-2:end,:) = f(:,end-2:end,:);
        end
    elseif dim==3
        fm(:,:, 1 ) = f(:,:, 1 );
        fp(:,:,end) = f(:,:,end);
        
        if (sten5)
            fmm (:,:,    1:2  ) = f(:,:,[  1,1  ]);
            fpp (:,:,end-1:end) = f(:,:,[end,end]);
            fppp(:,:,end-2:end) = f(:,:,end-2:end);
        end
    end
    
end
end

%% weno5 functions
function [dflux] = weno5 (f, v, dim, BC)
% calculate the flux difference through the cells

% define shifting dimension 
shift = circshift([1, 0, 0], [0, dim-1]);

% % Lax-Friedrichs flux splitting to ensure upwind-bias
% amp  = max(abs(v));
% fpos = 0.5*(f.*v + f.*amp);
% fneg = shiftflux(0.5*(f.*v - f.*amp), -shift, dim, BC);
% 
% % get positive fluxes
% [fm, fp, fmm, fpp] = makestencil(fpos, dim, BC);
% fluxpos = makeweno5poly(fmm, fm, fpos, fp, fpp);
% 
% % get negative fluxes
% % just change the upwind direction to use the same function
% [fm, fp, fmm, fpp] = makestencil(fneg, dim, BC);
% fluxneg = makeweno5poly(fpp, fp, fneg, fm, fmm);
% 
% dflux = fluxpos - shiftflux(fluxpos, shift, dim, BC) + ...
%         fluxneg - shiftflux(fluxneg, shift, dim, BC);

vpos = 0.5*(v + abs(v));    % positive velocity
vneg = 0.5*(v - abs(v));    % negative velocity

[fm, fp, fmm, fpp] = makestencil(f, dim, BC);
fpos    = makeweno5poly(fmm, fm, f, fp, fpp);
fluxpos = vpos.*fpos;

[fm, fp, ~, fpp, fppp] = makestencil(f, dim, BC);
fneg    = makeweno5poly(fppp, fpp, fp, f, fm);
fluxneg = vneg.*fneg;

dflux = fluxpos - shiftflux(fluxpos, shift, dim, BC) + ...
        fluxneg - shiftflux(fluxneg, shift, dim, BC);
end

function [fhalf] = makeweno5poly (fww, fw, fc, fe, fee)
% 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics

% 5th order polynomials
p1 = (2*fww - 7*fw + 11*fc )/6;
p2 = ( -fw  + 5*fc +  2*fe )/6;
p3 = (2*fc  + 5*fe -    fee)/6;

% smoothness measure (special from this paper)
b1 = 13/12*(fww - 2*fw + fc ).^2 + 1/4*(  fww - 4*fw + 3*fc ).^2; 
b2 = 13/12*(fw  - 2*fc + fe ).^2 + 1/4*(  fw  -          fe ).^2;
b3 = 13/12*(fc  - 2*fe + fee).^2 + 1/4*(3*fc  - 4*fe +   fee).^2;

% weights
g   = [1/10, 6/10, 3/10]; 
eps = 1e-16;    %stabiliser
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);
wp3 = g(3)./(b3.^2 + eps);

% get flux (normalize weights at the same time)
fhalf = (wp1.*p1 + wp2.*p2 + wp3.*p3) ./ (wp1 + wp2 + wp3) ;

end

%% tvd function

function [dflux] = tvd (f, v, dim, BC)
% total variation diminishing approach as described in Sramek 2010
% flux is split to ensure upwind-bias.
% vpos for positive fluxes; vneg for negative fluxes

% define shifting dimension 
shift = circshift([1, 0, 0], [0, dim-1]);

[fm, fp, ~, fpp] = makestencil(f, dim, BC);

vpos = 0.5*(v + abs(v));    % positive velocity
vneg = 0.5*(v - abs(v));    % negative velocity

Rpos = ( shiftflux(vpos, shift,dim,BC).*(f  -fm))./( vpos.*(fp-f) );
Rneg = ( shiftflux(vneg,-shift,dim,BC).*(fpp-fp))./( vneg.*(fp-f) );

% minmod approach
% lpos = max(0, min(1,Rpos));
% lneg = max(0, min(1,Rneg));

% superbee scheme
lpos = max( zeros(size(Rpos)), max(min(1,2*Rpos), min(2, Rpos)) );
lneg = max( zeros(size(Rneg)), max(min(1,2*Rneg), min(2, Rneg)) );

fpos = f  + 0.5*lpos.*(fp - f );
fneg = fp + 0.5*lneg.*(f  - fp);

Fpos = vpos.*fpos;
Fneg = vneg.*fneg;

dflux = Fpos - shiftflux(Fpos, shift, dim, BC) + ...
        Fneg - shiftflux(Fneg, shift, dim, BC);
end

