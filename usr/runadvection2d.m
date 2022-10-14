% script to test the different advection schemes, for a 2D problem
% YQW, 11 October 2022

clear all;
addpath('../src/');

%% options to change

D   = 10;     % domain size
N   = 51;     % number of grid cells
cfl = 0.4;    % courant flow number limiter

% constant velocities (can try both + and - velocities)
u0 = 1;
w0 = 0;

% initial condition type. options: {'square', 'gaussian'}
ftype = 'square';
BC    = 'closed';

% advection scheme 
% options: {'centr', 'upwd1', 'quick', 'fromm', 'weno3', 'weno5', 'tvdim'}
schm  = {'fromm','weno5','tvdim'};


%% initialize fields

Nschm= length(schm);

% space vectors
x   = linspace(-D,D,N);
z   = x';
dx  = x(2) - x(1);

% velocity field
u    = u0*ones(N,N+1);   
w    = w0*ones(N+1,N);
vmax = max(abs([u(:);w(:)]));

% time step
dt   = 0.5*cfl*dx/vmax;
if strcmp(BC, 'periodic'), tEnd = (2*D+dx)/vmax;
else, tEnd = (0.5*D+dx)/vmax;
end

% dimensions and BCs corresponding to z, x
dim = [1,2];
BC  = repmat({BC},1,2);

% initialize field
switch ftype
    case 'gaussian'
        sig = D/5;
        f0 = exp(-(x./sig).^2 - (z./sig).^2);
    case 'square'
        f0 = zeros(N,N) + 1*(x>-D/4 & x<D/4 & z>-D/4 & z<D/4);
        f0 = movmean(movmean(f0,3,1),3,2);
end

%% plot initial condition

figure; 
set(gcf,'Position',[500,400,400,700]);
tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
plotfield(x, f0, f0, [1,5,7], 'initial condition');


%% run advection

fmat = repmat(f0, 1, 1, Nschm);

f1 = figure; 
Nrow = 4; 
Ncol = 2*Nschm; 
pltvec = [0,2*Ncol,3*Ncol];

set(gcf,'Position',[200,500,Ncol/2*200,Nrow*120]);
tiledlayout(Nrow,Ncol,'TileSpacing','tight','Padding','tight');

% plot initial condition
for si = 1:Nschm
    plotfield(x, fmat(:,:,si), f0,  (2*si-1)+pltvec, schm{si});
end

t    = 0;
step = 0;

while t < tEnd
    t    = t + dt;
    step = step + 1;
      
    for si = 1:Nschm
        % loop through the different schemes
        fmat(:,:,si) = rk4(@(f) -1*advect(f,u,w,dx,{schm{si},''},dim,BC), fmat(:,:,si), dt);
    end
    
    % plot
    if mod(step,10)==0 || (t > tEnd)
        figure(f1);
        for si = 1:Nschm
            plotfield(x, fmat(:,:,si), f0,  (2*si-1)+pltvec, schm{si});
        end
        drawnow;
    end
end

%% time these different methods

fup = f0; fcd  = f0; fqk  = f0; ffm  = f0; fwn = f0; ftv = f0;
Niter = 2000;

tup = evaltime(@(f) -1*advect(f, u, w, dx, {'upwd1',''}, dim, BC), f0, dt, Niter);
tcd = evaltime(@(f) -1*advect(f, u, w, dx, {'centr',''}, dim, BC), f0, dt, Niter);
tqk = evaltime(@(f) -1*advect(f, u, w, dx, {'quick',''}, dim, BC), f0, dt, Niter);
tfm = evaltime(@(f) -1*advect(f, u, w, dx, {'fromm',''}, dim, BC), f0, dt, Niter);
tw3 = evaltime(@(f) -1*advect(f, u, w, dx, {'weno3',''}, dim, BC), f0, dt, Niter);
tw5 = evaltime(@(f) -1*advect(f, u, w, dx, {'weno5',''}, dim, BC), f0, dt, Niter);
ttv = evaltime(@(f) -1*advect(f, u, w, dx, {'tvdim',''}, dim, BC), f0, dt, Niter);

fprintf('\n\n time taken [seconds] by advection schemes (Niter = %d)\n\n', Niter);
fprintf('\t upwind  \t %.4f\n', tup); 
fprintf('\t central \t %.4f\n', tcd); 
fprintf('\t quick   \t %.4f\n', tqk); 
fprintf('\t fromm   \t %.4f\n', tfm); 
fprintf('\t weno3   \t %.4f\n', tw3); 
fprintf('\t weno5   \t %.4f\n', tw5); 
fprintf('\t tvdim   \t %.4f\n', ttv); 
fprintf('\n\n\n');

%% functions

function [tout] = evaltime (xfunc, f, dt, Nstep)

tic
for i = 1:Nstep
    f = rk4(xfunc, f, dt);
end
tout = toc;

end

function [f] = rk4 (xfunc, f, dt)
% runge kutta forward integration in time
k1 = xfunc(f);
k2 = xfunc(f + dt*k1/2);
k3 = xfunc(f + dt*k2/2);
k4 = xfunc(f + dt*k3);
f  = f + 1/6*dt*(k1 + 2*k2 + 2*k3 + k4);
end

function [] = plotfield (x, f, f0, iax, titletext)

nexttile(iax(1),[2,2]); 
imagesc(x,x,squeeze(f)); axis image; set(gca,'ydir','normal');
colorbar('Location','southoutside'); caxis([0,1]);
title(titletext);

nexttile(iax(2),[1,2]); 
plot(x, squeeze(f0(ceil(length(x)/2),:)), ':' ); hold on;
plot(x, squeeze( f(ceil(length(x)/2),:)), '-+'); hold off;
ylim([-0.2,1.2]);
xlabel('x, horizontal');

nexttile(iax(3),[1,2]); 
plot(x, squeeze(f0(:,ceil(length(x)/2))), ':' ); hold on;
plot(x, squeeze( f(:,ceil(length(x)/2))), '-+'); hold off;
ylim([-0.2,1.2]);
xlabel('z, vertical');
end