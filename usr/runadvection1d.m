% script to test the different advection schemes, for a 1D problem
% YQW, 11 October 2022

clear all;
addpath('../src/');

%% options to change

D   = 10;     % domain size
N   = 101;     % number of grid cells
cfl = 0.4;    % courant flow number limiter

% constant velocities (can try both + and - velocities)
u0 = 1;

% initial condition type. options: {'square', 'gaussian'}
ftype = 'triangle';

% advection scheme 
% options: {'centr', 'upwd1', 'quick', 'fromm', 'weno3', 'weno5', 'tvdim'}
schm  = {'centr', 'upwd1', 'quick', 'fromm', 'weno3', 'weno5', 'tvdim'};

%% initialize fields

Nschm = length(schm);

% space vectors
x   = linspace(-D,D,N)';
dx  = x(2) - x(1);

% velocity field - have to assign w 
u    = u0*ones(N+1,1);   
w    =  0*ones(N  ,2); 
vmax = max(abs([u(:);w(:)]));

% time step
dt   = 0.5*cfl*dx/vmax;

% dimensions and BCs corresponding to z, x
dim = [2,1];
BC  = {'periodic','periodic'};

% initialize field
f0 = zeros(size(x));
switch ftype
    case 'square'
        f0(x>-D/4 & x<D/4,:) = 1;
    case 'triangle'
        f0(x>-D/4 & x<0  ,:) =  1 + x(x>-D/4 & x<0  ,:)/(D/4);       % triangular wave
        f0(x>=0   & x<D/4,:) =  1 - x(x>=0   & x<D/4,:)/(D/4);  % triangular wave
end
f0   = movmean(movmean(f0,3),3);      % smooth

%% plot initial condition

figure; 
plot(x, f0, '+-'); 
title('initial condition');


%% run advection

fmat = repmat(f0, 1, 1, Nschm);

f1 = figure; 
Nrow = floor(sqrt(Nschm));
Ncol = ceil(Nschm/Nrow);

set(gcf,'Position',[200,500,Ncol*250,Nrow*200]);
tiledlayout(Nrow,Ncol,'TileSpacing','tight','Padding','tight');

% plot initial condition
for si = 1:Nschm
    nexttile(si); 
    plot(x, f0, ':'); 
    ylim([-0.2,1.2]); 
    title(schm{si});
end

t    = 0;
step = 0;

while t < (2*D+dx)/vmax
    t    = t + dt;
    step = step + 1;
      
    for si = 1:Nschm
        % loop through the different schemes
        fmat(:,:,si) = rk4(@(f) -1*advect(f,u,w,dx,{schm{si},''},dim,BC), fmat(:,:,si), dt);
    end
    
    % plot
    if mod(step,10)==0 || (t > (2*D+dx)/vmax)
        figure(f1);
        for si = 1:Nschm
            nexttile(si); 
            plot(x, f0, ':', x, fmat(:,:,si)); 
            ylim([-0.2,1.2]); 
            title(schm{si});
        end
        drawnow;
    end
end

%% time these different methods

fup = f0; fcd  = f0; fqk  = f0; ffm  = f0; fwn = f0; ftv = f0;
Niter = 10000;

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
