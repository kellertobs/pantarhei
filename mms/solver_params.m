


outdir = '../out/';         % directory to save output files
nop    = 1;               % plot and store output every [nop] time step
svop   = 1;                 % save output
restart= 0;

N      = 100;               % number of grid points in each direction
D      = 20;                % domain dimension in each direction [m]
h      = D/N;               % grid spacing [m]
BC     = 'periodic';        % boundary conditions: 'open', 'closed', 'periodic'
tend   = 1;                 % model run time [s]
NtMax  = 0;                 % max number of time steps

nupd   = 100;               % update residual and permissions every [nupd] iterations
atol   = 1e-5;              % absolute residual tolerance for convergence of iterative solver
rtol   = 1e-10;             % relative residual tolerance for convergence of iterative solver
minits = 500;               % minimum iteration count for iterative solver
maxits = 1e5;               % maximum iteration count for iterative solver
alpha  = 0.90;              % first-order iterative step size (reduce if not converging)
beta   = 0.60;            % second-order iterative step size (reduce if not converging)
cfl    = 0.99;              % Courant number to limit physical time step size
flim   = 1e-16;             % limit phase fractions in coefficient closures
thtlim = 1e+16;             % limit phase-internal permission contrasts
cfflim = 1e+16;             % limit inter-phase coefficient contrasts

grav = [-9.81,0];           % gravity in vertical and horizontal direction
Gmg  = [1;-1].*0e-4;      % impose gaussian-shaped mass transfer rate (unity sum!)
dfg  = dfr;               % initial guassian peak amplitude (unity sum!)