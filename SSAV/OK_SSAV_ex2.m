% Ohta-Kawasaki example 2

% clear;

% 8.7137e-0.6

model = 1;

Para.alpha = 0;
Para.beta = 1;
Para.epsilon = 0.06;
Para.err_tol = [1e-4, 1e-7, 1e-6, 1e-8, 1e-7, 1e-6];
Para.err_tol_AM3 = 1e-5;
% Para.err_tol_MEE = 1e-5;
% Para.err_tol_AMEE = 1e-8;
Para.rho_s = 0.9;
Para.sigma = 50;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 1;
Para.delta = 1e-8;

dim = 2;
L = [2*pi, 2*pi];
N = [128, 128];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

% % spatial discretization
% Grid.Lx = 2*pi;                              Grid.Ly = 2*pi;
% Grid.Nx = 128;                               Grid.Ny = 128;     
% 
% Grid.dx = Grid.Lx/Grid.Nx;      Grid.dy = Grid.Ly/Grid.Ny;
% 
% x = Grid.Lx*(1:Grid.Nx)' / Grid.Nx - Grid.Lx/2;             
% y = Grid.Ly*(1:Grid.Ny)' / Grid.Ny - Grid.Ly/2;
% 
% [xx,yy] = meshgrid(x,y);
% 
% % spectral discretization
% kx = [ 0:Grid.Nx/2-1, 0.0, -Grid.Nx/2+1:-1]' / (Grid.Lx/pi/2);
% ky = [ 0:Grid.Ny/2-1, 0.0, -Grid.Ny/2+1:-1]' / (Grid.Ly/pi/2);
% [kxx,kyy] = meshgrid(kx, ky);
% k = sqrt(kxx.^2 + kyy.^2);
% Grid.kxx = kxx;         Grid.kyy = kyy;
% Grid.k = k;
% 
% Grid.inv_k = 1./(Grid.k.^2);      Grid.inv_k(k == 0) = 1;

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 1E-5;        % minimum time step
Time.dt_max = 1E-1;        % maximum time step
Time.t0 = 0;
Time.tf = 30;
Time.adap = 5;

% Initial condition
u = 0.25*(-tanh(sqrt( (Grid.xx + 0.8).^2 + (Grid.yy).^2 - 1.4) / (1.5*Para.epsilon)) - ...
    tanh(sqrt( (Grid.xx - 1.7).^2 + (Grid.yy).^2 - 0.5) / (1.5*Para.epsilon))) + 0.75;
% u = 0.25*(-tanh(sqrt( (Grid.xx + 1.5).^2 + (Grid.yy).^2 - 1.5) / (1.5*Para.epsilon)) - ...
% tanh(sqrt( (Grid.xx - 2).^2 + (Grid.yy).^2 - 0.5) / (1.5*Para.epsilon))) + 0.75;
u = real(u);
u = (u-0.375)/0.125;

Para.m = mean(u(:));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 20);

% plot_SSAV
% plot_err