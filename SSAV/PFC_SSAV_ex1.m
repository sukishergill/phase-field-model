% Phase-field crystal example 1

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Model set up %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 2;

Para.alpha = 0;
Para.beta = 0.07;
Para.epsilon = 1;
Para.err_tol = 10E-4;
Para.m = 0.01;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty

% spatial discretization
Grid.Lx = 32*pi;        Grid.Ly = 32*pi;         % domain size
Grid.Nx = 128;          Grid.Ny = 128;           % # of grid points
     
Grid.dx = Grid.Lx/Grid.Nx;      Grid.dy = Grid.Ly/Grid.Ny;

x = Grid.Lx*(1:Grid.Nx)' / Grid.Nx - Grid.Lx/2;             
y = Grid.Ly*(1:Grid.Ny)' / Grid.Ny - Grid.Ly/2;

[xx,yy] = meshgrid(x,y);

% spectral discretization
kx = [ 0:Grid.Nx/2-1, 0.0, -Grid.Nx/2+1:-1]' / (Grid.Lx/pi/2);
ky = [ 0:Grid.Ny/2-1, 0.0, -Grid.Ny/2+1:-1]' / (Grid.Ly/pi/2);
[kxx,kyy] = meshgrid(kx, ky);
k = sqrt(kxx.^2 + kyy.^2);
Grid.kxx = kxx;         Grid.kyy = kyy;
Grid.k = k;

Grid.inv_k = 1./(Grid.k.^2);      Grid.inv_k(k == 0) = 1;

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_max = 0.5;        % minimum time step
Time.dt_min = 0.5;        % maximum time step
Time.tf = 20000;

% Initial condition
u = 0.2*rand(Grid.Nx,Grid.Ny);
u = Para.m + u - mean(u(:));


[u, Eu, Eu_SSAV, Em, mass, t_vals] = ...
    SSAV_2D(Grid, Time, Para, u, model);

% [u, Eu, Em, mass, m_est_vals, t_vals, dt_vals] = ...
%     SSAV_2D(Grid, Time, Para, u, model);

% plot_SSAV
% plot_err