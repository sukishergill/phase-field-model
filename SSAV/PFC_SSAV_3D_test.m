% Phase-field crystal example 1

% clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Model set up %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 2;

Para.alpha = 0;
Para.beta = 0.2;
Para.epsilon = 1;
% Para.err_tol = 10E-2;
Para.rho_s = 0.9;
Para.sigma = 500;
Para.m = 0.25;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 1;
Para.delta = 1e-8;
Para.err_tol = [1e-3, 1e-3, 1e-5, 1e-6, 1e-9, 1e-4];

% spatial discretization
dim = 3;
L = [8*pi, 8*pi, 8*pi];
% L = [128, 128, 128];
N = [64, 64, 64];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_max = 10;        % minimum time step
Time.dt_min = 0.1;        % maximum time step
Time.t0 = 0;
Time.tf = 1000;
Time.adap = 5;

% Initial condition
% u = 0.2*cos(pi*xx / 16).*cos(pi*yy / 16);
% u = 0.2*rand(Grid.Nx,Grid.Ny);
% u = Para.m + u - mean(u(:));

% u = 0.01*rand(Grid.Nx,Grid.Ny);
% u = Para.m + u;
% u = 0.01*rand(Grid.N(1), Grid.N(2), Grid.N(3));
u = 0.1*rand(Grid.N(1), Grid.N(2), Grid.N(3));
u = Para.m + u - mean(u(:));
% u = cos(Grid.xx).*cos(2*Grid.yy).*cos(4*Grid.zz);

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 10);

% [u, Eu, Eu_SSAV, Em, mass, t_vals] = ...
%     SSAV_2D(Grid, Time, Para, u, model);

% [u, Eu, Em, mass, m_est_vals, t_vals, dt_vals] = ...
%     SSAV_2D(Grid, Time, Para, u, model);

% plot_SSAV
% plot_err