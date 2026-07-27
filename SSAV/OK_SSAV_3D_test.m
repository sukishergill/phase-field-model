% Ohta-Kawasaki example 1

% clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Model set up %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 1;

dim = 3;

% model parameters
Para.alpha = 3.5^2;
Para.beta = 1;
Para.epsilon = 1/3.5;
Para.err_tol = [1e-3, 1e-3, 1e-5, 1e-6, 1e-9, 1e-4];
Para.rho_s = 0.9;
Para.sigma = 500;
Para.m = 0.01;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 1;
Para.delta = 1e-8;

dim = 3;

L = [4*pi, 4*pi, 4*pi];
N = [64, 64, 64];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 0.01;        % minimum time step
Time.dt_max = 0.01;        % maximum time step
Time.t0 = 0;
Time.tf = 100;
Time.adap = 5;

% Initial condition
% u = u_IC;
u = Para.m + 0.001*rand(Grid.N(1), Grid.N(2), Grid.N(3));
u = Para.m + u - mean(u, 'all');

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 10);

% [u, Eu, Em, mass, t_vals, rel_err_vals] = ...
%     SSAV_2D_AM(Grid, Time, Para, u, model);

% plot_SSAV
% plot_err