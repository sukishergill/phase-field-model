% Ohta-Kawasaki example 1

% clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Model set up %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 1;

% model parameters
Para.alpha = 250000;
Para.OK = 1;
Para.beta = 1;
Para.epsilon = 0.02;
Para.err_tol = [1e-4, 1e-4, 1e-6, 1e-8, 1e-7, 1e-6];
Para.err_tol_MEE = 1e-4;
Para.err_tol_AMEE = 1e-9;
Para.rho_s = 0.9;
Para.sigma = 500;
Para.m = 0.1;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 1;
Para.delta = 1e-8;

% spatial discretization
% 
dim = 2;
L = [2*pi, 2*pi];
N = [128, 128];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 1e-5;        % minimum time step
Time.dt_max = 0.01;        % maximum time step
Time.t0 = 0;
Time.tf = 5;
Time.adap = 5;

% Initial condition
% u = 0.001*rand(Grid.N(1), Grid.N(2));
u = 0.002*rand(Grid.N(1), Grid.N(2))-0.001;
% u = smoothdata2(u);
% u = Para.m + u_rand;
% u = Results_spot.uu;
% u = Results_01.uu;

% u = 0.01*sin(Grid.xx).*sin(Grid.yy);

u = Para.m + u - mean(u(:));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 50);

% [u, Eu, Em, mass, t_vals, rel_err_vals] = ...
%     SSAV_2D_AM(Grid, Time, Para, u, model);

% plot_SSAV
% plot_err