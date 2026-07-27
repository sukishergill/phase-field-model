% Allen-Cahn example 1

% clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Model set up %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 3;

% model parameters
Para.alpha = 0;
Para.beta = 1;
Para.epsilon = 0.1;
Para.err_tol = 10E-5;
Para.rho_s = 0.9;
Para.sigma = 500;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 1;

% spatial discretization
L = [2*pi, 2*pi];
N = [128, 128];
dim = 2;
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 0.01;        % minimum time step
Time.dt_max = 0.01;        % maximum time step
Time.tf = 8;
Time.adap = 3;
Time.t0 = 0;

% Initial condition
u = 0.05*sin(Grid.xx).*sin(Grid.yy);
% u = 0.05 * ones(Grid.N(1), Grid.N(2));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 8);

% plot_SSAV
% plot_err