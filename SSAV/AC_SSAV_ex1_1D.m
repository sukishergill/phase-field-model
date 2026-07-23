% Allen-Cahn 1D example

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
tau = 1e-4;
Para.err_tol = [tau, tau, 1e-6, 1e-8, 1e-7, 1e-6];
Para.err_tol_AM3 = 1e-3;
Para.rho_s = 0.9;
Para.sigma = 500;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 1;
Para.delta = 1e-8;

% spatial discretization
dim = 1;
L = 2*pi;
N = 128;
Grid = SSAV_helpers.generate_Grid(L, N, dim);

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 1e-5;        % minimum time step
Time.dt_max = 0.1;        % maximum time step
Time.dt_max = 1;
Time.t0 = 0;
Time.tf = 100;
Time.adap = 5;

% Initial condition
u = tanh((Grid.xx + pi*(Para.m + 1)/2) / Para.epsilon) - ...
    tanh((Grid.xx - pi*(Para.m + 1)/2) / Para.epsilon) - 1;
% u = 0.001*rand(size(Grid.xx));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 10);

% plot_SSAV
% plot_err