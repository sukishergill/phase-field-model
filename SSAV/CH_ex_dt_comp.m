% Cahn-Hilliard example 1

% clear all;

model = 1;

Para.alpha = 0;
Para.OK = 0;
Para.beta = 1;
Para.epsilon = 0.1;
tau = 1e-5;
Para.err_tol = [tau, tau, 1e-5, 1e-6, 1e-9, 1e-4];
% Para.err_tol_AM3 = 1e-5;
% Para.err_tol = 1e-5;
% Para.err_tol_MEE = 1e-4;
% Para.err_tol_AMEE = 1e-8;
Para.rho_s = 0.9;
Para.sigma = 10000;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 3;

Para.delta = 1e-8;

dim = 2;
L = [2*pi, 2*pi];
N = [128, 128];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 1E-5;        % minimum time step
Time.dt_max = 1E-5;        % maximum time step
Time.t0 = 0;    % starting time
Time.tf = 50;    % final time

% Choice of adaptive time stepper
%       - Third order Adams-Moulton: 'AM3'
%       - Variation in E_t: 'Et'
Time.adap = 5;

save_u = 1;

% Initial condition
cos_IC;
u = Para.m + u - mean(u(:));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, Time.tf);
% Results = SSAV_2D_AM (Grid, Time, Para, u, model);