% Cahn-Hilliard example 1

% clear all;

model = 1;

Para.alpha = 0;
Para.beta = 1;
Para.epsilon = 1;
% tau_1 = 1e-1/(2^12);
% tau_1 = 0.001;
Para.err_tol = [tau_1, tau_1, 1e-5, 1e-6, 1e-9, 1e-4];
% Para.err_tol_AM3 = 1e-8;
% Para.err_tol = 1e-5;
% Para.err_tol_MEE = 1e-4;
% Para.err_tol_AMEE = 1e-8;
Para.rho_s = 0.9;
Para.sigma = 2000;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.p = 3;

Para.delta = 1e-8;

dim = 2;
L = [2*pi, 2*pi];
N = [64, 64];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

dt = 0.01;
% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = tau_1 / 100;        % minimum time step
Time.dt_max = tau_1 * 100;        % maximum time step

% Time.dt_min = 1e-6;        % minimum time step
% Time.dt_max = 0.1;        % maximum time step
Time.t0 = 0;    % starting time
Time.tf = 1;    % final time

u_exact = sin(Grid.xx).*cos(Grid.yy).*cos(Time.tf);
% u_exact = cos(Grid.xx).*cos(Grid.yy).*exp(-Time.tf);

% Choice of adaptive time stepper
%       - Third order Adams-Moulton: 'AM3'
%       - Variation in E_t: 'Et'
Time.adap = 5;

save_u = 1;

% Initial condition
% u = cos(Grid.xx).*cos(Grid.yy);
u = sin(Grid.xx).*cos(Grid.yy);
% u = 0.25 + 0.4*rand(Grid.Nx, Grid.Ny);
% Para.m = mean(u(:));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 1);