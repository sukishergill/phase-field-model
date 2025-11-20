% Cahn-Hilliard example 1

% clear all;

model = 1;

Para.alpha = 0;
Para.OK = 0;
Para.beta = 1;
Para.epsilon = 0.1;
Para.err_tol = [1e-3, 1e-6, 1e-5, 1e-6, 1e-9, 1e-4];
% Para.err_tol_AM3 = 1e-6;
% Para.err_tol = 1e-5;
% Para.err_tol_MEE = 1e-4;
% Para.err_tol_AMEE = 1e-8;
Para.rho_s = 0.9;
Para.sigma = 500;
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

% spatial discretization
% % Grid.Lx = 2*pi;                              Grid.Ly = 2*pi;
% Grid.L = [2*pi, 2*pi];
% Grid.N = [128, 128];
% % Grid.Nx = 128;                               Grid.Ny = 128;     
% 
% % Grid.dx = Grid.Lx/Grid.Nx;      Grid.dy = Grid.Ly/Grid.Ny;
% Grid.d = Grid.L ./ Grid.N;
% % Grid.dxdy = Grid.dx * Grid.dy;
% 
% 
% % for i = 1:dim
% %     x(:,i) = Grid.L(i)*(1:Grid.N(i))' / Grid.N(i) - Grid.L(i)/2;  
% %     k(:,i) = [ 0:Grid.N(i)/2-1, 0.0, -Grid.N(i)/2+1:-1]' / (Grid.L(i)/pi/2);
% % end
% 
% x = Grid.L(1)*(1:Grid.N(1))' / Grid.N(1) - Grid.L(1)/2;  
% y = Grid.L(2)*(1:Grid.N(2))' / Grid.N(2) - Grid.L(2)/2; 
% [xx, yy] = meshgrid(x, y);
% 
% % 
% % if dim == 2
% %     [xx(:,:,1), xx(:,:,2)] = meshgrid(x(:,1), x(:,2));
% %     [kxx(:,:,1), kxx(:,:,2)] = meshgrid(x(:,1), x(:,2));
% % end
% % 
% % Grid.k = zeros(size(xx(:,:,1)));
% % Grid.k2 = zeros(size(xx(:,:,1)));
% % 
% % for i = 1:dim
% %     Grid.k2 = Grid.k2 + kxx(:,:,i).^2;
% % end
% 
% % spectral discretization
% kx = [ 0:Grid.N(1)/2-1, 0.0, -Grid.N(1)/2+1:-1]' / (Grid.L(1)/pi/2);
% ky = [ 0:Grid.N(2)/2-1, 0.0, -Grid.N(2)/2+1:-1]' / (Grid.L(2)/pi/2);
% [kxx,kyy] = meshgrid(kx, ky);
% Grid.k = sqrt(kxx.^2 + kyy.^2);
% Grid.kxx = kxx;         Grid.kyy = kyy;
% 
% Grid.inv_k = 1./(Grid.k.^2);      Grid.inv_k(Grid.k == 0) = 1;

% time discretization
% if dt_min = dt_max BDF2 will be implemented, otherwise an adaptive time
% stepping scheme will be used
Time.dt_min = 1e-6;        % minimum time step
Time.dt_max = 0.01;        % maximum time step
Time.t0 = 0;    % starting time
Time.tf = 8;    % final time

% Choice of adaptive time stepper
%       - Third order Adams-Moulton: 'AM3'
%       - Variation in E_t: 'Et'
Time.adap = 5;

save_u = 1;

% Initial condition
u = 0.05*sin(Grid.xx).*sin(Grid.yy);
% u = 0.25 + 0.4*rand(Grid.Nx, Grid.Ny);
% Para.m = mean(u(:));

Results = SSAV(Grid, Time, Para, u, model, dim, 0, 8);