% Simulation benchmark problem

Para.epsilon = 0.18;
Para.alpha = 0;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.err_tol = [1e-4, 1e-4];
Para.tau = 0;
Para.d = 1.3;
% m_override is set by run_ex5_m_sweep.m when this script is called from
% its parfor loop; standalone runs use the default below.
if exist('m_override', 'var')
    Para.m = m_override;
else
    Para.m = 0;
end
Para.q = 0;
Para.zeta = 1;
Para.delta = 1e-8;

dim = 2;
L = [2*pi, 2*pi];
% L = 2*pi;
% L = [L, L];
N = [64, 64];
% N = 2^7;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

Para.B = Grid.L(1)^5;
% Para.S = 3 * Para.epsilon^2;
Para.S = 1;

Time.dt_min = 1e-5;
Time.dt_max = 0.01;
Time.dt = 0.001;
Time.t0 = 0;
Time.tf = 150;

% Initial condition
% m = 0.9;
rng(10);
u = Para.m + (2*Para.m+0.05)*rand(Grid.N(1), Grid.N(2));

% Symmetrize so u is even in x and in y. Grid points are x_j = L*j/N - L/2,
% so x -> -x maps index j to N-j (with j = N, i.e. x = L/2 = -L/2, fixed).
[ny, nx] = size(u);
ix = [nx-1:-1:1, nx];
iy = [ny-1:-1:1, ny];
u = (u + u(:, ix) + u(iy, :) + u(iy, ix)) / 4;

u = u - mean(u(:)) + Para.m;

Results = FCH_BDF2_SAV_adap_dt(Grid, Time, Para, u, Time.tf/10);
Results.dt_vals = [Results.t_vals(2:end) - Results.t_vals(1:end-1) 0];
%

%{
figure;
plot(Results.t_vals(2:end), Results.Eu(2:end), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E_{FCH}$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

figure;
subplot(2,1,1)
plot(Results.t_vals(2:end), Results.Eu(2:end), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E_{FCH}$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
xlim([0 Time.tf])
subplot(2,1,2)
semilogy(Results.t_vals(2:end), (Results.t_vals(2:end) - Results.t_vals(1:end-1)), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$\delta t$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
%
if dim == 1
    figure;
    plot(Grid.x, Results.uu{end}, 'Linewidth', 3)

elseif dim == 2
    figure;
    surf(Results.uu{end}), shading interp, lighting phong, axis equal,axis off;
    % title(['$m = $', num2str(Para.m)], 'Interpreter', 'latex')
    view(2);
    % colormap(c);
    % clim([-1 1])
end
%}