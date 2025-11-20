% Benchmark simulation

Para.epsilon = 0.1;
Para.eta1 = 1.45 * Para.epsilon;
Para.eta2 = 3 * Para.epsilon;
Para.err_tol = [1e-4, 1e-4];
Para.tau = 0.3;
Para.d = 0.2;
Para.m = 0;
Para.q = 0;
Para.W = 1;
Para.zeta = 1;
Para.delta = 1e-8;

% Para.alpha = 0;

if Para.q == 0
    Para.alpha = 1.7;

elseif Para.q == 0.2
    Para.alpha = 5.1;

elseif Para.q == 0.5
    Para.alpha = 10.2;

end

dim = 2;
L = [4*pi, 4*pi];
N = [2^8, 2^8];
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

Para.B = Grid.L(1)^5;
Para.S = 3 * Para.epsilon^2;
Para.S = 0;

Time.dt_min = 1e-5;
Time.dt_max = 0.001;
Time.dt = 0.001;
Time.tf = 15;

% Initial condition

load('FCH_benchmark.mat', 'u');
u = u(1:4:end, 1:4:end);
u = u + Para.epsilon*Para.d / (Para.alpha^2);

Para.m = 0.5*sum(u+1, 'all')*prod(Grid.d);
% Para.m = 6.11;

Results = FCH_BDF2_SAV_v2(Grid, Time, Para, u, 15);

%%

figure;
subplot(2,1,1)
plot(Results.t_vals(2:end), Results.Eu(2:end)/prod(L), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E_{FCH}$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,1,2)
semilogy(Results.t_vals(2:end), (Results.t_vals(2:end) - Results.t_vals(1:end-1)), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$\delta t$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

figure;
surf(Results.uu{end}), shading interp, lighting phong, axis equal,axis off;
view(2);
% subplot(2,2,1)
% surf(Results.uu{1}), shading interp, lighting phong, axis equal,axis off;
% view(2);
% title(['t = ', num2str(Results.tt(1))], 'Interpreter', 'latex')
% % c = colorbar;
% % set(c, 'TickLabelInterpreter','latex')
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% subplot(2,2,2)
% surf(Results.uu{21}), shading interp, lighting phong, axis equal,axis off;
% view(2);
% title(['t = ', num2str(Results.tt(21))], 'Interpreter', 'latex')
% % c = colorbar;
% % set(c, 'TickLabelInterpreter','latex')
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% subplot(2,2,3)
% surf(Results.uu{41}), shading interp, lighting phong, axis equal,axis off;
% view(2);
% title(['t = ', num2str(Results.tt(41))], 'Interpreter', 'latex')
% % c = colorbar;
% % set(c, 'TickLabelInterpreter','latex')
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% subplot(2,2,4)
% surf(Results.uu{101}), shading interp, lighting phong, axis equal,axis off;
% view(2);
% title(['t = ', num2str(Results.tt(101))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
% set(c, 'Position', [0.9 0.168 0.022 0.7]);
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')

