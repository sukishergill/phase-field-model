% Simulation benchmark problem

Para.epsilon = 0.18;
Para.alpha = 0;
Para.eta1 = Para.epsilon^2;
% Para.eta1 = 1;
Para.eta2 = Para.eta1;
Para.err_tol = [1e-4, 1e-4];
Para.tau = 0;
Para.d = 1.3;
Para.m = 0;
Para.q = 0;
Para.zeta = 1;
Para.delta = 1e-8;

dim = 2;
L = [2*pi, 2*pi];
N = [2^7, 2^7];
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

Para.B = Grid.L(1)^5;
% Para.S = 3 * Para.epsilon^2;
Para.S = 10;

Time.dt_min = 1E-5;
Time.dt_max = 0.01;
Time.dt = 0.01;
Time.t0 = 0;
Time.tf = 50;

% Initial condition

u = 2*exp(sin(Grid.xx-pi) + sin(Grid.yy-pi) - 2) + ...
    2.2*exp(-sin(Grid.xx-pi) - sin(Grid.yy-pi) - 2) - 1;

Results = FCH_BDF2_SAV_adap_dt(Grid, Time, Para, u, Time.tf);

%%

figure;
plot(Results.t_vals(2:end), Results.Eu(2:end), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E_{FCH}$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

t1 = sum(Results.t_vals<1) + 1;
figure;
subplot(1,2,1)
plot(Results.t_vals, Results.Eu, 'Linewidth', 4)
set(gca, 'Fontsize', 30)
xlabel('Time', 'Interpreter','latex')
ylabel('$E(u)$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
xlim([0 Time.tf])
subplot(1,2,2)
semilogy(Results.t_vals(2:end-1), (Results.t_vals(2:end-1) - Results.t_vals(1:end-2)), 'Linewidth', 4)
set(gca, 'Fontsize', 30)
xlabel('Time', 'Interpreter','latex')
ylabel('$\Delta t$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
xlim([0,Time.tf])
%%

% figure;
% surf(Results.uu{end}), shading interp, lighting phong, axis equal,axis off;
% view(2);
figure;
subplot(1,4,1)
surf(Results.uu{1}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['$t = 0$'], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(1,4,2)
surf(Results.uu{21}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['$t = $ ', num2str(Results.t_vals( sum(Results.t_vals < 1)+1))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(1,4,3)
surf(Results.uu{51}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['$t = $ ', num2str(Results.t_vals( sum(Results.t_vals < 5)+1))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(1,4,4)
surf(Results.uu{end}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['$t= 20$'], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
% set(c, 'Position', [0.9 0.168 0.022 0.7]);
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')