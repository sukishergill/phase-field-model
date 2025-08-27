% Simulation benchmark problem

Para.epsilon = 0.18;
Para.alpha = 1.7;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;
Para.d = 1.3;
Para.m = 0;
Para.q = 0;

dim = 2;
L = [2*pi, 2*pi];
N = [2^7, 2^7];
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

Para.B = Grid.L(1)^5;
Para.S = 3 * Para.epsilon^2;

Time.dt_min = 0.01;
Time.dt_max = 0.01;
Time.dt = 0.001;
Time.tf = 20;

% Initial condition

u = 2*exp(sin(Grid.xx-pi) + sin(Grid.yy-pi) - 2) + ...
    2.2*exp(-sin(Grid.xx-pi) - sin(Grid.yy-pi) - 2) - 1;

Results = FCH_BDF2_SAV(Grid, Time, Para, u);

%%

figure(1)
plot(Results.tt(2:end), Results.Eu(2:end), 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E_{FCH}$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

figure(2)
subplot(2,2,1)
surf(Results.uu{1}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['t = ', num2str(Results.tt(1))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,2)
surf(Results.uu{21}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['t = ', num2str(Results.tt(21))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,3)
surf(Results.uu{41}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['t = ', num2str(Results.tt(41))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,4)
surf(Results.uu{101}), shading interp, lighting phong, axis equal,axis off;
view(2);
title(['t = ', num2str(Results.tt(101))], 'Interpreter', 'latex')
c = colorbar;
set(c, 'TickLabelInterpreter','latex')
set(c, 'Position', [0.9 0.168 0.022 0.7]);
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')