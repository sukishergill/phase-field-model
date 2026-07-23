% Simulation benchmark problem

Para.epsilon = 0.1;
Para.alpha = 0;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.err_tol = [1e-4, 1e-4];
Para.tau = 0;
Para.d = 1.3;
Para.m = 0;
Para.q = 0;
Para.zeta = 1;
Para.delta = 1e-8;

dim = 1;
% L = [4*pi, 4*pi];
L = 2*pi;
% L = [L, L];
% N = [2^7, 2^7];
N = 2^7;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

Para.B = Grid.L(1)^5;
Para.S = 3 * Para.epsilon^2;
% Para.S = 1;

Time.dt_min = 1E-5;
Time.dt_max = 0.001;
Time.dt = 0.001;
Time.t0 = 0;
Time.tf = 100;

% Initial condition
% m = 0.9;
% u = Para.m + 0.002*rand(Grid.N(1), Grid.N(2))-0.001;
% u = Para.m + 0.02*rand(Grid.N(1),1)-0.01;
% u = smoothdata2(u);
% u = 0.05*sin(Grid.xx).*sin(Grid.yy);
% u = 0.05*sin(Grid.x) + Para.m;

% a = [0.3 0.3 0.4];
% u = a(1)*cos(sqrt(2)*Grid.xx/Grid.L(1)) + a(2)*cos(-sqrt(1/2)*Grid.xx/Grid.L(1) + ...
%     sqrt(3/2)*Grid.yy / Grid.L(2)) + a(3)*cos(-sqrt(1/2)*Grid.xx/Grid.L(1) - ...
%     sqrt(3/2)*Grid.yy / Grid.L(2));

% Para.m = mean(u, 'all');
% u = Para.m * ones(size(Grid.x));

u = tanh((Grid.x + pi*(Para.m + 1)/2) / Para.epsilon) - ...
    tanh((Grid.x - pi*(Para.m + 1)/2) / Para.epsilon) - 1;

Results = FCH_BDF2_SAV_adap_dt(Grid, Time, Para, u, Time.tf/10);
Results.dt_vals = [Results.t_vals(2:end) - Results.t_vals(1:end-1) 0];
% Results.a = a;
%%

% figure;
% plot(Results.t_vals(2:end), Results.Eu(2:end), 'Linewidth', 4)
% set(gca, 'Fontsize', 40)
% xlabel('Time', 'Interpreter','latex')
% title('$E_{FCH}$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% 
% figure;
% subplot(2,1,1)
% plot(Results.t_vals(2:end), Results.Eu(2:end), 'Linewidth', 4)
% set(gca, 'Fontsize', 40)
% xlabel('Time', 'Interpreter','latex')
% title('$E_{FCH}$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% xlim([0 Time.tf])
% subplot(2,1,2)
% semilogy(Results.t_vals(2:end), (Results.t_vals(2:end) - Results.t_vals(1:end-1)), 'Linewidth', 4)
% set(gca, 'Fontsize', 40)
% xlabel('Time', 'Interpreter','latex')
% title('$\delta t$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
%%
% if dim == 1
    % figure;
    % plot(Grid.x, Results.uu{end}, 'Linewidth', 3)
% 
% elseif dim == 2
%     figure;
%     surf(Results.uu{end}), shading interp, lighting phong, axis equal,axis off;
%     % title(['$m = $', num2str(Para.m)], 'Interpreter', 'latex')
%     view(2);
%     % colormap(c);
%     % clim([-1 1])
% end