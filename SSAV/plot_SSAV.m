% plotting

figure;
subplot(2,1,1)
plot(Results.t_vals(2:end), Results.Eu(2:end), 'Linewidth', 3)
xlabel('Time', Interpreter='latex')
title('$E(u)$', Interpreter='latex')
set(gca, 'FontSize', 30)
set(gca, 'TickLabelInterpreter', 'latex')
subplot(2,1,2)
semilogy(Results.t_vals(2:end), Results.dt_vals(1:end-1), 'Linewidth', 3)
xlabel('Time', Interpreter='latex')
title('$\delta t$', Interpreter='latex')
set(gca, 'FontSize', 30)
set(gca, 'TickLabelInterpreter', 'latex')

% figure(1)
% plot(tt, Eu - Em, 'Linewidth', 4)
% set(gca, 'Fontsize', 40)
% xlabel('Time', 'Interpreter','latex')
% title('$E(u,w)-E(m)$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% 
% figure;
% 
% subplot(2,2,1)
% surf(uu{2}), shading interp, lighting phong, axis equal,axis off;
% view(2);
% title(['t = ', num2str(uu_times(2))], 'Interpreter', 'latex')
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% subplot(2,2,2)
% surf(uu{51}), shading interp, lighting phong, axis equal,axis off;
% view([-90 90]);
% title(['t = ', num2str(uu_times(51))], 'Interpreter', 'latex')
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% subplot(2,2,3)
% surf(uu{76}), shading interp, lighting phong, axis equal,axis off;
% view([-90 90]);
% title(['t = ', num2str(uu_times(76))], 'Interpreter', 'latex')
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% subplot(2,2,4)
% surf(uu{101}), shading interp, lighting phong, axis equal,axis off;
% view([-90 90]);
% title(['t = ', num2str(uu_times(101))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
% set(c, 'Position', [0.9 0.168 0.022 0.7]);
% caxis([-1 , 1])
% set(gca, 'Fontsize', 40)
% set(gca,'TickLabelInterpreter','latex')
% colormap('gray')
