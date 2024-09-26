% run FCH BDF2 SAV

epsilon = 0.18;
m = 0.2;
Lx = 2*pi;                              Ly = 2*pi;
nx = 256;                               ny = 256;                    
x = Lx*(1:nx)' / nx - Lx/2;             y = Ly*(1:ny)' / ny - Ly/2;
dt = 0.001;
tf = 10;

q = 0;
alpha = 1.7;
eta = epsilon^2;
tau = 0.125;
B = Lx^5;
S = 3*epsilon^2;


[xx,yy] = meshgrid(x,y);

% Initial condition

u = 2*exp(sin(xx) + sin(yy) - 2) + 2.2*exp(-sin(xx) - sin(yy) - 2) - 1;

[ xx, yy, k, tt, uu, Eu, Em, mass] = ...
    FCH_BDF2_SAV_v3 ( nx, ny, Lx, Ly, dt, tf, epsilon, eta, B, S, u);

% [ xx, yy, k, tt, uu, Eu, mass] = FCH_BDF2_SAV_v2 ( nx, ny,...
%     Lx, Ly, dt, tf, epsilon, eta1, eta2, tau, B, S, u);


figure(1)
plot(tt, Eu - Em, 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E_{FCH}$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

figure(2)
subplot(2,2,1)
surf(x,y,uu{26}), shading interp, lighting phong, axis equal,axis off;
view([-90 90]);
title(['t = ', num2str(tt(26))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,2)
surf(x,y,uu{51}), shading interp, lighting phong, axis equal,axis off;
view([-90 90]);
title(['t = ', num2str(tt(51))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,3)
surf(x,y,uu{76}), shading interp, lighting phong, axis equal,axis off;
view([-90 90]);
title(['t = ', num2str(tt(76))], 'Interpreter', 'latex')
% c = colorbar;
% set(c, 'TickLabelInterpreter','latex')
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,4)
surf(x,y,uu{101}), shading interp, lighting phong, axis equal,axis off;
view([-90 90]);
title(['t = ', num2str(tt(101))], 'Interpreter', 'latex')
c = colorbar;
set(c, 'TickLabelInterpreter','latex')
set(c, 'Position', [0.9 0.168 0.022 0.7]);
caxis([-1 , 1])
set(gca, 'Fontsize', 40)
set(gca,'TickLabelInterpreter','latex')
colormap('gray')