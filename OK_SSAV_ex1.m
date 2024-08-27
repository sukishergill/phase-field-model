Para.alpha = 250000;
Para.beta = 1;
Para.epsilon = 0.02;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.rho_s = 0.9;    % safety coeff


Lx = 2*pi;                              Ly = 2*pi;
nx = 256;                               ny = 256;                    
x = Lx*(1:nx)' / nx - Lx/2;             y = Ly*(1:ny)' / ny - Ly/2;
dt = 0.001;
tf = 5;


dt = 0.001;
tf = 5;

t_vals = linspace(dt, tf, round(tf / dt));


[xx,yy] = meshgrid(x,y);


% Initial condition
u = Para.m + 0.001*rand(nx, ny);


[k, tt, uu, Eu, Eu_SSAV, Em, mass, m_est_vals, t_vals, dt_vals] = ...
    SSAV_2D( nx, ny, Lx, Ly, dt, tf, Para, u, 1);


figure(1)
plot(tt, Eu - Em, 'Linewidth', 4)
set(gca, 'Fontsize', 40)
xlabel('Time', 'Interpreter','latex')
title('$E(u,w)-E(m)$', 'Interpreter','latex')
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