Para.alpha = 0;
Para.beta = 1;
Para.epsilon = 0.1;
Para.m = 0;
Para.M = 1;
Para.B = 1;          % const. that ensures positive radicand
Para.S = 2;          % positive stabilizing parameter S ~ ||f(u)||_\infty
Para.rho_s = 0.9;    % safety coeff


Grid.Lx = 2*pi;                              Grid.Ly = 2*pi;
Grid.Nx = 256;                               Grid.Ny = 256;     

Grid.dx = Grid.Lx/Grid.Nx;      Grid.dy = Grid.Ly/Grid.Ny;

x = Grid.Lx*(1:Grid.Nx)' / Grid.Nx - Grid.Lx/2;             
y = Grid.Ly*(1:Grid.Ny)' / Grid.Ny - Grid.Ly/2;

[xx,yy] = meshgrid(x,y);

kx = [ 0:nx/2-1, 0.0, -nx/2+1:-1]' / (Lx/pi/2);
ky = [ 0:ny/2-1, 0.0, -ny/2+1:-1]' / (Ly/pi/2);

[kkx,kky] = meshgrid(kx, ky);
Grid.k = sqrt(kkx.^2 + kky.^2);
Grid.k = k(:);

Grid.inv_k = 1./(k.^2);      Grid.inv_k(k == 0) = 1;


dt = 0.001;
tf = 8;

t_vals = linspace(dt, tf, round(tf / dt));




% Initial condition
u = 0.05*sin(xx).*sin(yy);


[tt, uu, Eu, Eu_SSAV, Em, mass, m_est_vals, t_vals, dt_vals] = ...
    SSAV_2D(Grid, dt, tf, Para, u, 1);


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