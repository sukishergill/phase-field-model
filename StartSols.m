function run = StartSols(IC,m, eta, eps0)

Lx = 2*pi;
Ly = 2*pi;
nx = 2^6;
ny = 2^6;
x = Lx*(1:nx)' / nx - Lx/2;
y = Ly*(1:ny)' / ny - Ly/2;
[xx,yy] = meshgrid(x,y);
dt = 0.001;
tf = 450;
dk = 1;

% Initial condition
switch IC
    case  1
        [xx,yy] = meshgrid(x,y);
        u = makerandu0(.35,xx,yy,12);
        tf = 200;
    case 2
        ny = 2;
        x = Lx*(1:nx)' / nx - Lx/2;
        y = Ly*(1:ny)' / ny - Ly/2;
        [xx,yy] = meshgrid(x,y);
        u = makerandu0(.25,xx,yy,10);
        u = mean(u,1);
        u = repmat(u,[ny,1]);
        tf = 50;
    case 3
        eps0 = .525;
        u = sin(yy).*cos(xx);
        u = u/max(abs(u(:)))*.5;
    case  4
        Lx = 2*pi;
        Ly = 2*pi;
        eps0 = .9;
        u = (rand(nx,ny)-0.5);
        tf = 100;
    case 5
        eps0 = .45;
        u = .5*cos(2*xx);
    case 6
        eps0 = .68;
        Lx = (2*pi*sqrt(2));
        Ly = (sqrt(2)*2*pi/sqrt(3));
        x = Lx*(1:nx)' / nx - Lx/2;
        y = Ly*(1:ny)' / ny - Ly/2;
        [xx,yy] = meshgrid(x,y);
        u = cos(xx*sqrt(2)) - cos(xx/sqrt(2)+sqrt(3)*yy/sqrt(2)) - ...
            cos(xx/sqrt(2)-sqrt(3)*yy/sqrt(2));
        u = u/max(abs(u(:)))*.9;
end

u = u - mean(u(:)) + m;
if max(abs(u(:))) > 1
    u = u-mean(u(:));
    u = u/(1-m) + m;
end

B = 10*Lx^5;
%B = (Lx*Ly)^2;
S = 1;

params.S = S;
params.B = B;
params.dk = dk;
params.eta = eta;
params.eps0 = eps0;
params.dt = dt;
params.m = m;
domain.nx = nx;
domain.ny = ny;
domain.Lx = Lx;
domain.Ly = Ly;

run = FCH_BDF2_SAV_JF_kappa(domain, params, u, tf);

end

%--------------------------------------------------------------------------
function u0 = makerandu0(uM,xx,yy,N)

u0 = zeros(size(xx));
[nM,nt] = size(u0);
if nM ~= nt
    nM = max(nM,nt);
    nM = floor(nM/4);
    kx1 = randi([-nM nM],1,N);
    kx2 = randi([-nM nM],1,N);
    amps1 = 4*(randn(1,N)-.5);
    amps2 = 4*(randn(1,N)-.5);
    for i=1:N
        u0 = u0 + amps1(i)*cos(xx*kx1(i)) + ...
            amps2(i)*sin(xx*kx2(i));
    end

else
    nM = floor(nM/4);
    kx1 = randi([-nM nM],1,N);
    kx2 = randi([-nM nM],1,N);
    amps1 = 4*(randn(1,N)-.5);
    amps2 = 4*(randn(1,N)-.5);
    ky1 = randi([-nM nM],1,N);
    ky2 = randi([-nM nM],1,N);

    for i=1:N
        u0 = u0 + amps1(i)*cos(xx*kx1(i)).*sin(yy*ky1(i)) + ...
            amps2(i)*cos(yy*ky2(i)).*sin(xx*kx2(i));
    end
end

u0 = uM*u0/max(abs(u0(:)));

end