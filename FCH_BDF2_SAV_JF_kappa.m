function DataOut = FCH_BDF2_SAV_JF_kappa (domain, params, u0, tf)

m = params.m;
u = u0 - mean(u0(:)) + m;
nx = domain.nx;
ny = domain.ny;
Lx = domain.Lx;
Ly = domain.Ly;
eta = params.eta;
S = params.S;
B = params.B;
dk = params.dk;
eps0 = params.eps0;
dt = params.dt;

% FCH BDF2 SAV
%
% This function solves the functionalized Cahn-Hilliard energy using the
% BDF2-SAV scheme

% spatial grid
dx = Lx/nx;                         dy = Ly/ny; dA = dx*dy;
x = Lx*(1:nx)' / nx - Lx/2;         y = Ly*(1:ny)' / ny - Ly/2;
A = Lx*Ly;
[xx, yy] = meshgrid(x, y);

kx = [ 0:nx/2-1, 0.0, -nx/2+1:-1]' / (Lx/pi/2);
ky = [ 0:ny/2-1, 0.0, -ny/2+1:-1]' / (Ly/pi/2);
if nx == 2
    kx = 0*kx;
elseif ny == 2
    ky = 0*ky;
end

[kxx,kyy] = meshgrid(kx, ky);
k2 = kxx.^2 + kyy.^2;
k4 = k2.^2;
k6 = k2.^3;

% define number of time steps and time steps for plotting purposes
nmax = round(tf / dt);
nplt = floor( (tf / 100) / dt);
tt = zeros(1, 101);
tt(1) = 0;

uu = cell(101, 1);
uu{1} = u;

mass = zeros(1, 101);
mass(1) = mean(u(:));
epsilon = mass;
epsilon(1) = eps0;

Eu = zeros(1, 101);
kappa = eps0^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f, fp] = compute_F(u);
v = fft2(u);
ux = real(ifft2(1i*kxx.*v));
uy = real(ifft2(1i*kyy.*v));

G = compute_G(u);

w_old = compute_w;

lapu = real(ifft2(-k2.*v));

H = compute_H(u, w_old);

r = w_old - 0.5 * sum(H.*u,'all')*dA;

H_hat = fft2(H);

r_hat = u/dt + (kappa*(2 + eta) + S) * ...
    real(ifft2(k4 .* v)) + ...
    r * real(ifft2(-k2 .* H_hat));


% define linear operator for BDF1
P = LinOp(dt);

r_hat = fft2(r_hat);
psi_r = P.*r_hat(:);
psi_r = real(ifft2(reshape(psi_r, ny, nx)));
psi_H = P .* (-k2(:) .* H_hat(:));
psi_H = real(ifft2(reshape(psi_H, ny, nx)));

innprod_Hu = compute_ip(H);

w = innprod_Hu + r;

u_old = u;

u = 0.5*innprod_Hu * psi_H + psi_r;
v = fft2(u);
[F, f, fp] = compute_F(u);
G = compute_G(u);
w = compute_w;
Eu(1) = compute_E(v);

if nplt == 1

    uu{2} = u;

    tt(2) = dt;

    w_vals(2) = w;

    Eu(2) = compute_E(v);

    mass(2) = mean(u(:));

end

Em = .5*(m^3-m)^2 -eta*0.25*(m^2 - 1)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux = real(ifft2(1i*kxx.*v));
uy = real(ifft2(1i*kyy.*v));
ES = Eu;
for j = 2 : nmax
    P = LinOp((2*dt)/3);

    t = j * dt;

    u_snew = compute_u_snew(u, u_old);

    [F, f, fp] = compute_F(u_snew);
    G = compute_G(u_snew);
    w_snew = compute_w;

    H_snew = compute_H(u_snew, w_snew);
    H_hat = fft2(H_snew);
    u_snew = fft2(u_snew);
    w_star = (4*w - w_old)/3;
    u_star = (4*u - u_old)/3;

    [r, psi_r, psi_H] = psiFUNC;

    innprod_Hu = compute_ip(H_snew);

    w_new = 0.5*innprod_Hu + r;

    u_new = 0.5*innprod_Hu * psi_H + psi_r;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;

    v = fft2(u);
    ux = real(ifft2(1i*kxx.*v));
    uy = real(ifft2(1i*kyy.*v));
    lapu = real(ifft2(-k2.*v));

    if j > 500 && dk > 0
        dkappa = mean(eta/2*(ux.^2+uy.^2)+lapu.*f-kappa*lapu.^2,'all');
        if kappa + dk*dt*dkappa > 0
            kappa = max(.9*kappa,min(kappa + dk*dt*dkappa,1.1*kappa));
        else
            kappa = max(.1,.75*kappa);
        end
    end

    if (mod(j, nplt) == 0)

        uu{j/nplt + 1} = u;

        [F, f, fp] = compute_F(u);
        G = compute_G(u);
        w = compute_w;

        Eu(j/nplt + 1) = compute_E(v);

        ES(j/nplt + 1) = compute_E2;

        mass(j/nplt + 1) = mean(u(:));

        tt(j/nplt + 1) = t;

        w_vals(j/nplt + 1) = w;

        epsilon(j/nplt + 1) = sqrt(kappa);
    end

end

run.u = uu;
run.t = tt;
run.Eu = Eu;
run.Em = Em;
run.Es = ES;
run.epsilon = epsilon;
domain.Lx = Lx;
domain.Ly = Ly;
domain.x = xx;
domain.y = yy;

DataOut.Domain = domain;
DataOut.Run = run;
DataOut.Params = params;

%%
    function P = LinOp(DT)

        P = 1/DT + kappa^2*k6(:) + S*k4(:);
        P = 1./P;

    end


    function u_snew = compute_u_snew(u, u_old)

        u_snew = 2*u - u_old;

    end

    function [F, f, fp] = compute_F(u)

        F = 0.25 * (u.^2 - 1).^2;

        f = u.^3 - u;

        fp = 3 * u.^2 - 1;

    end

    function G_loc = compute_G(U)

        G_loc = sum(3*kappa * U.^2 .* (ux.^2 + uy.^2) + ...
            0.5 * f.^2 - eta * F,'all') *  dA;

    end

    function W = compute_w

        W = sqrt(G + B);

    end

    function H_loc = compute_H(U, W)

        v2 = fft2(U.^2);
        ux2 = real(ifft2(1i*kxx.*v2));
        uy2 = real(ifft2(1i*kyy.*v2));
        u_grad1 = ux.^2 + uy.^2;

        H_loc = (6*kappa * (u.* u_grad1 -  ux.*ux2 - uy.*uy2 - ...
            U.^2 .* real(ifft2(-k2 .* v))) + f .* fp - eta * f) / W;

    end


    function innprod_Hu = compute_ip(HH)

        innprod_Hu = sum(HH .* psi_r,'all')*dA / ...
            (1 - 0.5*sum(HH .* psi_H,'all') * dA);

    end

    function   [R, Psi_r, Psi_H] = psiFUNC

        R = w_star - 0.5 * sum(H_snew .* u_star,'all') ...
            * dA;

        r_hat = u_star*3/(2*dt) + (kappa*(2 + eta) + S) * ...
            real(ifft2(k4 .* u_snew)) + ...
            r * real(ifft2(-k2 .* H_hat));

        r_hat = fft2(r_hat);

        Psi_r = P .* r_hat(:);
        Psi_r = reshape(Psi_r, ny, nx);
        Psi_r = real(ifft2(Psi_r));

        Psi_H = P.*(-k2(:).*H_hat(:));
        Psi_H = reshape(Psi_H, ny, nx);
        Psi_H = real(ifft2(Psi_H));

    end

    function E = compute_E(V)

        du = ux.^2 + uy.^2;
        lapu = real(ifft2(-k2.*V));
        Et = .5*(-kappa*lapu+f).^2 - eta*(kappa*du/2 + F);
        E = sum(Et*dA,'all');

        %E = sum(0.5*kappa^2 * (real(ifft2(-k2 .* V))).^2 - ...
        %    (eta/2 + 1)*kappa * du,'all') + w^2 - B;
        E = E/A;
    end

    function E2 = compute_E2

        V = fft2(u);
        Vo = fft2(u_old);

        du = ux.^2 + uy.^2;
        lapu = real(ifft2(-k2.*V));

        ux_old = real(ifft2(1i*kxx.*Vo));
        uy_old = real(ifft2(1i*kyy.*Vo));
        lapu_old = real(ifft2(-k2.*Vo));

        E2 = ...
            kappa^2/2*(sum(lapu.^2,'all')*dA + sum((2*lapu-lapu_old).^2,'all'))/2*dA + ...
            (S+kappa*(2+eta))/2*sum((ux-ux_old).^2+(uy-uy_old).^2,'all')*dA - ...
            kappa*(2+eta)/2*(sum(du,'all') + sum((2*ux-ux_old).^2+(2*uy-uy_old).^2,'all'))*dA/2 + ...
            (w^2+(2*w-w_old)^2)/2 - B;
        E2 = E2/A;

    end

end