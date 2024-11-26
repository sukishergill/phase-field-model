function [ xx, yy, kxx, kyy, tt, uu, Eu, Em, mass] = FCH_BDF2_SAV_JF ( nx, ny,...
    Lx, Ly, dt, tf, epsilon, eta, B, S, u)

% FCH BDF2 SAV
%
% This function solves the functionalized Cahn-Hilliard energy using the
% BDF2-SAV scheme

% spatial grid
dx = Lx/nx;                         dy = Ly/ny; dA = dx*dy;
x = Lx*(1:nx)' / nx - Lx/2;         y = Ly*(1:ny)' / ny - Ly/2;
[xx, yy] = meshgrid(x, y);

kx = [ 0:nx/2-1, 0.0, -nx/2+1:-1]' / (Lx/pi/2);
ky = [ 0:ny/2-1, 0.0, -ny/2+1:-1]' / (Ly/pi/2);

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

Eu = zeros(1, 101);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f, fp] = compute_F(u);
v = fft2(u);
ux = real(ifft2(1i*kxx.*v));
uy = real(ifft2(1i*kyy.*v));

G = compute_G(u);

w_old = compute_w;

Eu(1) = compute_E(v);
lapu = real(ifft2(-k2.*v));

H = compute_H(u, w_old);

r = w_old - 0.5 * sum(H.*u,'all')*dA;

H_hat = fft2(H);

r_hat = u/dt + (epsilon^2*(2 + eta) + S) * ...
    real(ifft2(k4 .* v)) + ...
    r * real(ifft2(-k2 .* H));


% define linear operator for BDF1
P = LinOp(dt);

r_hat = fft2(r_hat);
psi_r = P.*r_hat(:);
psi_r = real(ifft2(reshape(psi_r, nx, ny)));
psi_H = P .* (-k2(:) .* H_hat(:));
psi_H = real(ifft2(reshape(psi_H, nx, ny)));

innprod_Hu = compute_ip(H);

w = innprod_Hu + r;

u_old = u;

u = 0.5*innprod_Hu * psi_H + psi_r;
v = fft2(u);

if nplt == 1

    uu{2} = u;

    tt(2) = dt;

    w_vals(2) = w;

    Eu(2) = compute_E(v);

    mass(2) = mean(u(:));

end

Em = .5*(mass(1)^3-mass(1))^2 -eta*0.25 * (mass(1)^2 - 1)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define P for BDF2
P = LinOp((2*dt)/3);
ux = real(ifft2(1i*kxx.*v));
uy = real(ifft2(1i*kyy.*v));

for j = 2 : nmax

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

    if (mod(j, nplt) == 0)

        Eu(j/nplt + 1) = compute_E(v);

        mass(j/nplt + 1) = mean(u(:));

        uu{j/nplt + 1} = u;

        tt(j/nplt + 1) = t;

        w_vals(j/nplt + 1) = w;
    end

end

% Eu = Eu / (Lx*Ly);

%%
    function P = LinOp(DT)

        P = 1/DT + epsilon^4*k6(:) + S*k4(:);
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

    function G = compute_G(U)

        G = sum(3*epsilon^2 * U.^2 .* (ux.^2 + uy.^2) + ...
            0.5 * f.^2 - eta * F,'all') *  dA;

    end

    function w = compute_w

        w = sqrt(G + B);

    end

    function H = compute_H(U, W)

        v2 = fft2(U.^2);
        ux2 = real(ifft2(1i*kxx.*v2));
        uy2 = real(ifft2(1i*kyy.*v2));
        u_grad1 = ux.^2 + uy.^2;

        H = (6*epsilon^2 * (u.* u_grad1 -  ux.*ux2 - uy.*uy2 - ...
            U.^2 .* real(ifft2(-k2 .* v))) + f .* fp - eta * f) / W;

    end


    function innprod_Hu = compute_ip(HH)

        innprod_Hu = sum(HH .* psi_r,'all')*dA / ...
            (1 - 0.5*sum(HH .* psi_H,'all') * dA);

    end

    function   [R, Psi_r, Psi_H] = psiFUNC

        R = w_star - 0.5 * sum(H_snew .* u_star,'all') ...
            * dA;

        r_hat = u_star*3/(2*dt) + (epsilon^2*(2 + eta) + S) * ...
            real(ifft2(k4 .* u_snew)) + ...
            r * real(ifft2(-k2 .* H_hat));

        r_hat = fft2(r_hat);

        Psi_r = P .* r_hat(:);
        Psi_r = reshape(Psi_r, nx, ny);
        Psi_r = real(ifft2(Psi_r));

        Psi_H = P.*(-k2(:).*H_hat(:));
        Psi_H = reshape(Psi_H, nx, ny);
        Psi_H = real(ifft2(Psi_H));

    end

    function E = compute_E(V)

        du = ux.^2 + uy.^2;
        lapu = real(ifft2(-k2.*V));
        Et = .5*(-epsilon^2*lapu+f).^2 - eta*(epsilon^2*du/2 + F);
        E = mean(Et,'all');

        % E = sum(0.5*epsilon^4 * (real(ifft2(-k2 .* V))).^2 - ...
        %     (eta/2 + 1)*epsilon^2 * du,'all')*...
        %     dA + W^2 - B;

    end

end