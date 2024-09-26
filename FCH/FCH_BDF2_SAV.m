function [ xx, yy, k, tt, uu, Eu, Em, mass] = FCH_BDF2_SAV ( nx, ny,...
    Lx, Ly, dt, tf, epsilon, eta, B, S, u)

% FCH BDF2 SAV
%
% This function solves the functionalized Cahn-Hilliard energy using the
% BDF2-SAV scheme

% spatial grid
dx = Lx/nx;                         dy = Ly/ny;
x = Lx*(1:nx)' / nx - Lx/2;         y = Ly*(1:ny)' / ny - Ly/2;      

[xx, yy] = meshgrid(x, y);

kx = [ 0:nx/2-1, 0.0, -nx/2+1:-1]' / (Lx/pi/2);
ky = [ 0:ny/2-1, 0.0, -ny/2+1:-1]' / (Ly/pi/2);

[kkx,kky] = meshgrid(kx, ky);
k = sqrt(kkx.^2 + kky.^2);
k = k(:);

% define number of time steps and time steps for plotting purposes
nmax = round(tf / dt);
nplt = floor( (tf / 100) / dt);
tt = zeros(1, 101);
tt(1) = 0;

uu = cell(101, 1);
uu{1} = u;

mass = zeros(1, 101);
mass(1) = (sum(u(:))*dx*dy)/(Lx*Ly);

Eu = zeros(1, 101);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f, fp] = compute_F(u);

G = compute_G(u, epsilon, eta, F, f, k, dx, dy, nx, ny);

w_old = compute_w(G, B);

Eu(1) = compute_E(u, w_old, epsilon, eta, B, dy, dy, k);

H = compute_H(u, w_old, epsilon, eta, f, fp, x, y, xx, yy, nx, ny, k);

% r = w_old - 0.5 * sum(sum(H.*u))*dx*dy;

H2 = fftn(H);

% rhat = u/dt + (epsilon^2*(2 + eta) + S) * ...
%     real(ifftn(reshape(k.^4 .* u(:), nx, ny))) + ...
%     r * real(ifftn( reshape(-k.^2 .* H(:), nx, ny)));

r = u/dt + (epsilon^2*(2 + eta) + S) * ...
    real(ifftn(reshape(k.^4 .* u(:), nx, ny))) + ...
    (w_old - 0.5 * sum(sum(H.*u))*dx*dy) * ...
    real(ifftn( reshape(-k.^2 .* H(:), nx, ny)));

P = spdiags(1/dt + epsilon^4*k.^6 + S.*k.^4, 0, nx*ny, nx*ny);

r = fftn(r);
psi_r = P \ r(:);
psi_r = real(ifftn(reshape(psi_r, nx, ny)));

psi_H = P \ (-k.^2 .* H2(:));         
psi_H = real(ifftn(reshape(psi_H, nx, ny)));

innprod_Hu = compute_ip(H, psi_r, psi_H, dx, dy);

w = innprod_Hu + w_old - 0.5 * sum(sum(H.*u))*dx*dy;

u_old = u;

u = 0.5*innprod_Hu * psi_H + psi_r;

% define P for BDF2
P = spdiags(3/(2*dt) + epsilon*k.^6 + S*k.^4, 0, nx*ny, nx*ny);

if nplt == 1
    
    uu{2} = u;    
    
    tt(2) = dt;  

    w_vals(2) = w;
    
    Eu(2) = compute_E(u, w, epsilon, eta, B, dx, dy, k);
    
    mass(2) = sum(sum(u))*dx*dy/(Lx*Ly);
    
end

Em = sum(sum(0.25 * (mass(1)*ones(size(u)).^2 - 1).^2))*dx*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2 : nmax

    t = j * dt;

    u_snew = compute_u_snew(u, u_old);

    [F, f, fp] = compute_F(u_snew);

    G = compute_G(u_snew, epsilon, eta, F, f, k, dx, dy, nx, ny);

    w_snew = compute_w(G, B);

    H_snew = compute_H(u_snew, w_snew, epsilon, eta, f, fp, x, y,...
        xx, yy, nx, ny, k);

    H = fftn(H_snew);
    u_snew = fftn(u_snew);

    r = (4*w - w_old)/3 - 0.5 * sum(sum(H_snew .* (4*u - u_old)/3)) ...
        * dx*dy;

    r2 = (4*u - u_old)/(2*dt) + (epsilon^2*(2 + eta) + S) * ...
        real(ifftn(reshape(-k.^2 .* u_snew(:), nx, ny))) + ...
        r * real(ifftn( reshape(-k.^2 .* H(:), nx, ny)));

    r2 = fftn(r2);
    H2 = fftn(H_snew);
    
    psi_r = P \ r2(:);           psi_r = reshape(psi_r, nx, ny);
    psi_r = real(ifftn(psi_r));
    psi_H = P \ (-k.^2.*H2(:));     psi_H = reshape(psi_H, nx, ny);
    psi_H = real(ifftn(psi_H));

    innprod_Hu = compute_ip(H_snew, psi_r, psi_H, dx, dy);

    w_new = 0.5*innprod_Hu + r; 

    u_new = 0.5*innprod_Hu * psi_H + psi_r;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;

    if (mod(j, nplt) == 0)                  

        Eu(j/nplt + 1) = compute_E(u, w, epsilon, eta, B, dx, dy, k);
        
        mass(j/nplt + 1) = sum(sum(u))*dx*dy/(Lx*Ly);
        
        uu{j/nplt + 1} = u;
        
        tt(j/nplt + 1) = t;

        w_vals(j/nplt + 1) = w;
    end

end

Eu = Eu / (Lx*Ly);
Em = Em / (Lx*Ly);

end

function u_snew = compute_u_snew(u, u_old)

u_snew = 2*u - u_old;

end

function [F, f, fp] = compute_F(u)

F = 0.25 * (u.^2 - 1).^2;

f = u.^3 - u;

fp = 3 * u.^2 - 1;

end

function G = compute_G(u, epsilon, eta, F, f, k, dx, dy, nx, ny)

k = reshape(k, nx, ny);

v = fftn(u);

G = sum(sum(3*epsilon^2 * u.^2 .* real(ifftn(-1i*k.*v)).^2 + ...
    0.5 * f.^2 - eta * F)) *  dx*dy;

end

function w = compute_w(G, B)

w = sqrt(G + B);

end

function H = compute_H(u, w, epsilon, eta, f, fp, x, y, xx, yy, nx, ny, k)

k = reshape(k, nx, ny);
v = fftn(u);

[ux, uy] = gradient(u, x, y);

div_u = divergence(xx, yy, u.^2.*ux, u.^2.*uy);

H = (6*epsilon^2 * (u.* real(ifftn(-1i*k.*v)).^2 - div_u) + f .* fp - ...
    eta * f) / w;

end


function innprod_Hu = compute_ip(H_snew, psi_r, psi_H, dx, dy)

innprod_Hu = sum(sum(H_snew .* psi_r))*dx*dy / ...
    (1 - 0.5*sum(sum(H_snew .* psi_H)) * dx*dy);

end

function E = compute_E(u, w, epsilon, eta, B, dx, dy, k)

u = fftn(u);

E = sum(sum(0.5*epsilon^2 * (real(ifftn(-k.^2 .* u(:)))).^2 - ...
    (eta/2 + 1)*epsilon^2 * real(ifftn(-1i*k.*u(:))).^2))*...
    dx*dy + w^2 - B;

end