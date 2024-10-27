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

[kxx,kyy] = meshgrid(kx, ky);
k = sqrt(kxx.^2 + kyy.^2);
% k = k(:);

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

G = compute_G(u, epsilon, eta, F, f, kxx, kyy, dx, dy, nx, ny, x, y);

w_old = compute_w(G, B);

Eu(1) = compute_E(u, w_old, epsilon, eta, B, dy, dy, k, kxx, kyy, nx, ny);

H = compute_H(u, w_old, epsilon, eta, f, fp, x, y, xx, yy, nx, ny, kxx, kyy, k);
H_old = H;

r = w_old - 0.5 * sum(sum(H.*u))*dx*dy;


H2 = fft2(H);

v = fft2(u);

r_hat = u/dt + (epsilon^2*(2 + eta) + S) * ...
    real(ifft2(k.^4 .* v)) + ...
    r * real(ifft2(-k.^2 .* H));


% define linear operator for BDF1
P = spdiags(1/dt + epsilon^4*k(:).^6 + S*k(:).^4, 0, nx*ny, nx*ny);

r_hat = fft2(r_hat);
psi_r = P \ r_hat(:);
psi_r = real(ifft2(reshape(psi_r, nx, ny)));

psi_H = P \ (-k(:).^2 .* H2(:));         
psi_H = real(ifft2(reshape(psi_H, nx, ny)));

innprod_Hu = compute_ip(H, psi_r, psi_H, dx, dy, x, y);

w = innprod_Hu + r;

u_old = u;

u = 0.5*innprod_Hu * psi_H + psi_r;

% define P for BDF2
P = spdiags(3/(2*dt) + epsilon^4*k(:).^6 + S*k(:).^4, 0, nx*ny, nx*ny);

if nplt == 1
    
    uu{2} = u;    
    
    tt(2) = dt;  

    w_vals(2) = w;
    
    Eu(2) = compute_E(u, w, epsilon, eta, B, dx, dy, k, kxx, kyy, nx, ny);
    
    mass(2) = mean(u(:));
    
end

Em = 0.25 * (mass(1)^2 - 1)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2 : nmax

    t = j * dt;

    u_snew = compute_u_snew(u, u_old);

    [F, f, fp] = compute_F(u_snew);

    G = compute_G(u_snew, epsilon, eta, F, f, kxx, kyy, dx, dy, nx, ny, x, y);

    w_snew = compute_w(G, B);

    H_snew = compute_H(u_snew, w_snew, epsilon, eta, f, fp, x, y,...
        xx, yy, nx, ny, kxx, kyy, k);

    % H_1 = compute_H(u, w, epsilon, eta, f, fp, x, y, xx, yy, nx, ny, k);

    % H_snew = 2*H_1 - H_old;

    % H_old = H_1;

    H = fft2(H_snew);
    u_snew = fft2(u_snew);

    r = (4*w - w_old)/3 - 0.5 * sum(sum(H_snew .* (4*u - u_old)/3)) ...
        * dx*dy;


    r_hat = (4*u - u_old)/(2*dt) + (epsilon^2*(2 + eta) + S) * ...
        real(ifft2(k.^4 .* u_snew)) + ...
        r * real(ifft2(-k.^2 .* H));

    r_hat = fft2(r_hat);    
    
    psi_r = P \ r_hat(:);           psi_r = reshape(psi_r, nx, ny);
    psi_r = real(ifft2(psi_r));
    psi_H = P \ (-k(:).^2.*H(:));     psi_H = reshape(psi_H, nx, ny);
    psi_H = real(ifft2(psi_H));

    innprod_Hu = compute_ip(H_snew, psi_r, psi_H, dx, dy, x, y);

    w_new = 0.5*innprod_Hu + r; 

    u_new = 0.5*innprod_Hu * psi_H + psi_r;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;

    if (mod(j, nplt) == 0)                  

        Eu(j/nplt + 1) = compute_E(u, w, epsilon, eta, B, dx, dy, k, kxx, kyy, nx, ny);
        
        mass(j/nplt + 1) = mean(u(:));
        
        uu{j/nplt + 1} = u;
        
        tt(j/nplt + 1) = t;

        w_vals(j/nplt + 1) = w;
    end

end

Eu = Eu / (Lx*Ly);

end

function u_snew = compute_u_snew(u, u_old)

u_snew = 2*u - u_old;

end

function [F, f, fp] = compute_F(u)

F = 0.25 * (u.^2 - 1).^2;

f = u.^3 - u;

fp = 3 * u.^2 - 1;

end

function G = compute_G(u, epsilon, eta, F, f, kxx, kyy, dx, dy, nx, ny, x, y)

v = fft2(u);

ux = ifft2(-1i*kxx.*v);         
uy = ifft2(-1i*kyy.*v);

% [ux, uy] = gradient(u, x, y);

G = sum(sum(3*epsilon^2 * u.^2 .* (ux.^2 + uy.^2) + ...
    0.5 * f.^2 - eta * F)) *  dx*dy;
% 
% G = sum(sum(3*epsilon^2 * u.^2 .* real(ifftn(-1i*k.*v)).^2 + ...
%     0.5 * f.^2 - eta * F)) *  dx*dy;


end

function w = compute_w(G, B)

w = sqrt(G + B);

end

function H = compute_H(u, w, epsilon, eta, f, fp, x, y, xx, yy, nx, ny, kxx, kyy, k)

v = fft2(u);
v2 = fft2(u.^2);

ux = real(ifft2(-1i*kxx.*v));         
uy = real(ifft2(-1i*kyy.*v));

ux2 = real(ifft2(-1i*kxx.*v2));         
uy2 = real(ifft2(-1i*kyy.*v2));

% [ux, uy] = gradient(u, x, y);
% [ux2, uy2] = gradient(u.^2, x, y);

% div_u = divergence(xx, yy, u.^2.*ux, u.^2.*uy);
% u_div_u = u.^2 .* (divergence(xx, yy, ux, uy));

% H = (6*epsilon^2 * (u.* real(ifftn(-1i*k.*v)).^2 - div_u) + f .* fp - ...
%     eta * f) / w;

u_grad1 = ux.^2 + uy.^2;

H = (6*epsilon^2 * (u.* u_grad1 -  ux.*ux2 - uy.*uy2 - ...
    u.^2 .* real(ifft2(-k.^2 .* v))) + f .* fp - eta * f) / w;

% H = (6*epsilon^2 * (u.* real(ifftn(-1i*k.*v)).^2 - (u_div_u + ux.*ux2 + ...
%     uy .* uy2)) + f .* fp - eta * f) / w;

% H = (6*epsilon^2 * (u.^2 .* real(ifftn(-k.^2.*v)) -...
%     2*u.* real(ifftn(-1i*k.*v)).^2) + f .* fp - eta * f) / w;

end


function innprod_Hu = compute_ip(H_snew, psi_r, psi_H, dx, dy, x, y)

innprod_Hu = sum(sum(H_snew .* psi_r))*dx*dy / ...
    (1 - 0.5*sum(sum(H_snew .* psi_H)) * dx*dy);


end

function E = compute_E(u, w, epsilon, eta, B, dx, dy, k, kxx, kyy, nx, ny)

u = fft2(u);
ux = real(ifft2(-1i*kxx.*u));         
uy = real(ifft2(-1i*kyy.*u));

du = ux.^2 + uy.^2;

E = sum(sum(0.5*epsilon^4 * (real(ifft2(-k.^2 .* u))).^2 - ...
    (eta/2 + 1)*epsilon^2 * du))*...
    dx*dy + w^2 - B;

end