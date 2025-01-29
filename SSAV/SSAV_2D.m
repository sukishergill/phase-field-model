function [u, Eu, Eu_SSAV, Em, mass, t_vals]...
    = SSAV_2D (Grid, Time, Para, u, model)

% This function solves the following PDE
%
%           u_t = G(-epsilon*Du + f(u) + alpha*epsilon v)
%
% where f(u) = u^3 - beta * u
%
% We implemented the stabilized scalar auxiliary variable (SSAV) method.
%
% We use BDF2 or adaptive time stepping scheme for the time
% discretization and Fourier method for the spatial discretization.
%
% References:
%
%           "The scalar auxiliary variable (SAV) approach for gradient
%           flows"
%           Jie Shen, Jie Xu and Jiang Yang
%
%           "Efficient and energy stable method for the Cahn-Hilliard
%           phase-field model for diblock copolymers"
%           Jun Zhang, Chuanjun Chen and Xiaofeng Yang
%
%           "Benchmark computation of morphological complexity in the
%           functionalized Cahn-Hilliard gradient flow"
%           Andrew Christlieb, Keith Promislow, Zengqiang Tan, Sulin Wang,
%           Brian Wetton, and Steven M. Wise
% 
%
% The input variable will indicate which eqn to solve
%           1. Ohta-Kawasaki or Cahn-Hilliard
%           2. Phase-field-crystals
%           3. Allen-Cahn

% err_tol = 10E-5;

dt = Time.dt_min;

k2 = Grid.kxx.^2 + Grid.kyy.^2;
k4 = k2.^2;         k6 = k2.^3.;


if model == 2

    D = (-k2 + 1);    % PFC
    
else
    
    D = -1i*Grid.k;          % AC, CH and OK
end

% For all 3 models the operator G is the same
if model == 3
    
    G = -1;             % AC
    
else
    G = -k2;          % PFC, CH and OK
    
end    

% define number of time steps and time steps for plotting purposes
nplt = floor( (Time.tf / 100) / Time.dt_min);
nmax = round(Time.tf / Time.dt_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f] = compute_F(u);         int_F = sum(F(:)) * Grid.dx*Grid.dy;
% f = compute_f(u, Para.beta);

% define the following as w_old b/c we will need this in the for loop
w_old = compute_w(int_F);   

t = 0;

mass = zeros(1, 101);
mass(1) = (sum(u(:))*Grid.dx*Grid.dy)/(Grid.Lx*Grid.Ly);

Em = Grid.Lx*Grid.Ly * compute_F(Para.m);

Eu = zeros(1, nmax); 
Eu_SSAV = zeros(1, nmax);
% Eu_SSAV_test = zeros(1, nmax);
t_vals = zeros(1, nmax);
% uu = cell(101, 1);
% uu{1} = u;

Eu(1) = compute_En(u, w_old);
% Eu_SSAV_test(1) = (Eu(1))/2 + (compute_En(2*u - u, 2*w_old - w_old)) / 2 + ...
%     sum(Para.S/2 * (u - u).^2, 'all')*Grid.dx*Grid.dy;

H = compute_H(f, w_old);

r = -0.5 * sum(H.*u, 'all') * Grid.dx*Grid.dy + w_old;
H2 = fft2(H);

ufft = fft2(u);

r_tilde = u/dt - Para.S*real(ifft2(G .* ufft)) - r*real(ifft2(G .* H2));

r_hat = (Para.alpha*Para.epsilon^2*dt)/(Grid.Lx*Grid.Ly) * ...
    sum(r_tilde(:))*Grid.dx*Grid.dy + r_tilde;

if model == 2

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G(:) - Para.epsilon^2*G(:).*D(:).^2;

else

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G(:) + Para.epsilon^2*G(:).*D(:).^2;

end

P = spdiags(1/dt + P1, 0, Grid.Nx*Grid.Ny, Grid.Nx*Grid.Ny);

r_hat = fft2(r_hat);
psi_r = P \ r_hat(:);
psi_r = real(ifft2(reshape(psi_r, Grid.Nx, Grid.Ny)));
  
psi_H = P \ (G(:) .* H2(:)); 
psi_H = real(ifft2(reshape(psi_H, Grid.Nx, Grid.Ny)));

innprod_Hu = compute_ip(H, psi_r, psi_H);

w = 0.5*innprod_Hu + r;

u_old = u;

u = 0.5*innprod_Hu * psi_H + psi_r;

t = t + dt;

Eu(2) = compute_En(u, w);
% Eu_SSAV_test(2) = (Eu(2))/2 + (compute_En(2*u - u, 2*w - w)) / 2 + ...
%     sum(sum(Para.S/2 * (u(:) - u(:)).^2))*Grid.dx*Grid.dy;
Eu_SSAV(1) = (Eu(2))/2 + (compute_En(2*u - u_old, 2*w - w_old)) / 2 + ...
    sum(sum(Para.S/2 * (u(:) - u_old(:)).^2))*Grid.dx*Grid.dy;

E_old = Eu(1);       E = Eu(2);       t_vals(2) = t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;


while t < Time.tf

    j = j + 1;

    dt_new = dt;
    
    t = t + dt_new;

    [u_new, w_new, gamma] = compute_unew(u, u_old, w, w_old, ...
        dt, dt_new);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_new = compute_En(u_new, w_new);

    E_mod = E_new/2 + (compute_En(2*u_new - u, 2*w_new - w)) / 2 + ...
        sum(sum(Para.S/2 * (u(:) - u_old(:)).^2))*Grid.dx*Grid.dy;

    rel_err = abs(E_new - E_mod);

    adap = 0.5*dt*(Para.err_tol / rel_err);

    dt_new = max(Time.dt_min, min(adap, Time.dt_max));

    l = 1;

    while (rel_err > Para.err_tol) && (dt_new ~= Time.dt_min)

        [u_new, w_new, gamma] = compute_unew(u, u_old, w, w_old, ...
            dt, dt_new);

        E_new = compute_En(u_new, w_new);
        
   
        E_mod = E_new/2 + (compute_En(2*u_new - u, 2*w_new - w)) / 2 + ...
            sum(sum(Para.S/2 * (u(:) - u_old(:)).^2))*Grid.dx*Grid.dy;
    
        rel_err = abs(E_new - E_mod);

        adap = 0.5*dt*(Para.err_tol / rel_err);

        dt_new = max(Time.dt_min, min(adap, Time.dt_max));

        l = l + 1;

        if l == 10
            dt_new = Time.dt_min;
        end
    end

    % compute new time step
    dt = dt_new;
 
    Eu(j + 1) = E_new;

    % Eu_SSAV_test(j+1) = (Eu(j+1))/2 + (compute_En(2*u_new - u_new,...
    %     2*w_new - w_new)) / 2 + ...
    %     sum(sum(Para.S/2 * (u(:) - u(:)).^2))*Grid.dx*Grid.dy;
    Eu_SSAV(j) = (Eu(j+1))/2 + (compute_En(2*u_new - u, 2*w_new - w)) / 2 + ...
        sum(sum(Para.S/2 * (u(:) - u_old(:)).^2))*Grid.dx*Grid.dy;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;

    % uu{j + 1} = u;
    % if (mod(j, nplt) == 0)                  
    % 
    %     uu{j/nplt + 1} = u;
    % 
    % end

    t_vals(j + 1) = t;

end

nt = sum(t_vals > 0) + 1;

% Em = Em; % / (Grid.Lx*Grid.Ly);
% Eu = Eu(1:nt); %./ (Grid.Lx*Grid.Ly);
% t_vals = t_vals(1:nt);
% Eu_SSAV = Eu_SSAV(1:nt);

function u_snew = compute_u_snew(u, u_old)

u_snew = 2*u - u_old;

end

function [F, f] = compute_F(u)

f = u.^3 - Para.beta*u;

F = 0.25 * (u.^2 - Para.beta).^2;

end


function w = compute_w(int_F)

w = sqrt(int_F + Para.B);

end

function H = compute_H(f, w)

H = f / w;

end

function [r, r_hat] = compute_r(u_snew, u, u_old, w, w_old, ... 
    H_snew, dt, a, b, c)

r = -(b*w + c*w_old)/a + 0.5 * sum(sum((H_snew .* (b*u + c*u_old)/a)))...
    * Grid.dx*Grid.dy;

rH = fft2(r .* H_snew);
u_snew = fft2(u_snew);

r_tilde = -(b*u + c*u_old) - Para.S*real(ifft2(G .* u_snew)) + ...
    real(ifft2(G .* rH));

m_est = (Para.alpha*Para.epsilon^2)/(a*Grid.Lx*Grid.Ly) * ...
    sum(r_tilde(:)) * Grid.dx*Grid.dy;

r_hat = r_tilde + m_est;

end


function innprod_Hu = compute_ip(H_snew, psi_r, psi_H)
% This fuction computes the inner product of H^{*,n+1} and u^{n+1}

innprod_Hu = sum(sum(H_snew .* psi_r)) * Grid.dx*Grid.dy...
    / (1 - 0.5*sum(sum(H_snew .* psi_H)) * Grid.dx*Grid.dy);

end

function E_n = compute_En(u, w)

e = Para.epsilon^2/2;

um = fft2(u - Para.m);

v = -Grid.inv_k .* um;
v(Grid.k == 0) = 0;

vx = real(ifft2(-1i*Grid.kxx.*v));         
vy = real(ifft2(-1i*Grid.kyy.*v));

dv = vx.^2 + vy.^2;

u = fft2(u);

if model == 2

    du = D.*u;
else

    ux = real(ifft2(-1i*Grid.kxx.*u));         
    uy = real(ifft2(-1i*Grid.kyy.*u));

    du = ux.^2 + uy.^2;
end

E_n = sum(sum(e *du + Para.alpha*e * dv))*Grid.dx*Grid.dy +  w^2 - Para.B;

end

function [a, b, c, gamma] = compute_dt_coeffs(dt, dt_new)

gamma = dt_new / dt;

% a = (1 + 2*gamma) / (1 + gamma);            a = a / dt_new;

a = 1 / dt_new + 1 / (dt_new + dt);

% b = -(1 + gamma)^2 / (1 + gamma);           b = b / dt_new;
b = - 1 / dt_new - 1 / dt;

% c = gamma^2 / (1 + gamma);                  c = c / dt_new;
c = 1 / dt - 1 / (dt_new + dt);

end


function A_dp = compute_A_dp(rel_err, dt)

A_dp = rho_s * dt * (Para.err_tol / rel_err) ^ (1/3);

end

function [u_new, w_new, gamma] = compute_unew(u, u_old, w, w_old,...
    dt, dt_new)

u_snew = compute_u_snew(u, u_old);
    
[F, f] = compute_F(u_snew);
% f = compute_f(u_snew, Para.beta);

int_F = sum(F(:)) * Grid.dx*Grid.dy;

%w(u^{*,n+1})
w_snew = compute_w(int_F);

% H^{*,n+1}
H_snew = compute_H(f, w_snew);

[a, b, c, gamma] =  compute_dt_coeffs(dt, dt_new);

[r, r_hat] = compute_r(u_snew, u, u_old, w, w_old, H_snew, ...
    dt_new, a, b, c);

% define P for BDF2    
P = spdiags(a + P1, 0, Grid.Nx*Grid.Ny, Grid.Nx*Grid.Ny);

r_hat = fft2(r_hat);
psi_r = P \ r_hat(:);           psi_r = reshape(psi_r, Grid.Nx, Grid.Ny);
psi_r = real(ifft2(psi_r));

H2 = fft2(H_snew);
psi_H = P \ (G(:).*H2(:));     psi_H = reshape(psi_H, Grid.Nx, Grid.Ny);
psi_H = real(ifft2(psi_H));

innprod_Hu = compute_ip(H_snew, psi_r, psi_H);

w_new = 0.5*innprod_Hu + r;      

u_new = 0.5*innprod_Hu * psi_H + psi_r;

end

end