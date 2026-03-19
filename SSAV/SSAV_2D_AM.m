function Results = SSAV_2D_AM (Grid, Time, Para, u, model)

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

num_fft = 0;

rho_s = 0.9;

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

nmax = round(Time.tf / Time.dt_min);

plt = linspace(0, Time.tf, 101);
plt_idx = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f] = compute_F(u);         int_F = sum(F(:)) * prod(Grid.d);
% f = compute_f(u, Para.beta);

u_fft = fft2(u);        num_fft = num_fft + 1;

F_tilde_old = compute_F_tilde(u, u_fft, f);

% define the following as w_old b/c we will need this in the for loop
w_old = compute_w(int_F);   

t = 0;

mass = zeros(1, 101);
mass(1) = (sum(u(:))*prod(Grid.d))/(prod(Grid.L));

Em = prod(Grid.L) * compute_F(Para.m);

Eu = zeros(1, nmax);
Eu_SSAV = zeros(1, nmax);
t_vals = zeros(1, nmax);
l_vals = zeros(1, nmax);
uu = cell(101, 1);
uu{1} = u;
u_times = zeros(1, 101);

rel_err_vals = zeros(1,nmax);

Eu(1) = compute_En(u_fft, w_old);

H = compute_H(f, w_old);

r = -0.5 * sum(H.*u, 'all') * prod(Grid.d) + w_old;
H2 = fft2(H);           num_fft = num_fft + 1;

r_tilde = u/dt - Para.S*real(ifft2(G .* u_fft)) - r*real(ifft2(G .* H2));

r_hat = (Para.alpha*Para.epsilon^2*dt)/(prod(Grid.L)) * ...
    sum(r_tilde(:))*prod(Grid.d) + r_tilde;

if model == 2

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G - Para.epsilon^2*G.*D.^2;

else

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G + Para.epsilon^2*G.*D.^2;

end

P = 1/dt + P1;

r_hat = fft2(r_hat);        num_fft = num_fft + 1;
psi_r = P .\ r_hat;
psi_r = real(ifft2(psi_r));
  
psi_H = P .\ (G .* H2); 
psi_H = real(ifft2(psi_H));

innprod_Hu = compute_ip(H, psi_r, psi_H);

w = 0.5*innprod_Hu + r;

u_old = u;
u_old_fft = u_fft;

u = 0.5*innprod_Hu * psi_H + psi_r;

t = t + dt;

u_fft = fft2(u);        num_fft = num_fft + 1;

Eu(2) = compute_En(u_fft, w);
Eu_SSAV(1) = (Eu(2))/2 + (compute_En(2*u_fft - u_old_fft, 2*w - w_old)) / 2 + ...
    sum(Para.S/2 * (u - u_old).^2, 'all')*prod(Grid.d);

E_old = Eu(1);       E = Eu(2);       t_vals(2) = t;

[F, f] = compute_F(u);

F_tilde = compute_F_tilde(u, u_fft, f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;

dt_new = dt;

while t < Time.tf

    j = j + 1;

    [u_new, w_new, gamma] = compute_unew(u, u_fft, u_old, u_old_fft, w, w_old, ...
        dt, dt_new);

    [F, f] = compute_F(u_new);

    u_new_fft = fft2(u_new);        num_fft = num_fft + 1;

    F_tilde_new = compute_F_tilde(u_new, u_new_fft, f);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute third order accurate approximation
    u_p = AM3_up(u, gamma, dt_new, F_tilde_old, F_tilde, F_tilde_new);

    rel_err = sum((u_new - u_p).^2, 'all') / sum((u_p).^2, 'all');

    A_dp = compute_Adp(rel_err, dt_new); 

    l = 1;

    while (rel_err > Para.err_tol_AM3) && (dt_new ~= Time.dt_min) && (l < 10)

        dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));
        [u_new, w_new, gamma] = compute_unew(u, u_fft, u_old, u_old_fft, ...
            w, w_old, dt, dt_new);

        [F, f] = compute_F(u_new);

        u_new_fft = fft2(u_new);        num_fft = num_fft + 1;

        F_tilde_new = compute_F_tilde(u_new, u_new_fft, f);
    
        u_p = AM3_up(u, gamma, dt_new, F_tilde_old, F_tilde, F_tilde_new);

        rel_err = sum((u_new - u_p).^2, 'all') / sum((u_p).^2, 'all');

        A_dp = compute_Adp(rel_err, dt_new);

        l = l + 1;

        if l == 10
            dt_new = 0.5*dt_new;
        end
    end

    dt = dt_new;
    dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));

    t = t + dt;

    rel_err_vals(j - 1) = sum((u_new - u_p).^2, 'all') / (prod(Grid.N));
 
    Eu(j + 1) = compute_En(u_new_fft, w_new);

    Eu_SSAV(j) = (Eu(j+1))/2 + (compute_En(2*u_new_fft - u_fft, 2*w_new - w)) / 2 + ...
        sum(Para.S/2 * (u - u_old).^2, 'all')*prod(Grid.d);

    l_vals(j) = l;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;
    u_old_fft = u_fft;              u_fft = u_new_fft;
    F_tilde_old = F_tilde;          F_tilde = F_tilde_new;

    if t >= plt(plt_idx)

        uu{plt_idx} = u;
        u_times(plt_idx) = t;
        mass(plt_idx) = (sum(u(:))*prod(Grid.d))/(prod(Grid.L));
        plt_idx = plt_idx + 1;
        
    end

    t_vals(j + 1) = t;

    if (t + dt_new) > Time.tf

        dt_new = Time.tf - t;

        if dt_new < Time.dt_min
            break
        end

    end

end

if Time.dt_max ~= Time.dt_min

    nt = sum(t_vals > 0) + 1;
    Eu = Eu(1:nt);
    Eu_SSAV = Eu_SSAV(1:nt);
    t_vals = t_vals(1:nt);
    l_vals = l_vals(1:nt);

end

Results.uu = uu;
Results.u_times = u_times;
Results.Eu = Eu;
Results.Eu_SSAV = Eu_SSAV;
Results.Em = Em;
Results.mass = mass;
Results.t_vals = t_vals;
% Results.rel_err_vals = rel_err_vals;
Results.l_vals = l_vals;
Results.num_fft = num_fft;

function u_snew = compute_u_snew(u, u_old)

u_snew = 2*u - u_old;

end

function [F, f] = compute_F(u)

f = u.^3 - Para.beta*u;

F = 0.25 * (u.^2 - Para.beta).^2;

end

% function F = compute_F(u, beta)
% 
% F = 0.25 * (u.^2 - beta).^2;
% 
% end

function w = compute_w(int_F)

w = sqrt(int_F + Para.B);

end

function H = compute_H(f, w)

H = f / w;

end

function [r, r_hat] = compute_r(u_snew_fft, u, u_old, w, w_old, ... 
    H_snew, dt, a, b, c)

r = -(b*w + c*w_old)/a + 0.5 * sum(sum((H_snew .* (b*u + c*u_old)/a)))...
    * prod(Grid.d);

rH = fft2(r .* H_snew);
% u_snew = fft2(u_snew);
num_fft = num_fft + 1;

r_tilde = -(b*u + c*u_old) - Para.S*real(ifft2(G .* u_snew_fft)) + ...
    real(ifft2(G .* rH));

m_est = (Para.alpha*Para.epsilon^2)/(a*prod(Grid.L)) * ...
    sum(r_tilde(:)) * prod(Grid.d);

r_hat = r_tilde + m_est;

end


function innprod_Hu = compute_ip(H_snew, psi_r, psi_H)
% This fuction computes the inner product of H^{*,n+1} and u^{n+1}

innprod_Hu = sum(sum(H_snew .* psi_r)) * prod(Grid.d)...
    / (1 - 0.5*sum(sum(H_snew .* psi_H)) * prod(Grid.d));

end

function E_n = compute_En(u_fft, w)

e = Para.epsilon^2/2;

% um = fft2(u - Para.m);      num_fft = num_fft + 1;
% um = u_fft - Para.m * Grid.Nx * Grid.Ny;
um_fft = u_fft;       um_fft(1,1) = 0;

v = -Grid.inv_k .* um_fft;
v(Grid.k == 0) = 0;

vx = real(ifft2(-1i*Grid.kxx.*v));         
vy = real(ifft2(-1i*Grid.kyy.*v));

dv = vx.^2 + vy.^2;


if model == 2

    du = D.*u_fft;
else

    ux = real(ifft2(-1i*Grid.kxx.*u_fft));         
    uy = real(ifft2(-1i*Grid.kyy.*u_fft));

    du = ux.^2 + uy.^2;
end

E_n = sum(sum(e *du + Para.alpha*e * dv))*prod(Grid.d) +  w^2 - Para.B;

end

function [a, b, c, gamma] = compute_dt_coeffs(dt, dt_new)

gamma = dt_new / dt;

a = (1 + 2*gamma) / (1 + gamma);            a = a / dt_new;

b = -(1 + gamma)^2 / (1 + gamma);           b = b / dt_new;

c = gamma^2 / (1 + gamma);                  c = c / dt_new;

end


function u_AM = AM3_up(u, gamma, dt_new, ...
        F_tilde_old, F_tilde, F_tilde_new)

u_AM = u + (dt_new / (6*(1 + gamma))) * ((3 + 2*gamma)*F_tilde_new + ...
    (3 + gamma)*(1 + gamma)*F_tilde - (gamma^2)*F_tilde_old);

end

function F_tilde = compute_F_tilde(u, u_fft, f)

f_fft = fft2(f);        num_fft = num_fft + 1;

F_tilde = -Para.epsilon^2 * real(ifft2(k4 .* u_fft)) - ...
    real(ifft2(k2 .* f_fft)) + ...
    Para.alpha * Para.epsilon * (u - Para.m);

end

function A_dp = compute_Adp(rel_err, dt)

A_dp = rho_s * dt * (Para.err_tol_AM3 / rel_err) ^ (1/3);

end

function [u_new, w_new, gamma] = compute_unew(u, u_fft, u_old, u_old_fft, w, w_old,...
    dt, dt_new)

u_snew = compute_u_snew(u, u_old);
    
[F, f] = compute_F(u_snew);
% f = compute_f(u_snew, Para.beta);

int_F = sum(F(:)) * prod(Grid.d);

%w(u^{*,n+1})
w_snew = compute_w(int_F);

% H^{*,n+1}
H_snew = compute_H(f, w_snew);

[a, b, c, gamma] =  compute_dt_coeffs(dt, dt_new);

[r, r_hat] = compute_r(2*u_fft - u_old_fft, u, u_old, w, w_old, H_snew, ...
    dt_new, a, b, c);

% define P for BDF2    
P = a + P1;

r_hat = fft2(r_hat);        num_fft = num_fft + 1;
psi_r = P .\ r_hat;   
psi_r = real(ifft2(psi_r));

H_fft = fft2(H_snew);       num_fft = num_fft + 1;
psi_H = P .\ (G.*H_fft);     
psi_H = real(ifft2(psi_H));

innprod_Hu = compute_ip(H_snew, psi_r, psi_H);

w_new = 0.5*innprod_Hu + r;      

u_new = 0.5*innprod_Hu * psi_H + psi_r;

end

end