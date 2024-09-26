function [u, Eu, Em, mass, m_est_vals, t_vals, dt_vals]...
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

err_tol = 1E-5;        % error tolerance

dt = Time.dt_min;


if model == 2

    D = (-Grid.k.^2 + 1);    % PFC
    
else
    
    D = -1i*Grid.k;          % AC, CH and OK
end

% For all 3 models the operator G is the same
if model == 3
    
    G = -1;             % AC
    
else
    G = -Grid.k.^2;          % PFC, CH and OK
    
end
    
v = fftn(u);

% define number of time steps and time steps for plotting purposes
nplt = floor( (Time.tf / 100) / Time.dt_min);
nmax = round(Time.tf / Time.dt_min);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = compute_F(u, Para.beta);         int_F = sum(F(:)) * Grid.dx*Grid.dy;
f = compute_f(u, Para.beta);

% define the following as w_old b/c we will need this in the for loop
w_old = compute_w(int_F, Para.B);   

t = 0;


mass = zeros(1, 101);
mass(1) = (sum(u(:))*Grid.dx*Grid.dy)/(Grid.Lx*Grid.Ly);

m_est_vals = zeros(1, 100);

Em = sum(sum((compute_F(Para.m*ones(size(u)), Para.beta))))*Grid.dx*Grid.dy;

Eu = zeros(1, nmax);
t_vals = zeros(1, nmax);

Eu(1) = compute_En(u, w_old, Para, Grid.dx, Grid.dy, Grid.inv_k,...
    Grid.k, Em, G, D, model);

H = compute_H(f, w_old);

r = -0.5 * sum(sum(H.*u)) * Grid.dx*Grid.dy + w_old;
H2 = fftn(H);

r_tilde = u/dt - Para.S*real(ifftn(reshape(G .* v(:),Grid.Nx,Grid.Ny))) - ...
    r*real(ifftn(reshape(G .* H2(:),Grid.Nx,Grid.Ny)));

r_hat = (Para.alpha*Para.epsilon^2*dt)/(Grid.Lx*Grid.Ly) * ...
    sum(r_tilde(:))*Grid.dx*Grid.dy + r_tilde;


m_est_vals(1) = (Para.alpha*Para.epsilon^2*dt)/(Grid.Lx*Grid.Ly) * ...
    sum(r_tilde(:))*Grid.dx*Grid.dy;


if model == 2

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G - Para.epsilon^2*G.*D.^2;

else

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G + Para.epsilon^2*G.*D.^2;

end

P = spdiags(1/dt + P1, 0, Grid.Nx*Grid.Ny, Grid.Nx*Grid.Ny);

r_hat = fftn(r_hat);
psi_r = P \ r_hat(:);
psi_r = real(ifftn(reshape(psi_r, Grid.Nx, Grid.Ny)));
  
psi_H = P \ (G .* H2(:)); 
psi_H = real(ifftn(reshape(psi_H, Grid.Nx, Grid.Ny)));

innprod_Hu = compute_ip(H, psi_r, psi_H, Grid.dx, Grid.dy);

w = 0.5*innprod_Hu + r;

u_old = u;

u = 0.5*innprod_Hu * psi_H + psi_r;

dt_vals = dt;

t = t + dt;


Eu(2) = compute_En(u, w, Para, Grid.dx, Grid.dy, Grid.inv_k, Grid.k, Em, G, D, model);
E_old = Eu(1);       E = Eu(2);       t_vals(2) = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;


while t < Time.tf

    j = j + 1;

    dt_new = dt;
    
    t = t + dt_new;

    [u_new, w_new, m_est, gamma] = compute_unew(u, u_old, w, w_old, ...
        Grid, dt, dt_new, Para, G, D,...
        model, P1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_new = compute_En(u, w, Para, Grid.dx, Grid.dy, Grid.inv_k, Grid.k,...
        Em, G, D, model);

    ut = sum(sum(((u_new - u_old) / (dt_new + dt)).^2))*Grid.dx * Grid.dy;

    Et = (E_new - E_old) / (dt_new + dt);

    Ett = ((dt/dt_new)*E_new - (dt/dt_new + 1)*E + E_old) / ...
        (0.5 * (dt^2 + dt*dt_new));

    rel_err = 0.5*dt_new*abs(Et) + dt_new^2 * (Para.S/2 * ut + ...
        (Ett < 0)*abs(Ett));
    rel_err = rel_err / abs(E_new);

    adap = 0.5*dt*(err_tol / rel_err);


    dt_new = max(Time.dt_min, min(adap, Time.dt_max));

    while (rel_err > err_tol) && (dt_new ~= Time.dt_min)

        [u_new, w_new, m_est, gamma] = compute_unew(u, u_old, w, w_old, ...
        Grid, dt, dt_new, Para, G, D,...
        model, P1);

        Eu(j + 1) = compute_En(u, w, Para, Grid.dx, Grid.dy, Grid.inv_k, Grid.k, Em, G, D, model);
        
        E_new = Eu(j + 1);
    
        ut = sum(sum(((u_new - u_old) / (dt_new + dt)).^2))*Grid.dx * Grid.dy;
    
        Et = (E_new - E_old) / (dt_new + dt);
    
        Ett = ((dt/dt_new)*E_new - (dt/dt_new + 1)*E + E_old) / ...
            (0.5 * (dt^2 + dt*dt_new));
    
        rel_err = 0.5*dt_new*abs(Et) + dt_new^2 * (Para.S/2 * ut + ...
            (Ett < 0)*abs(Ett));
        rel_err = rel_err / abs(E_new);

        adap = 0.5*dt*(err_tol / rel_err);

        dt_new = max(Time.dt_min, min(adap, Time.dt_max));
    end

    % compute new time step
    dt = dt_new;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    Eu(j + 1) = E_new;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;
    E_old = E;                      E = E_new;

    t_vals(j + 1) = t;

end

nt = sum(t_vals > 0) + 1;

Em = Em / (Grid.Lx*Grid.Ly);
Eu = Eu(1:nt) ./ (Grid.Lx*Grid.Ly);
t_vals = t_vals(1:nt);

end


function u_snew = compute_u_snew(u, u_old)

u_snew = 2*u - u_old;

end

function f = compute_f(u, beta)

f = u.^3 - beta*u;

end

function F = compute_F(u, beta)

F = 0.25 * (u.^2 - beta).^2;

end

function w = compute_w(int_F, B)

w = sqrt(int_F + B);

end

function H = compute_H(f, w)

H = f / w;

end

function [r, r_hat, m_est] = compute_r(u_snew, u, u_old, w, w_old, ... 
    H_snew, S, Grid, dt, alpha, epsilon, a, b, c, G)

r = -(b*w + c*w_old)/a + 0.5 * sum(sum((H_snew .* (b*u + c*u_old)/a)))...
    * Grid.dx*Grid.dy;

rH = fftn(r .* H_snew);
u_snew = fftn(u_snew);

r_tilde = -(b*u + c*u_old) - S*real(ifftn(reshape(G .* u_snew(:), ...
    Grid.Nx, Grid.Ny))) + ...
    real(ifftn( reshape(G .* rH(:), Grid.Nx, Grid.Ny)));

m_est = (alpha*epsilon^2)/(a*Grid.Lx*Grid.Ly) * ...
    sum(r_tilde(:)) * Grid.dx*Grid.dy;

r_hat = r_tilde + m_est;

end


function innprod_Hu = compute_ip(H_snew, psi_r, psi_H, dx, dy)
% This fuction computes the inner product of H^{*,n+1} and u^{n+1}

innprod_Hu = sum(sum(H_snew .* psi_r)) * dx*dy...
    / (1 - 0.5*sum(sum(H_snew .* psi_H)) * dx*dy);

end

function E_n = compute_En(u, w, Para, dx, dy, inv_k, k, Em, G, D, model)

e = Para.epsilon^2/2;

um = fftn(u - Para.m);

v = -inv_k .* um(:);
v(spdiags(k) == 0) = 0;

u = fftn(u);

E_n = sum(sum(e * real(ifftn(D.*u(:))).^2 + ...
    Para.alpha*e * real(ifftn(D.*v)).^2))*dx*dy + ...
    w^2 - Para.B;

end

function [a, b, c, gamma] = compute_dt_coeffs(dt, dt_new)

gamma = dt_new / dt;

a = (1 + 2*gamma) / (1 + gamma);            a = a / dt_new;

b = -(1 + gamma)^2 / (1 + gamma);           b = b / dt_new;

c = gamma^2 / (1 + gamma);                  c = c / dt_new;

end


function u_p = AM3_up(u_new, u, u_old, epsilon, alpha, gamma, k, ...
    dt_new, m, Nx, Ny, beta)

unew_fft = fftn(u_new);     fnew_fft = fftn(compute_f(u_new, beta));   

u_fft = fftn(u);            f_fft = fftn(compute_f(u, beta));                        

uold_fft = fftn(u_old);     fold_fft = fftn(compute_f(u_old, beta));

F_tilde_unew = compute_F_tilde(unew_fft, u_new, fnew_fft, k, epsilon, ...
    alpha, m, Nx, Ny);

F_tilde_u = compute_F_tilde(u_fft, u, f_fft, k, epsilon, ...
    alpha, m, Nx, Ny);

F_tilde_uold = compute_F_tilde(uold_fft, u_old, fold_fft, k, epsilon, ...
    alpha, m, Nx, Ny);

u_p = u + (dt_new / 6) * (((3 + 2*gamma)/(1 + gamma))*F_tilde_unew + ...
    (3 + gamma)*F_tilde_u - (gamma^2 / (1 + gamma))*F_tilde_uold);

end

function F_tilde = compute_F_tilde(u_fft, u, f_fft, k, epsilon, alpha, m, Nx, Ny)

F_tilde = -epsilon^2 * real(ifftn(reshape(k.^4 .* u_fft(:), Nx, Ny))) + ...
    real(ifftn(reshape(-k.^2 .* f_fft(:), Nx, Ny))) + ...
    alpha * epsilon * (u - m);

end

function A_dp = compute_A_dp(rho_s, err_tol, err, dt)

A_dp = rho_s * dt * (err_tol / err) ^ (1/3);

end

function [u_new, w_new, m_est, gamma] = compute_unew(u, u_old, w, w_old,...
    Grid, dt, dt_new, Para, G, D, model, P1)

u_snew = compute_u_snew(u, u_old);
    
F = compute_F(u_snew, Para.beta);
f = compute_f(u_snew, Para.beta);

int_F = sum(F(:)) * Grid.dx*Grid.dy;

%w(u^{*,n+1})
w_snew = compute_w(int_F, Para.B);

% H^{*,n+1}
H_snew = compute_H(f, w_snew);

[a, b, c, gamma] =  compute_dt_coeffs(dt_new, dt);

[r, r_hat, m_est] = compute_r(u_snew, u, u_old, w, w_old, H_snew, Para.S, ...
    Grid, dt_new, Para.alpha, Para.epsilon, a, b, c, G);

% define P for BDF2    
P = spdiags(a + P1, 0, Grid.Nx*Grid.Ny, Grid.Nx*Grid.Ny);

r_hat = fftn(r_hat);
psi_r = P \ r_hat(:);           psi_r = reshape(psi_r, Grid.Nx, Grid.Ny);
psi_r = real(ifftn(psi_r));

H2 = fftn(H_snew);
psi_H = P \ (G.*H2(:));     psi_H = reshape(psi_H, Grid.Nx, Grid.Ny);
psi_H = real(ifftn(psi_H));

innprod_Hu = compute_ip(H_snew, psi_r, psi_H, Grid.dx, Grid.dy);

w_new = 0.5*innprod_Hu + r;      

u_new = 0.5*innprod_Hu * psi_H + psi_r;

end