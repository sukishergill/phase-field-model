function Results = SSAV(Grid, Time, Para, u, model, dim, save_u, plt_save)

% This function solves the following PDE
%
%           u_t = G(-epsilon*Du + f(u) + alpha*epsilon v)
%
% where f(u) = u^3 - beta * u
%
% We implemented the stabilized scalar auxiliary variable (SSAV) method.
%
% We use BDF2 with variable step size for the time discretization and
% Fourier method for the spatial discretization.
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
% Input variables:
%     - Grid: Field with variables for the grid discretization
%     - Time: Field with variables for the time discretization
%     - Para: Field with parameter values store
%     - u: initial condition
%     - model: indicates which eqn to solve
%               1. Ohta-Kawasaki or Cahn-Hilliard
%               2. Phase-field-crystals
%               3. Allen-Cahn
%     - dim: dimensions (1, 2, or 3)
%     - save_u: indicate if u values are stored or not

t = Time.t0;

num_fft = 0;        % start counter for number of FFTs
reject_steps = 0;   % counter for number of rejected steps for time stepper

dt = Time.dt_min;   % initialize the time step

% Define operator D
if model == 2

    D = (-Grid.k2 + 1);          % PFC
    
else
    
    D = -1i*Grid.k;          % AC, CH and OK
end

% Define operator G
if model == 3
    
    G = -1;                % AC
    
else
    G = -Grid.k2;          % PFC, CH and OK
    
end    

nmax = round(Time.tf / Time.dt_min);

plt = linspace(0, Time.tf, plt_save + 1);
plt_idx = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f] = SSAV_helpers.compute_F(u , Para.beta);         
int_F = sum(F, 'all') * prod(Grid.d);

% define the following as w_old b/c we will need this in the for loop
w_old = SSAV_helpers.compute_w(int_F, Para.B);   

mass = zeros(1, 101);
mass(1) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));

Em = prod(Grid.L) * SSAV_helpers.compute_F(Para.m, Para.beta);

Eu = zeros(1, nmax); 
Eu_SSAV = zeros(1, nmax);
Eu_sym = zeros(1, nmax);
Et_vals = zeros(1, nmax);
t_vals = zeros(1, nmax);
u_times = zeros(1, plt_save + 1);
dt_vals = zeros(1, nmax);
dt_idx = zeros(1, nmax);
dt_prop_true = zeros(1, nmax);
w_vals = zeros(1, nmax);
l_vals = zeros(1,nmax);

w_vals(1) = w_old;


if save_u == 1
    uu = cell(nmax, 1);
    uu{1} = u;
end

if save_u == 0
    uu = cell(plt_save, 1);
    uu{1} = u;

    for i = 2:plt_save + 1
        uu{i} = zeros(size(uu{1}));
    end
end

M_vals = zeros(1, nmax); 


u_fft = fftn(u);         num_fft = num_fft + 1;
Eu(1) = SSAV_helpers.compute_En(u_fft, w_old, Para, Grid, D, model, dim);

H = SSAV_helpers.compute_H(f, w_old);

r = -0.5 * sum(H.*u, 'all') * prod(Grid.d) + w_old;
H2 = fftn(H);           num_fft = num_fft + 1;

forcing = SSAV_helpers.compute_forcing(Para.epsilon, Grid.xx, Grid.yy, dt);

r_tilde = u/dt - Para.S*(ifftn(G .* u_fft, 'symmetric')) - ...
    r*(ifftn(G .* H2, 'symmetric')) + forcing;

r_hat = (Para.alpha*Para.epsilon^2*dt)/(prod(Grid.L)) * ...
    sum(r_tilde, 'all')*prod(Grid.d) + r_tilde;


if model == 2

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G - Para.epsilon^2*G.*D.^2;

else

    P1 = Para.alpha*Para.epsilon^2 - Para.S*G + Para.epsilon^2*G.*D.^2;

end


P = 1/dt + P1;

r_hat = fftn(r_hat);        num_fft = num_fft + 1;
psi_r = P .\ r_hat;
psi_r = ifftn(psi_r, 'symmetric');

  
psi_H = P .\ (G .* H2);
psi_H = ifftn(psi_H, 'symmetric');


innprod_Hu = SSAV_helpers.compute_ip(H, psi_r, psi_H, Grid);

w = 0.5*innprod_Hu + r;
w_vals(2) = w;

u_old = u;
u_old_fft = u_fft;

u = 0.5*innprod_Hu * psi_H + psi_r;
u_fft = fftn(u);        num_fft = num_fft + 1;


if Time.adap == 1

    [F_tilde_old, num_fft] = SSAV_helpers.compute_F_tilde(u_old, ...
        u_old_fft, f, Grid, Para, num_fft);
    
    [~, f] = SSAV_helpers.compute_F(u, Para.beta);
    [F_tilde, num_fft] = SSAV_helpers.compute_F_tilde(u, u_fft, f, Grid,...
        Para, num_fft);

end

t = t + dt;

Eu(2) = SSAV_helpers.compute_En(u_fft, w, Para, Grid, D, model, dim);

Eu_SSAV(1) = (Eu(2))/2 + (SSAV_helpers.compute_En(2*u_fft - u_old_fft, ...
    2*w - w_old, Para, Grid, D, model, dim)) / 2 + ...
    sum(Para.S/2 * (u - u_old).^2, 'all')*prod(Grid.d);

Eu_sym(1) = (Eu(2))/2 + (SSAV_helpers.compute_En(2*u_fft - u_old_fft, ...
    2*w - w_old, Para, Grid, D, model, dim)) / 2;

E_t = (Eu(2) - Eu(1)) / dt;

if Time.adap == 5

    dt_1 = Para.err_tol(1) / (abs(E_t) + Para.delta);

    dt_2 = Para.err_tol(2) / sqrt(abs(E_t)+ Para.delta^2);

    [dt_prop, dt_idx(1)] = min([dt_1, dt_2]);
    % dt_prop = dt_1;

    [min_dt, dt_prop_true(1)] = min([Time.dt_max, 2*dt, dt_prop]);

    dt_new = max([Time.dt_min, min_dt, 0.1*dt]);

else
    dt_new = dt;
end

Et_vals(2) = E_t;

t_vals(2) = t;

E_old = Eu(1);      E = Eu(2);

mass(2) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;

% dt_new = dt;
% dt_new = 1.10*dt;
dt_vals(1) = dt_new;

while t < Time.tf

    gamma = dt_new / dt;

    j = j + 1;

    dt_vals(j) = dt_new;

    [u_new, w_new, num_fft] = SSAV_helpers.compute_unew(u, u_fft, u_old, ...
        u_old_fft, w, w_old, dt, dt_new, Para, Grid, G, P1, num_fft, t);

    u_new_fft = fftn(u_new);        num_fft = num_fft + 1;

    w_vals(j+1) = w_new;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E_new = SSAV_helpers.compute_En(u_new_fft, w_new, Para, Grid, D, model, dim);

    Eu_sym(j + 1) = 0.5*(E_new + SSAV_helpers.compute_En(2*u_new_fft - u_fft, ...
            2*w_new - w, Para, Grid, D, model, dim));

    E_mod = Eu_sym(j+1) + sum(Para.S/2 * (u_new - u).^2, 'all')*prod(Grid.d);

    E_t = (Eu_sym(j+1) - Eu_sym(j)) / dt_new;
    Et_vals(j) = E_t;

    l = 0;
    if Time.adap == 5
    while E_t > 1e-8 && dt_new > Time.dt_min && l <= 5

        dt_floor = min(Time.dt_min, 0.1*dt);

        reject_steps = reject_steps + 1;

        % dt_new = max(gamma * dt_new / 4, Time.dt_min);
        % dt_new = max(0.5*dt_new, Time.dt_min);
        dt_new = max(0.5*dt_new, dt_floor);
        gamma = dt_new / dt;

        dt_vals(j) = dt_new;

        [u_new, w_new, num_fft] = SSAV_helpers.compute_unew(u, u_fft, u_old, ...
            u_old_fft, w, w_old, dt, dt_new, Para, Grid, G, P1, num_fft, t);

        u_new_fft = fftn(u_new);        num_fft = num_fft + 1;

        w_vals(j+1) = w_new;

        E_new = SSAV_helpers.compute_En(u_new_fft, w_new, Para, Grid, D, model, dim);

        Eu_sym(j + 1) = 0.5*(E_new + SSAV_helpers.compute_En(2*u_new_fft - u_fft, ...
            2*w_new - w, Para, Grid, D, model, dim));

        E_mod = Eu_sym(j+1) + sum(Para.S/2 * ...
            (u_new - u).^2, 'all')*prod(Grid.d);

        E_t = (Eu_sym(j+1) - Eu_sym(j)) / dt_new;

        l = l + 1;

    end
    end

    l_vals(j+1) = l;
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Time.dt_max == Time.dt_min
        dt_new = Time.dt_min;
    
    elseif Time.adap == 1

        [A_dp, rel_err, F_tilde_new, num_fft] = SSAV_helpers.compute_Adp(u_new, ...
            u, u_new_fft, F_tilde_old, F_tilde, dt_new, gamma, Para, ...
            Grid, num_fft);

        l = 1;

        while (rel_err > Para.err_tol_AM3) && (dt_new ~= Time.dt_min) && (l < 10)
            
            dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));
            gamma = dt_new / dt;

            [u_new, w_new, num_fft] = SSAV_helpers.compute_unew(u, u_fft, ...
                u_old, u_old_fft, w, w_old, dt, dt_new, Para, ...
                Grid, G, P1, num_fft, t);

            u_new_fft = fftn(u_new);        num_fft = num_fft + 1;

            [A_dp, rel_err, F_tilde_new, num_fft] = ...
                SSAV_helpers.compute_Adp(u_new, u, u_new_fft, F_tilde_old, ...
                F_tilde, dt, dt_new, Para, Grid, num_fft);

            l = l + 1;

            if l == 10
                dt_new = 0.5*dt_new;
            end

        end

        dt = dt_new;
        dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));

        F_tilde_old = F_tilde;      F_tilde = F_tilde_new;

    else

        dt = dt_new;

        if Time.adap == 3
            
            dt_new = max(Time.dt_min, Time.dt_max ...
                / sqrt(1 + Para.sigma^2 * abs(E_t)^2));

        elseif Time.adap == 4

            dt_new = max(Time.dt_min, Time.dt_max/M_vals(j));

        else
           

            dt_1 = Para.err_tol(1) / (abs(E_t) + Para.delta);

            dt_2 = Para.err_tol(2) / sqrt(abs(E_t)+ Para.delta^2);

            [dt_prop, dt_idx(j)] = min([dt_1, dt_2]);
            % dt_prop = dt_1;
        
            % [min_dt, dt_prop_true(j)] = min([Time.dt_max, 2*dt, dt_prop]);

            min_dt = min(Time.dt_max, max(Time.dt_min, dt_prop));
            dt_new = min(2*dt, max(0.1*dt, min_dt));

            % dt_new = max([Time.dt_min, min_dt, 0.1*dt]);
        end
    end

    t = t + dt;
 
    Eu(j + 1) = E_new;
    Eu_SSAV(j + 1) = E_mod;
    Et_vals(j + 1) = E_t;
    mass(j) = (sum(u_new, 'all')*prod(Grid.d))/(prod(Grid.L));

    if save_u == 1 
        uu{j+1} = u_new;
    end


    E_old = E;                      E = E_new;
    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;
    u_old_fft = u_fft;              u_fft = u_new_fft;

    mass(j + 1) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));

    if t >= plt(plt_idx)

        uu{plt_idx} = u;
        u_times(plt_idx) = t;
        plt_idx = plt_idx + 1;

    end

    t_vals(j + 1) = t;

    if (t + dt_new) > Time.tf && Time.dt_max ~= Time.dt_min

        dt_new = Time.tf - t;

        if dt_new < Time.dt_min
            uu{end} = u_new;
            break
        end

    end


end

if Time.dt_max ~= Time.dt_min

    nt = sum(t_vals > 0) + 1;
    Eu = Eu(1:nt);
    Eu_SSAV = Eu_SSAV(1:nt);
    Eu_sym = Eu_sym(1:nt);
    Et_vals = Et_vals(1:nt);
    t_vals = t_vals(1:nt);
    l_vals = l_vals(1:nt);
    dt_idx = dt_idx(1:nt-1);
    dt_vals = dt_vals(1:nt-1);
    dt_prop_true = dt_prop_true(1:nt-1);
    mass = mass(1:nt);
    w_vals = w_vals(1:nt);

end

Results.uu = uu;
Results.u_times = u_times;
Results.Eu = Eu;
Results.Eu_SSAV = Eu_SSAV;
Results.Eu_sym = Eu_sym;
Results.Et_vals = Et_vals;
Results.Em = Em;
Results.mass = mass;
Results.t_vals = t_vals;
Results.l_vals = l_vals;
Results.num_fft = num_fft;
Results.reject_steps = reject_steps;
Results.dt_idx = dt_idx;
Results.dt_vals = dt_vals;
Results.dt_prop_true = dt_prop_true;
Results.w_vals = w_vals;

end