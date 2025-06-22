function Results = SSAV(Grid, Time, Para, u, model, dim)

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

num_fft = 0;

dt = Time.dt_min;

% k2 = Grid.kxx.^2 + Grid.kyy.^2;
k4 = Grid.k2.^2;         k6 = Grid.k2.^3.;


if model == 2

    D = (-Grid.k2 + 1);          % PFC
    
else
    
    D = -1i*Grid.k;          % AC, CH and OK
end

% For all 3 models the operator G is the same
if model == 3
    
    G = -1;             % AC
    
else
    G = -Grid.k2;          % PFC, CH and OK
    
end    

nmax = round(Time.tf / Time.dt_min);

plt = linspace(0, Time.tf, 101);
plt_idx = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F, f] = SSAV_helpers.compute_F(u , Para.beta);         
int_F = sum(F, 'all') * prod(Grid.d);

% define the following as w_old b/c we will need this in the for loop
w_old = SSAV_helpers.compute_w(int_F, Para.B);   

t = 0;

mass = zeros(1, 101);
mass(1) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));

Em = prod(Grid.L) * SSAV_helpers.compute_F(Para.m, Para.beta);

Eu = zeros(1, nmax); 
Eu_SSAV = zeros(1, nmax);
Et_vals = zeros(1, nmax);
t_vals = zeros(1, nmax);
l_vals = zeros(1, nmax);
% uu = cell(101, 1);
% uu = cell(nmax, 1);
% uu{1} = u;
% u_times = zeros(1, 101);

M_vals = zeros(1, nmax); 


u_fft = fftn(u);         num_fft = num_fft + 1;
Eu(1) = SSAV_helpers.compute_En(u_fft, w_old, Para, Grid, D, model, dim);

H = SSAV_helpers.compute_H(f, w_old);

r = -0.5 * sum(H.*u, 'all') * prod(Grid.d) + w_old;
H2 = fftn(H);           num_fft = num_fft + 1;


r_tilde = u/dt - Para.S*real(ifftn(G .* u_fft)) - r*real(ifftn(G .* H2));

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
psi_r = real(ifftn(psi_r));

  
psi_H = P .\ (G .* H2);
psi_H = real(ifftn(psi_H));


innprod_Hu = SSAV_helpers.compute_ip(H, psi_r, psi_H, Grid);

w = 0.5*innprod_Hu + r;

u_old = u;
u_old_fft = u_fft;

u = 0.5*innprod_Hu * psi_H + psi_r;
u_fft = fftn(u);        num_fft = num_fft + 1;

uu{2} = u;

t = t + dt;


Eu(2) = SSAV_helpers.compute_En(u_fft, w, Para, Grid, D, model, dim);

Eu_SSAV(1) = (Eu(2))/2 + (SSAV_helpers.compute_En(2*u_fft - u_old_fft, ...
    2*w - w_old, Para, Grid, D, model, dim)) / 2 + ...
    sum(Para.S/2 * (u - u_old).^2, 'all')*prod(Grid.d);

Et_vals(1) = (Eu(2) - Eu(1)) / dt;

t_vals(2) = t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;

dt_new = dt;

while t < Time.tf

    j = j + 1;

    [u_new, w_new] = SSAV_helpers.compute_unew(u, u_fft, u_old, ...
        u_old_fft, w, w_old, dt, dt_new, Para, Grid, G, P1);

    u_new_fft = fftn(u_new);        num_fft = num_fft + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E_new = SSAV_helpers.compute_En(u_new_fft, w_new, Para, Grid, D, model, dim);

    E_mod = E_new/2 + (SSAV_helpers.compute_En(2*u_new_fft - u_fft, ...
        2*w_new - w, Para, Grid, D, model, dim)) / 2 + ...
        sum(Para.S/2 * (u_new - u).^2, 'all')*prod(Grid.d);

    E_t = (E_new - Eu(j)) / dt_new;

    Et_vals(j) = E_t;

    M_vals(j) = sqrt(1 + 2* Para.sigma /dt * (E_new - E_mod)^2);

    if Time.adap < 3

        if Time.adap == 1
            
            rel_err = SSAV_helpers.compute_rel_err(u, u_p);

        else

            rel_err = abs(E_new - E_mod);

        end
    else
        dt = dt_new;

        if Time.adap == 3
            
            dt_new = max(Time.dt_min, Time.dt_max ...
                / sqrt(1 + Para.sigma * abs(E_t)^2));

        else
           
            dt_new = max(Time.dt_min, Time.dt_max / ...
                sqrt(1 + Para.sigma*(2*abs(E_new - E_mod)/dt)^(Para.p)));
        end
    end



    % rel_err = abs(E_new - E_mod);
    % 
    % A_dp = compute_Adp(rel_err, dt_new);
    % 
    % l = 1;
    % 
    % while ((rel_err > Para.err_tol) && (dt_new ~= Time.dt_min)) && (l < 10)
    % 
    %     dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));
    % 
    %     [u_new, w_new] = compute_unew(u, u_fft, u_old, u_old_fft, w, w_old, ...
    %         dt, dt_new);
    % 
    %     u_new_fft = fftn(u_new);        num_fft = num_fft + 1;
    % 
    %     E_new = compute_En(u_new_fft, w_new);
    % 
    %     E_mod = E_new/2 + (compute_En(2*u_new_fft - u_fft, 2*w_new - w)) / 2 + ...
    %         sum(Para.S/2 * (u_new - u).^2, 'all')*Grid.dx*Grid.dy;
    % 
    %     rel_err = abs(E_new - E_mod);
    % 
    %     A_dp = compute_Adp(rel_err, dt_new);
    % 
    %     l = l + 1;
    % 
    %     if l == 10
    %         dt_new = 0.5*dt_new;
    %     end
    % end

    % compute new time step
    dt = dt_new;
    % dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));
    % dt_new = max(Time.dt_min, Time.dt_max / sqrt(1 + Para.sigma * abs(E_t)^2));
    dt_new = max(Time.dt_min, Time.dt_max / ...
        sqrt(1 + Para.sigma*(2*abs(E_new - E_mod)/dt)^(Para.p)));

    t = t + dt;
 
    Eu(j + 1) = E_new;

    Eu_SSAV(j) = E_mod;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;
    u_old_fft = u_fft;              u_fft = u_new_fft;

    % if t >= plt(plt_idx)
    % 
    %     uu{plt_idx} = u;
    %     u_times(plt_idx) = t;
    %     mass(plt_idx) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));
    %     plt_idx = plt_idx + 1;
    % 
    % end

    t_vals(j + 1) = t;

    if (t + dt_new) > Time.tf && Time.dt_max ~= Time.dt_min

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
    Et_vals = Et_vals(1:nt);
    l_vals = l_vals(1:nt);
    t_vals = t_vals(1:nt);

end

% Results.uu = uu(~cellfun('isempty', uu));
% Results.uu = uu;
% Results.u_times = u_times;
Results.u = u;
Results.Eu = Eu;
Results.Eu_SSAV = Eu_SSAV;
Results.Et_vals = Et_vals;
Results.Em = Em;
Results.mass = mass;
Results.t_vals = t_vals;
Results.l_vals = l_vals;
Results.num_fft = num_fft;
Results.M_vals = M_vals;

end