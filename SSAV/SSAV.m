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
Ett_vals = zeros(1, nmax);
t_vals = zeros(1, nmax);
u_t_vals = zeros(1, 9);
dt_vals = zeros(1, nmax);
dt_idx = zeros(1, nmax);
dt_prop_true = zeros(1, nmax);
% l_vals = zeros(1, nmax);
w_vals = zeros(1, nmax);
wt_vals = zeros(1, nmax);
wtt_vals = zeros(1, nmax);
wttt_vals = zeros(1, nmax);

w_vals(1) = w_old;

u_ttt_pts = cell(4, 1);
w_ttt_pts = cell(4, 1);
w_ttt_pts{1} = w_old;
w_ttt_vals = zeros(3,1);

if save_u == 1
    uu = cell(nmax, 1);
    uu{1} = u;
end
% uu{1} = u;
u_ttt_pts{1} = u;
u_ttt_vals = cell(3, 1);
% u_times = zeros(1, 101);

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

w_ttt_pts{2} = w;
u_ttt_pts{2} = u;

if Time.adap == 1

    [F_tilde_old, num_fft] = SSAV_helpers.compute_F_tilde(u_old, ...
        u_old_fft, f, Grid, Para, num_fft);
    
    [~, f] = SSAV_helpers.compute_F(u, Para.beta);
    [F_tilde, num_fft] = SSAV_helpers.compute_F_tilde(u, u_fft, f, Grid,...
        Para, num_fft);

end
ut_old = (u - u_old) / dt;
wt_old = (w - w_old) / dt;

wt_vals(2) = wt_old;


% uu{2} = u;

t = t + dt;


Eu(2) = SSAV_helpers.compute_En(u_fft, w, Para, Grid, D, model, dim);

Eu_SSAV(1) = (Eu(2))/2 + (SSAV_helpers.compute_En(2*u_fft - u_old_fft, ...
    2*w - w_old, Para, Grid, D, model, dim)) / 2 + ...
    sum(Para.S/2 * (u - u_old).^2, 'all')*prod(Grid.d);

Eu_sym(1) = (Eu(2))/2 + (SSAV_helpers.compute_En(2*u_fft - u_old_fft, ...
    2*w - w_old, Para, Grid, D, model, dim)) / 2;

Et_vals(1) = (Eu(2) - Eu(1)) / dt;

t_vals(2) = t;

E_old = Eu(1);      E = Eu(2);

mass(2) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;

dt_new = dt;
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

    % if j == 2
    %     w_ttt_pts{3} = w_new;
    %     u_ttt_pts{3} = u_new;
    % 
    % elseif j == 3
    % 
    %     w_ttt_pts{4} = w_new;
    %     u_ttt_pts{4} = u_new;
    % 
    % else
    %     w_ttt_pts{1} = [];
    %     w_ttt_pts = w_ttt_pts(~cellfun('isempty', w_ttt_pts));
    %     w_ttt_pts{4} = w_new;
    % 
    %     u_ttt_pts{1} = [];
    %     u_ttt_pts = u_ttt_pts(~cellfun('isempty', u_ttt_pts));
    %     u_ttt_pts{4} = u_new;
    % end

    if j == 2
        % w_ttt_pts{3} = w_new;
        % u_ttt_pts{3} = u_new;

        ut = (u_new - u) / dt_new;
        utt = (ut - ut_old) / (dt_new + dt);

        ut_old = ut;

        wt = (w_new - w) / dt_new;
        wt_vals(j) = wt;

        wtt_old = (wt - wt_old) / (dt_new + dt);
        wtt_vals(j) = wtt_old;

        wt_old = wt;

    else%if j == 3
        % w_ttt_pts{4} = w_new;
        % u_ttt_pts{4} = u_new;
        utt_old = utt;

        ut = (u_new - u) / dt_new;
        utt = (ut - ut_old) / (dt_new + dt);

        wt = (w_new - w) / dt_new;
        wt_vals(j) = wt;

        wtt = (wt - wt_old) / (dt_new + dt);
        wtt_vals(j) = wtt;

        u_ttt = (utt - utt_old) / (sum(dt_vals(j-2:j)));
        w_ttt = (wtt - wtt_old) / (sum(dt_vals(j-2:j)));

        wttt_vals(j) = w_ttt;

        ut_old = ut;        wt_old = wt;
        utt_old = utt;      wtt_old = wtt;
    % else
    %     w_ttt_pts{1} = [];
    %     w_ttt_pts = w_ttt_pts(~cellfun('isempty', w_ttt_pts));
    %     w_ttt_pts{4} = w_new;
    % 
    %     u_ttt_pts{1} = [];
    %     u_ttt_pts = u_ttt_pts(~cellfun('isempty', u_ttt_pts));
    %     u_ttt_pts{4} = u_new;


    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E_new = SSAV_helpers.compute_En(u_new_fft, w_new, Para, Grid, D, model, dim);

    E_mod = E_new/2 + (SSAV_helpers.compute_En(2*u_new_fft - u_fft, ...
        2*w_new - w, Para, Grid, D, model, dim)) / 2 + ...
        sum(Para.S/2 * (u_new - u).^2, 'all')*prod(Grid.d);

    drift = abs(E_new - E_mod);

    E_t = (E_new - Eu(j)) / dt_new;
    Et_vals(j) = E_t;
    % E_tt = (E_new - 2*E + E_old)/(dt_new^2);

    M_vals(j) = sqrt(1 + Para.sigma*(2*abs(drift)/dt_new));
    % M_vals(j) = (E_new - E_mod)/dt_new;


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
                / sqrt(1 + Para.sigma * abs(E_t)^2));

        elseif Time.adap == 4

            dt_new = max(Time.dt_min, Time.dt_max/M_vals(j));

        else
           

            % E_tt = SSAV_helpers.compute_deriv2(E_new, E, E_old, ...
            %     dt_new, dt);
            % E_tt = (Et_vals(j) - Et_vals(j-1)) / (dt_new + dt);
            E_tt = ((E_new - E)/dt_new - (E - E_old)/dt) / (dt_new + dt);
            Ett_vals(j) = E_tt;

            % u_tt = SSAV_helpers.compute_deriv2(u_new_fft, u_fft, u_old_fft, ...
            %     dt_new, dt);
            % u_tt = ((u_new - u)/dt_new - (u - u_old)/dt) / (dt_new + dt);

            Du = ifftn(D.^2 .* u_new_fft, 'symmetric');
            
            delu = ifftn(G .* utt, 'symmetric');

            dt_1 = 2 * Para.err_tol(1) / (abs(E_t) + Para.delta);

            dt_2 = sqrt(2 * Para.err_tol(2) / (abs(E_t)+ Para.delta));

            dt_3 = sqrt(2 * Para.err_tol(3) / (abs(E_tt) + Para.delta));

            dt_5 = (gamma*Para.err_tol(5) / ((gamma + 1) * ...
                abs(sum(Du .* delu, 'all')*prod(Grid.d)/prod(Grid.N)) + Para.delta)) ^ (1/3);

            if isreal(dt_5) == 0
                    dt_5 = Time.dt_max;
            end

            if j > 3

                C_BDF2 = (gamma + 1) / (6*gamma);

                % u_ttt = SSAV_helpers.compute_deriv3(u_ttt_pts, ...
                %     dt_new, dt, dt_vals(j-2));
                % 
                % w_ttt = SSAV_helpers.compute_deriv3(w_ttt_pts, ...
                %     dt_new, dt, dt_vals(j-2));

                % if j == 4
                %     u_ttt_vals{1} = u_ttt;
                %     w_ttt_vals(1) = w_ttt;
                % 
                % elseif j == 5
                %     u_ttt_vals{2} = u_ttt;
                %     w_ttt_vals(2) = w_ttt;
                % 
                % else
                %     if j == 6
                %         u_ttt_vals{3} = u_ttt;
                %         w_ttt_vals(3) = w_ttt;
                % 
                %     else
                %         u_ttt_vals{1} = [];
                %         u_ttt_vals = u_ttt_vals(~cellfun('isempty', u_ttt_vals));
                %         u_ttt_vals{3} = u_ttt;
                % 
                %         w_ttt_vals = [w_ttt_vals(2:3); w_ttt];
                %     end
                % 
                %     u_ttt = (u_ttt_vals{3} + 0.5*u_ttt_vals{2} + ...
                %         0.25*u_ttt_vals{1}) / (1 + 0.5 + 0.25);
                % 
                %     w_ttt = (w_ttt_vals(3) + 0.5*w_ttt_vals(2) + ...
                %         0.25*w_ttt_vals(1)) / (1 + 0.5 + 0.25);
                % end

                dt_4 = (6*gamma*Para.err_tol(4) / ((gamma + 1) * ...
                abs(sum(Du .* u_ttt, 'all')*prod(Grid.d)/prod(Grid.N)) + Para.delta)) ^ (1/3);

                if isreal(dt_4) == 0
                    dt_4 = Time.dt_max;
                end

                dt_6 = (Para.err_tol(6) / (2 * C_BDF2 * abs(w_new * w_ttt) ...
                    + Para.delta)) ^ (1/3);

                if isreal(dt_6) == 0
                    dt_6 = Time.dt_max;
                end

                % [dt_prop, dt_idx(j)] = min([dt_1, dt_2, dt_3, dt_4, dt_5, dt_6]);
                % dt_prop = dt_2;
                [dt_prop, dt_idx(j)] = min([dt_1, dt_2]);
                % dt_prop = dt_4;

            else
                % [dt_prop, dt_idx(j)] = min([dt_1, dt_2, dt_3, Time.dt_max, ...
                %     dt_5, Time.dt_max]);
                [dt_prop, dt_idx(j)] = min([dt_1, dt_2]);

                % dt_prop = Time.dt_min;
                % dt_prop = dt_2;
            end

            % dt_prop = dt_5;

            min_dt = dt_prop;
            % [min_dt, dt_prop_true(j)] = min([Time.dt_max, 1.1*dt, 0.9*dt_prop]);
            % [min_dt, dt_prop_true(j)] = min([Time.dt_max, 1.5*dt, Time.dt_max*dt_prop]);

            dt_new = max([Time.dt_min, min_dt, 0.6*dt]);
        end
    end
    % dt_new = 1.10*dt;
    dt_new = min(Time.dt_max, dt_new);

    % compute new time step
    % dt = dt_new;
    % dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));
    % dt_new = max(Time.dt_min, Time.dt_max / sqrt(1 + Para.sigma * abs(E_t)^2));
    % dt_new = max(Time.dt_min, Time.dt_max / ...
    %     sqrt(1 + Para.sigma*(2*abs(E_new - E_mod)/dt)^(Para.p)));

    t = t + dt;
 
    Eu(j + 1) = E_new;
    Eu_SSAV(j) = E_mod;
    Et_vals(j) = E_t;
    mass(j) = (sum(u_new, 'all')*prod(Grid.d))/(prod(Grid.L));
    % Ett_vals(j) = E_tt;
    if save_u == 1 
        uu{j+1} = u_new;
    end

    Eu_sym(j) = 0.5*(E_new + (SSAV_helpers.compute_En(u_new_fft + ...
        gamma*(u_new_fft - u_fft), w_new + gamma*(w_new - w), ...
        Para, Grid, D, model, dim)));

    E_old = E;                      E = E_new;
    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;
    u_old_fft = u_fft;              u_fft = u_new_fft;

    if t >= plt(plt_idx)

        uu{plt_idx} = u;
        u_t_vals(plt_idx) = t;
        % mass(plt_idx) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));
        plt_idx = plt_idx + 1;

    end

    t_vals(j + 1) = t;

    if (t + dt_new) > Time.tf && Time.dt_max ~= Time.dt_min

        dt_new = Time.tf - t;

        if dt_new < Time.dt_min
            break
        end

    end

    % if (t + dt_new) > Time.tf && Time.dt_max == Time.dt_min
    % 
    %     dt_new = Time.tf - t;
    % 
    % end

end

if Time.dt_max ~= Time.dt_min

    nt = sum(t_vals > 0) + 1;
    Eu = Eu(1:nt);
    Eu_SSAV = Eu_SSAV(1:nt);
    Eu_sym = Eu_sym(1:nt);
    Et_vals = Et_vals(1:nt);
    Ett_vals = Ett_vals(1:nt);
    % l_vals = l_vals(1:nt);
    t_vals = t_vals(1:nt);
    % M_vals = M_vals(1:nt);
    dt_idx = dt_idx(1:nt);
    dt_vals = dt_vals(1:nt);
    dt_prop_true = dt_prop_true(1:nt);
    w_vals = w_vals(1:nt);
    wt_vals = wt_vals(1:nt);
    wtt_vals = wtt_vals(1:nt);
    wttt_vals = wttt_vals(1:nt);

end

% if save_u == 1
%     Results.uu = uu(~cellfun('isempty', uu));
% else
%     Results.uu = u;
% end
Results.uu = uu;
% Results.uu = u;
Results.u_t_vals = u_t_vals;
Results.Eu = Eu;
Results.Eu_SSAV = Eu_SSAV;
Results.Eu_sym = Eu_sym;
Results.Et_vals = Et_vals;
Results.Ett_vals = Ett_vals;
Results.Em = Em;
% Results.mass = mass;
Results.t_vals = t_vals;
% Results.l_vals = l_vals;
Results.num_fft = num_fft;
% Results.M_vals = M_vals;
Results.dt_idx = dt_idx;
Results.dt_vals = dt_vals;
Results.dt_prop_true = dt_prop_true;
Results.w_vals = w_vals;
Results.wt_vals = wt_vals;
Results.wtt_vals = wtt_vals;
Results.wttt_vals = wttt_vals;

end