function Results = FCH_BDF2_SAV_v2(Grid, Time, Para, u, plt_save)

dt = Time.dt_min;

t = Time.t0;

% define number of time steps and time steps for plotting purposes
nmax = round(Time.tf / Time.dt_min);

plt = linspace(0, Time.tf, plt_save + 1);
plt_idx = 2;

uu = cell(plt_save + 1, 1);
uu{1} = u;
u_t_vals = zeros(1, plt_save);

for i = 2:plt_save + 1
    uu{i} = zeros(size(uu{1}));
end

u_fft = fftn(u);

Eu = zeros(1, nmax); 
t_vals = zeros(1, nmax);
dt_idx = zeros(1, nmax);
dt_prop_true = zeros(1, nmax);


mass = zeros(1, 101);
mass(1) = (sum(u, 'all')*prod(Grid.d))/(prod(Grid.L));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Compute u^1 using backward Euler %%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Wq, Wq1, Wq2] = SSAV_FCH_helpers.compute_Wq(u, Para);

G = SSAV_FCH_helpers.compute_G_2(u, Para, Wq, Wq1, Grid);

w_old = SSAV_FCH_helpers.compute_w(G, Para.B);

Eu(1) = SSAV_FCH_helpers.compute_E(u_fft, w_old, Para, Grid);

H = SSAV_FCH_helpers.compute_H(u, w_old, Para, Wq1, Wq2, Grid);
H_old = H;

r = w_old - 0.5 * sum(H.*u ,'all')*prod(Grid.d);

H2 = fftn(H);


% r_hat = u/dt + (Para.epsilon^2*(2 + Para.eta1 + 2*Para.tau^2/3) + Para.S) * ...
%     real(ifftn(Grid.k4 .* v)) + r * real(ifftn(-Grid.k2 .* H));

% r_hat = u/dt + r * real(ifftn(-Grid.k2 .* H)) + ...
%     (Para.epsilon^2*(2 + Para.eta1 + 2*Para.tau^2/3) + Para.S + ...
%     2*Para.epsilon^2*Para.alpha) * real(ifftn(Grid.k4 .* v))  ...
%     - Para.alpha^2*real(ifftn(-Grid.k2 .* v));

r_hat = u/dt + r * ifftn(-Grid.k2 .* H2, 'symmetric') + ...
    (Para.epsilon^2*(Para.eta1 + 2*Para.zeta) + Para.S + ...
    2*Para.epsilon^2*Para.alpha) * ifftn(Grid.k4 .* u_fft, 'symmetric')  ...
    - Para.alpha^2*ifftn(-Grid.k2 .* u_fft, 'symmetric');



% define linear operator P for BDF1
% P = 1/dt + Para.epsilon^4 * Grid.k6 + Para.S * Grid.k4;

P1 = Para.epsilon^4 * Grid.k6 + ...
    (Para.S + 2*Para.epsilon^2*Para.alpha) * Grid.k4 + ...
    (Para.alpha^2)*Grid.k2;

P = 1/dt + P1;

r_hat = fftn(r_hat);
psi_r = P .\ r_hat;         psi_r = ifftn(psi_r, 'symmetric');

psi_H = P .\ (G .* H2);     psi_H = ifftn(psi_H, 'symmetric');

innprod_Hu = SSAV_FCH_helpers.compute_ip(H, psi_r, psi_H, prod(Grid.d));

w = 0.5*innprod_Hu + r;

u_old = u;

% if Time.adap == 1
% 
%     [F_tilde_old, num_fft] = SSAV_helpers.compute_F_tilde(u_old, ...
%         Wq, Grid, Para, num_fft);
% 
%     [~, f] = SSAV_helpers.compute_F(u, Para.beta);
%     [F_tilde, num_fft] = SSAV_helpers.compute_F_tilde(u, u_fft, f, Grid,...
%         Para, num_fft);
% 
% end


u = 0.5 * innprod_Hu * psi_H + psi_r;
u_fft = fftn(u);

% % define linear operator P for BDF2
% P = P + 1/(2*dt);

t = t + dt;

Eu(2) = SSAV_FCH_helpers.compute_E(u_fft, w, Para, Grid);

t_vals(2) = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;
dt_new = dt;
% dt_vals(1) = dt_new;


while t < Time.tf

    j = j + 1;

    u_snew = SSAV_FCH_helpers.compute_u_snew(u, u_old);

    [Wq, Wq1, Wq2] = SSAV_FCH_helpers.compute_Wq(u_snew, Para);

    G = SSAV_FCH_helpers.compute_G_2(u, Para, Wq, Wq1, Grid);

    w_snew = SSAV_FCH_helpers.compute_w(G, Para.B);

    H_snew = SSAV_FCH_helpers.compute_H(u_snew, w_snew, Para, Wq1, Wq2,...
        Grid);

    [a, b, c] = SSAV_FCH_helpers.compute_dt_coeffs(dt, dt_new);

    H = fftn(H_snew);
    u_snew = fftn(u_snew);

    r = -(b*w + c*w_old)/a + 0.5 * sum(H_snew .* (b*u + c*u_old)/a, 'all') ...
        * prod(Grid.d);

    % r_hat = (4*u - u_old)/(2*dt) + (Para.epsilon^2*(2 + Para.eta1 + ...
    %     2*Para.tau^2/3) + Para.S) * real(ifftn(Grid.k4 .* u_snew)) + ...
    %     r * real(ifftn(-Grid.k2 .* H));

    % r_hat = (4*u - u_old)/(2*dt) + r * real(ifftn(-Grid.k2 .* H)) + ...
    %     (Para.epsilon^2*(2 + Para.eta1 + 2*Para.tau^2/3) + Para.S + ...
    %     2*Para.epsilon^2*Para.alpha) * real(ifftn(Grid.k4 .* u_snew))  ...
    %     - Para.alpha^2*real(ifftn(-Grid.k2 .* u_snew));

    r_hat = -(b*u + c*u_old) + r * ifftn(-Grid.k2 .* H, 'symmetric') + ...
        (Para.epsilon^2*(Para.eta1 + 2*Para.zeta) + Para.S + ...
        2*Para.epsilon^2*Para.alpha) * ifftn(Grid.k4 .* u_snew, 'symmetric')  ...
        - Para.alpha^2*ifftn(-Grid.k2 .* u_snew, 'symmetric');


    r_hat = fftn(r_hat);

    P = a + P1;

    psi_r = P .\ r_hat;         psi_r = ifftn(psi_r, 'symmetric');

    psi_H = P .\ (G .* H);     psi_H = ifftn(psi_H, 'symmetric');

    innprod_Hu = SSAV_FCH_helpers.compute_ip(H_snew, psi_r, psi_H, ...
        prod(Grid.d));

    w_new = 0.5*innprod_Hu + r;

    u_new = 0.5*innprod_Hu * psi_H + psi_r;

    % [u_new, w_new] = SSAV_FCH_helpers.compute_unew(u, u_old, w, w_old,...
    %     dt, dt_new, Para, Grid, P1);
    u_new_fft = fftn(u_new);

    E_new = SSAV_FCH_helpers.compute_E(u_new_fft, w_new, Para, Grid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Time.dt_max == Time.dt_min

        dt_new = Time.dt_min;

    % elseif Time.adap == 1
    %     [A_dp, rel_err, F_tilde_new, num_fft] = SSAV_helpers.compute_Adp(u_new, ...
    %         u, u_new_fft, F_tilde_old, F_tilde, dt_new, gamma, Para, ...
    %         Grid, num_fft);
    % 
    %     l = 1;
    % 
    %     while (rel_err > Para.err_tol_AM3) && (dt_new ~= Time.dt_min) && (l < 10)
    % 
    %         dt_new = max(Time.dt_min, min(A_dp, Time.dt_max));
    % 
    %         [u_new, w_new, num_fft] = SSAV_helpers.compute_unew(u, u_fft, ...
    %             u_old, u_old_fft, w, w_old, dt, dt_new, Para, ...
    %             Grid, G, P1, num_fft);
    % 
    %         u_new_fft = fftn(u_new);        num_fft = num_fft + 1;
    % 
    %         [A_dp, rel_err, F_tilde_new, num_fft] = ...
    %             SSAV_helpers.compute_Adp(u_new, u, u_new_fft, F_tilde_old, ...
    %             F_tilde, dt, dt_new, Para, Grid, num_fft);
    % 
    %         l = l + 1;
    % 
    %         if l == 10
    %             dt_new = 0.5*dt_new;
    %         end
    %
    % F_tilde_old = F_tilde;      F_tilde = F_tilde_new;
    else

        dt = dt_new;

        E_t = (E_new - Eu(j)) / dt_new;

        dt_1 = Para.err_tol(1) / (abs(E_t) + Para.delta);

        dt_2 = Para.err_tol(2) / sqrt(abs(E_t)+ Para.delta);

        [dt_prop, dt_idx(j)] = min([dt_1, dt_2]);
        % dt_prop = dt_2;

        [min_dt, dt_prop_true(j)] = min([Time.dt_max, 2*dt, dt_prop]);
        
        dt_new = max([Time.dt_min, min_dt, 0.1*dt]);

    end

    % dt_new = min(Time.dt_max, dt_new);
    % 
    % 
    % dt_new = dt;
    t = t + dt;

    Eu(j + 1) = E_new;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;


    if t >= plt(plt_idx)

        uu{plt_idx} = u;
        u_t_vals(plt_idx) = t;
        plt_idx = plt_idx + 1;

    end

    t_vals(j + 1) = t;

    if (t + dt_new) > Time.tf && Time.dt_max ~= Time.dt_min

        dt_new = Time.tf - t;

        if dt_new < Time.dt_min
            break
        end

    end

end


% Eu = Eu / prod(Grid.L);

if Time.dt_max ~= Time.dt_min

    nt = sum(t_vals > 0) + 1;
    Eu = Eu(1:nt);
    t_vals = t_vals(1:nt);
    dt_idx = zeros(1, nt);
    dt_prop_true = zeros(1, nt);

end

Results.t_vals = t_vals;
Results.Eu = Eu;
Results.uu = uu;
Results.dt_idx = dt_idx;
Results.dt_prop_true = dt_prop_true;

end
