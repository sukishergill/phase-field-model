function Results = FCH_BDF2_SAV(Grid, Time, Para, u, plt_save)

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

[F, F1, F2] = SSAV_FCH_helpers.compute_F(u, Para.tau);

G = SSAV_FCH_helpers.compute_G(u, Para, F, F1, Grid);

w_old = SSAV_FCH_helpers.compute_w(G, Para.B);

Eu(1) = SSAV_FCH_helpers.compute_E(u_fft, w_old, Para, Grid);

H = SSAV_FCH_helpers.compute_H(u, w_old, Para, F1, F2, Grid);
H_old = H;

r = w_old - 0.5 * sum(H.*u, 'all')*prod(Grid.d);
H2 = fftn(H);

r_hat = u/dt + (Para.epsilon^2*(2 + Para.eta1 + 2*Para.tau^2/3) + Para.S) * ...
    ifftn(Grid.k4 .* u_fft, 'symmetric') + r * ifftn(-Grid.k2 .* H2, 'symmetric');

% define linear operator P for BDF1
P1 = Para.epsilon^4 * Grid.k6 + Para.S * Grid.k4;
P = 1/dt + P1;

r_hat = fftn(r_hat);
psi_r = P .\ r_hat;         psi_r = ifftn(psi_r, 'symmetric');

psi_H = P .\ (G .* H2);     psi_H = ifftn(psi_H, 'symmetric');

innprod_Hu = SSAV_FCH_helpers.compute_ip(H, psi_r, psi_H, prod(Grid.d));

% w = innprod_Hu + r;
w = 0.5*innprod_Hu + r;

u_old = u;
u_old_fft = u_fft;

u = 0.5 * innprod_Hu * psi_H + psi_r;
u_fft = fftn(u);

t = t + dt;

Eu(2) = SSAV_FCH_helpers.compute_E(u_fft, w, Para, Grid);

t_vals(2) = t;
% define linear operator P for BDF2
% P = P + 1/(2*dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main time loop %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;
dt_new = dt;

while t < Time.tf

    gamma = dt_new / dt;

    j = j + 1;

    u_snew = SSAV_FCH_helpers.compute_u_snew(u, u_old);

    [F, F1, F2] = SSAV_FCH_helpers.compute_F(u_snew, Para.tau);

    G = SSAV_FCH_helpers.compute_G(u, Para, F, F1, Grid);

    w_snew = SSAV_FCH_helpers.compute_w(G, Para.B);

    H_snew = SSAV_FCH_helpers.compute_H(u_snew, w_snew, Para, F1, F2,...
        Grid);

    [a, b, c] = SSAV_FCH_helpers.compute_dt_coeffs(dt, dt_new);

    H = fftn(H_snew);
    u_snew = fftn(u_snew);

    r = -(b*w + c*w_old)/a + 0.5 * sum(H_snew .* (b*u + c*u_old)/a, 'all') ...
        * prod(Grid.d);

    r_hat = -(b*u + c*u_old) + (Para.epsilon^2*(2 + Para.eta1 + ...
        2*Para.tau^2/3) + Para.S) * ifftn(Grid.k4 .* u_snew, 'symmetric') + ...
        r * ifftn(-Grid.k2 .* H, 'symmetric');

    r_hat = fftn(r_hat);

    P = a + P1;

    psi_r = P .\ r_hat;         psi_r = ifftn(psi_r, 'symmetric');

    psi_H = P .\ (G .* H);     psi_H = ifftn(psi_H, 'symmetric');

    innprod_Hu = SSAV_FCH_helpers.compute_ip(H_snew, psi_r, psi_H, ...
        prod(Grid.d));

    w_new = 0.5*innprod_Hu + r;

    u_new = 0.5*innprod_Hu * psi_H + psi_r;
    u_new_fft = fftn(u_new);

    E_new = SSAV_FCH_helpers.compute_E(u_new_fft, w_new, Para, Grid);

    E_t = (E_new - Eu(j)) / dt_new;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Adaptive time stepping %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Time.dt_max == Time.dt_min
        dt_new = Time.dt_min;

    else

        dt = dt_new;

        dt_1 = Para.err_tol(1) / (abs(E_t) + Para.delta);

        dt_2 = Para.err_tol(2) / sqrt(abs(E_t)+ Para.delta);

        [dt_prop, dt_idx(j)] = min([dt_1, dt_2]);

        [min_dt, dt_prop_true(j)] = min([Time.dt_max, 2*dt, dt_prop]);

        dt_new = max([Time.dt_min, min_dt, 0.1*dt]);

    end

    t = t + dt;
    Eu(j + 1) = E_new;

    w_old = w;                      w = w_new;
    u_old = u;                      u = u_new;
    u_old_fft = u_fft;              u_fft = u_new_fft;

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