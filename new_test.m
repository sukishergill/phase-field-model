% new_test.m
%
% Self-convergence test for the SSAV_helpers adaptive (Time.adap == 5,
% energy-equidistribution) BDF2 controller on a GENUINE, UNFORCED
% Cahn-Hilliard gradient flow. See earlier discussion for why this
% (self-convergence on real unforced coarsening dynamics) is a more
% faithful test than comparing against a forced manufactured solution.
%
% Averages over several random-seed realizations of the initial
% condition, since a single realization's coarsening events (domain
% mergers etc.) happen at seed-dependent times and can make a one-off
% rate estimate noisy or occasionally badly behaved at the coarse end
% of the sweep.
%
% Two rates are reported:
%   rate_u  - self-convergence rate of the solution u itself (theory:
%             should be second order in tau_1).
%   rate_E  - rate of max_n |Esym^n - E(u^n)| (theory: proven first
%             order in tau_1 by construction -- Prop. 2.1). Indices 1-2
%             of Eu_sym/Eu are excluded from this max: they are
%             startup-indexing artifacts in SSAV.m (Eu_sym(2) is never
%             assigned and stays at its preallocated 0; Eu_sym(1)
%             compares against the wrong time level), not real
%             per-step discrepancies, and are O(1) regardless of tau_1.

model = 1;              % Cahn-Hilliard
dim   = 2;

Para.alpha      = 0;
Para.beta       = 1;
Para.epsilon    = 0.08;
Para.rho_s      = 0.9;
Para.sigma      = 2000;
Para.m          = .5;
Para.M          = 1;
Para.B          = 1;    % const. that ensures positive radicand
Para.S          = 2;    % positive stabilizing parameter S ~ ||f(u)||_inf
Para.p          = 3;
Para.delta      = 1e-8;
Para.forcing_on = 0;    % genuine unforced gradient flow (see SSAV.m)

L = [2*pi, 2*pi];
N = [128, 128];
Grid = SSAV_helpers.generate_Grid(L, N, dim, model);

Time.t0   = 0;
Time.tf   = 5;          % adjust if this doesn't cover real coarsening activity
Time.adap = 5;

% tau_1 sweep, plus one much finer reference run used as ground truth
% for self-convergence. tau_ref tracks tau_vals so it stays strictly
% finer than every swept value even if the sweep is retuned later.
tau_vals = [1.6e-1 8e-2 4e-2, 2e-2, 1e-2, 5e-3, 2.5e-3]/256;
% tau_vals = [1.6e-1 8e-2 4e-2, 2e-2, 1e-2, 5e-3, 2.5e-3];

tau_ref  = min(tau_vals) / 8;
tau_all  = [tau_vals, tau_ref];
n_sweep  = numel(tau_vals);
n_runs   = numel(tau_all);

seed_list = 1:5;
n_seeds   = numel(seed_list);

u_err_all    = zeros(n_seeds, n_sweep);
Ediscrep_all = zeros(n_seeds, n_sweep);
reject_all   = zeros(n_seeds, n_sweep);

% Band-limiting mask for the initial data (see below) -- depends only on
% Grid, so build it once outside the seed loop.
k_max = max(Grid.k(:));
k_cut = k_max / 4;
mask  = Grid.k <= k_cut;

for s = 1:n_seeds
    seed = seed_list(s);
    fprintf('\n########## seed %d/%d (rng=%d) ##########\n', s, n_seeds, seed);

    % Small-amplitude, BAND-LIMITED random initial data (spinodal
    % decomposition). Fully random (white) noise puts equal energy in
    % every mode, including ones right at the grid's Nyquist frequency
    % that are barely resolved -- filtering those out first gives a
    % well-resolved initial condition instead. Fourier-filter by
    % zeroing everything with |k| above k_cut = k_max/2, then
    % renormalize since discarding half the spectrum (by radius)
    % removes most of the original variance.
    rng(seed);
    noise = rand(size(Grid.xx)) - 0.5;
    noise_hat = fftn(noise) .* mask;
    noise_filt = ifftn(noise_hat, 'symmetric');
    noise_filt = noise_filt / max(abs(noise_filt(:)));      % restore unit amplitude

    u0 = Para.m + 0.1*noise_filt;
    u0 = u0 - mean(u0(:)) + Para.m;    % force exact mass = Para.m

    u_final = cell(1, n_runs);

    for i = 1:n_runs
        tau_1 = tau_all(i);
        fprintf('\n=== seed %d, run %d/%d: tau_1 = %g ===\n', seed, i, n_runs, tau_1);

        Para.err_tol = [tau_1, tau_1, 1e-5, 1e-6, 1e-9, 1e-4];

        Time.dt_min = tau_1 / 100;
        Time.dt_max = tau_1 * 100;

        Results = SSAV(Grid, Time, Para, u0, model, dim, 0, 1);

        u_final{i} = Results.uu{end};

        if i <= n_sweep
            Ediscrep_all(s, i) = max(abs(Results.Eu_sym(3:end) - Results.Eu(3:end)));
            reject_all(s, i)   = Results.reject_steps;
        end
    end

    u_ref = u_final{end};
    for i = 1:n_sweep
        u_err_all(s, i) = sqrt(sum((u_final{i} - u_ref).^2, 'all') / prod(Grid.N));
    end
end

% Average over seeds. Errors span several decades across the sweep, so
% the geometric mean (computed as exp(mean(log(.)))  -- no toolbox
% needed) is the more natural average: it keeps one seed's occasional
% large/small outlier from dominating the way an arithmetic mean would.
u_err     = exp(mean(log(u_err_all), 1));
Ediscrep  = exp(mean(log(Ediscrep_all), 1));
reject_ct = mean(reject_all, 1);

rate_u = zeros(1, n_sweep - 1);
rate_E = zeros(1, n_sweep - 1);
for i = 1:n_sweep - 1
    rate_u(i) = log(u_err(i+1) / u_err(i)) / log(tau_vals(i+1) / tau_vals(i));
    rate_E(i) = log(Ediscrep(i+1) / Ediscrep(i)) / log(tau_vals(i+1) / tau_vals(i));
end

fprintf('\n=== Self-convergence vs. reference (tau_1 = %g), averaged over %d seeds ===\n', ...
    tau_ref, n_seeds);
fprintf('%10s %14s %10s %16s %10s %10s\n', ...
    'tau_1', 'u_err', 'rate_u', 'max|Esym-E|', 'rate_E', 'rejects');
for i = 1:n_sweep
    if i < n_sweep
        fprintf('%10g %14.6e %10.3f %16.6e %10.3f %10.1f\n', ...
            tau_vals(i), u_err(i), rate_u(i), Ediscrep(i), rate_E(i), reject_ct(i));
    else
        fprintf('%10g %14.6e %10s %16.6e %10s %10.1f\n', ...
            tau_vals(i), u_err(i), '--', Ediscrep(i), '--', reject_ct(i));
    end
end

figure;
loglog(tau_vals, u_err, 'o-', 'LineWidth', 2, 'DisplayName', 'self-convergence error in u');
hold on;
loglog(tau_vals, Ediscrep, 's-', 'LineWidth', 2, 'DisplayName', 'max|Esym-E|');
% Anchor the reference slopes at the FINEST tau, not the coarsest: the
% coarse end of the sweep can show pre-asymptotic behaviour (large step
% sizes, occasional branch-selection sensitivity -- e.g. rates like
% 0, 28, then settling to ~2), so a guide line forced through the
% coarse end can be badly offset from the trend the fine end actually
% follows. Anchoring at tau_vals(end) makes the guide line agree with
% the (trustworthy) asymptotic regime by construction.
loglog(tau_vals, u_err(end) * (tau_vals/tau_vals(end)).^2, 'k--', 'DisplayName', 'O(\tau^2) reference');
loglog(tau_vals, Ediscrep(end) * (tau_vals/tau_vals(end)), 'k:', 'DisplayName', 'O(\tau) reference');
hold off;
xlabel('\tau_1');
ylabel('error');
legend('Location', 'best');
set(gca, 'FontSize', 16);
title(sprintf('Self-convergence, unforced adaptive CH coarsening run (%d-seed avg)', n_seeds));

save('new_test_results.mat', 'tau_vals', 'tau_ref', 'u_err', 'rate_u', ...
    'Ediscrep', 'rate_E', 'reject_ct', 'u_err_all', 'Ediscrep_all', ...
    'reject_all', 'seed_list', 'Para', 'Time', 'Grid');
