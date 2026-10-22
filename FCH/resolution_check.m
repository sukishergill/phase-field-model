% Grid-resolution check of an energy comparison at fixed mass: the sweep
% ("other") state at m_target against the best n-stripe state, each
% Fourier-interpolated to finer grids and re-polished at the same m.
% If dE of each converges and the gap keeps its sign and size, the
% comparison is resolved.
%
% Settings (override before running): m_target (0.35), i_sweep (index of
% the sweep state, 4 = m 0.35), n_stripe (6), net_factors ([1 2]: 128^2 and
% 256^2), stripe_factors ([1 2 4 8]).

if ~exist('m_target', 'var'), m_target = 0.35; end
if ~exist('i_sweep', 'var'), i_sweep = 4; end
if ~exist('n_stripe', 'var'), n_stripe = 6; end
if ~exist('net_factors', 'var'), net_factors = [1 2]; end
if ~exist('stripe_factors', 'var'), stripe_factors = [1 2 4 8]; end

pol_files = dir('ex5_m_sweep_*_polished.mat');
[~, newest] = max([pol_files.datenum]);
in_file = pol_files(newest).name;
PL = load(in_file, 'u_polished', 'm_vals', 'Para');
SW = load(strrep(in_file, '_polished', ''), 'Grids');
Para = PL.Para;
popts = struct('even', true, 'verbose', false, 'maxit', 3000);

%% sweep ("other") state
u0 = PL.u_polished{i_sweep};
G0 = SW.Grids{i_sweep};
assert(abs(mean(u0(:)) - m_target) < 1e-12, 'sweep state %d has m = %.4f, not %.4f', i_sweep, mean(u0(:)), m_target);
net = struct('N', {}, 'dE', {}, 'res', {}, 'time', {});
for f = net_factors
    N = f * G0.N(1);
    G = SSAV_FCH_helpers.generate_Grid(G0.L, [N N], 2);
    u = SSAV_FCH_helpers.fourier_resample(u0, [N N]);
    t0 = tic;
    [u, inf_n] = newton_polish_FCH(u, G, Para, popts);
    net(end+1) = struct('N', N, 'dE', inf_n.dE, 'res', inf_n.res_hist(end), 'time', toc(t0)); %#ok<SAGROW>
    fprintf('other  (start m %.2f) %4d^2: dE = %+.10f  |R| %.1e  converged %d  (%.0f s)\n', ...
        PL.m_vals(i_sweep), N, inf_n.dE, inf_n.res_hist(end), inf_n.converged, net(end).time);
end

%% n identical stripes at the same m: nearest point on the lowest n-stripe branch
Rf = load_reference_patterns();
rows = find(Rf.n_stripes == n_stripe);
cand = Rf.stripes(rows, :);
[~, w] = lowest_dE_at_m(cand, m_target, true);
b = cand{w};
[~, k] = min(abs(b.m - m_target) + 1e3*(b.n_neg ~= 0));
us0 = b.u{k} + (m_target - b.m(k));          % shift to exactly m_target
Gs0 = b.Grid;
stripes_res = struct('N', {}, 'dE', {}, 'res', {});
for f = stripe_factors
    N = f * Gs0.N;
    G = SSAV_FCH_helpers.generate_Grid(Gs0.L, N, 1);
    u = SSAV_FCH_helpers.fourier_resample(us0, N);
    [u, inf_s] = newton_polish_FCH(u, G, Para, popts);
    stripes_res(end+1) = struct('N', N, 'dE', inf_s.dE, 'res', inf_s.res_hist(end)); %#ok<SAGROW>
    fprintf('stripes n = %d, cell %4d points (dx %.4f): dE = %+.10f  |R| %.1e  converged %d\n', ...
        n_stripe, N, Gs0.L/N, inf_s.dE, inf_s.res_hist(end), inf_s.converged);
end

%% gap at matching resolution (same dx factor)
fprintf('\ngap other - stripes at m = %.3f (negative: the other state is lower)\n', m_target);
for j = 1:numel(net_factors)
    js = find(stripe_factors == net_factors(j), 1);
    if isempty(js), continue, end
    fprintf('  resolution x%d: %+.3e\n', net_factors(j), net(j).dE - stripes_res(js).dE);
end
fprintf('  finest each:    %+.3e\n', net(end).dE - stripes_res(end).dE);
