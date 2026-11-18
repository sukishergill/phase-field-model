% SymmetricBranchesTest.m
%
% Stripe and hexagonal-spot branches with their symmetry enforced, so that
% neither can turn into the other or into anything else (ovals, dumbbells,
% wavy stripes):
%   hex spots  dark hex lattice (u ~ -1 spots in u ~ +1), spacing d, on one
%              exact hex cell (width d, height d*sqrt(3), grid spacing
%              ~h), continued by continuation_FCH restricted to the full
%              hexagonal (D6) symmetry (SSAV_FCH_helpers.hex_symmetry,
%              opts.restrict)
%   stripes    one stripe of u ~ -1 per period d, as a 1D problem (u
%              independent of y), on the same x grid as the hex cell
% Both start at mass m0 and are continued in both directions, through folds
% and through the uniform state, until m = 0. E(-u) = E(u), so the m < 0
% half of each branch is the mirror image (u -> -u, m -> -m) of the m > 0
% half; each full branch is a closed loop over -1 < m < 1. Energies are per
% unit area, so they compare directly.
%
% Stability (solid lines) is against all even perturbations of the hex
% cell, for both families: for the stripes it is recomputed on the 2D cell,
% so it includes perturbations that bend the stripes or break them into
% spots.
%
% Settings (set any of these before running the script to override):
%   sbt_m0      starting mass (0.57)
%   sbt_d       spot spacing = stripe period (pi)
%   sbt_h       target grid spacing (0.08)
%   sbt_replot  true: skip the runs and replot from SymmetricBranchesTest_results.mat

if ~exist('sbt_m0', 'var'), sbt_m0 = 0.57; end
if ~exist('sbt_d', 'var'), sbt_d = pi; end
if ~exist('sbt_h', 'var'), sbt_h = 0.08; end
if ~exist('sbt_replot', 'var'), sbt_replot = false; end

out_file = 'SymmetricBranchesTest_results.mat';

if ~sbt_replot
    Para.epsilon = 0.18;
    Para.eta1 = Para.epsilon^2;
    Para.eta2 = Para.eta1;
    Para.tau = 0;
    copts = struct('even', true, 'stop_on_uniform', false, 'ds_max', 0.05, ...
        'max_steps', 3000, 'm_range', [0 0.99], 'verbose', false);

    % ---- hex spots, D6 symmetry enforced ----
    [u0, Grid] = SSAV_FCH_helpers.hex_lattice_state(sbt_m0, sbt_d, 'dark', Para, sbt_h);
    [S, sym_info] = SSAV_FCH_helpers.hex_symmetry(Grid);
    fprintf('hex cell %.3f x %.3f, grid %d x %d, symmetric subspace dimension %d\n', ...
        Grid.L(1), Grid.L(2), Grid.N(1), Grid.N(2), sym_info.dim);
    o = copts;  o.symmetry = S;  o.restrict = true;
    t_run = tic;
    o.direction = 1;   b_up = continuation_FCH(u0, Grid, Para, o);
    o.direction = -1;  b_down = continuation_FCH(u0, Grid, Para, o);
    hex = join_branches(b_down, b_up);
    hex.stable = hex.n_neg_full == 0;
    hex.run_time = toc(t_run);
    report('hex spots', hex);

    % ---- stripes: 1D, period d, same x grid ----
    G1 = SSAV_FCH_helpers.generate_Grid(Grid.L(1), Grid.N(1), 1);
    layer = @(s) (1 - tanh(s / (sqrt(2)*Para.epsilon))) / 2;
    w = (1 - sbt_m0) / 2 * Grid.L(1);                % stripe width for area fraction (1 - m0)/2
    u1 = 1 - 2*layer(abs(G1.x) - w/2);
    u1 = u1 - mean(u1) + sbt_m0;
    o = copts;
    t_run = tic;
    o.direction = 1;   s_up = continuation_FCH(u1, G1, Para, o);
    o.direction = -1;  s_down = continuation_FCH(u1, G1, Para, o);
    stripe = join_branches(s_down, s_up);
    % stability on the 2D hex cell (all even perturbations)
    Z = SSAV_FCH_helpers.even_basis([Grid.N(2), Grid.N(1)]);
    stripe.n_neg_cell = zeros(size(stripe.m));
    for k = 1:numel(stripe.m)
        stripe.n_neg_cell(k) = n_neg_even(repmat(stripe.u{k}(:)', Grid.N(2), 1), Z, Grid, Para);
    end
    stripe.stable = stripe.n_neg_cell == 0;
    stripe.run_time = toc(t_run);
    report('stripes', stripe);

    save(out_file, 'hex', 'stripe', 'Grid', 'G1', 'Para', 'sbt_m0', 'sbt_d', 'sbt_h');
else
    load(out_file, 'hex', 'stripe', 'Grid', 'G1', 'Para', 'sbt_m0', 'sbt_d', 'sbt_h');
end

%% energy diagram (m >= 0; the m < 0 half is the mirror image)
mm = linspace(0, 1, 401);
Ec = arrayfun(@(m) SSAV_FCH_helpers.steady_energy(m*ones(Grid.N(2), Grid.N(1)), Para, Grid), mm);
fam = {stripe, hex};
names = {'stripes', 'hex spots'};
colors = [0.00 0.45 0.74; 0.85 0.33 0.10];

figure('Position', [100 100 1100 480], 'Color', 'w');
for panel = 1:2
    subplot(1, 2, panel);
    hold on;
    plot(mm, Ec, 'k-.', 'DisplayName', '$u = m$');
    for q = 1:2
        b = fam{q};
        plot_by_stability(b.m, b.E, b.stable, colors(q, :));
        plot(NaN, NaN, '-', 'Color', colors(q, :), 'LineWidth', 2, 'DisplayName', names{q});
        plot(b.m(b.i_start), b.E(b.i_start), 'o', 'Color', colors(q, :), ...
            'MarkerFaceColor', colors(q, :), 'MarkerSize', 7, 'HandleVisibility', 'off');
    end
    xlabel('$m$', 'Interpreter', 'latex');
    ylabel('$E$', 'Interpreter', 'latex');
    xlim([0 1]);
    if panel == 1
        title(sprintf('symmetry enforced, $d = %.2f$, $\\epsilon = %.2f$', sbt_d, Para.epsilon), ...
            'Interpreter', 'latex');
        lgd = legend('Location', 'north');
        set(lgd, 'Interpreter', 'latex');
    else
        E_lo = min([hex.E, stripe.E]);
        ylim([E_lo - 3e-4, E_lo + 9e-3]);
        ax = gca;   ax.YAxis.Exponent = 0;   ytickformat('%.3f');
        title('zoom (solid: stable, dashed: unstable)', 'Interpreter', 'latex');
    end
    set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
    box on;
    hold off;
end
exportgraphics(gcf, 'SymmetricBranchesTest_E_vs_m.png', 'Resolution', 150);

%% profiles along each branch (2 x 2 hex cells each), evenly spaced in arclength
figure('Position', [50 50 1500 640], 'Color', 'w');
for q = 1:2
    b = fam{q};
    s = [0, cumsum(sqrt(diff(b.m).^2 + diff(b.amp).^2))];      % approximate arclength
    k_show = arrayfun(@(t) find(s >= t - 1e-12, 1), linspace(0, s(end), 8));
    for j = 1:numel(k_show)
        k = k_show(j);
        u = b.u{k};
        if isvector(u), u = repmat(u(:)', Grid.N(2), 1); end
        subplot(2, 8, (q - 1)*8 + j);
        imagesc([0 2*Grid.L(1)], [0 2*Grid.L(2)], repmat(u, 2, 2), [-1 1]);
        axis image off;
        title(sprintf('%s\nm = %+.3f, %s', names{q}, b.m(k), ternary(b.stable(k), 'stable', 'unstable')), ...
            'FontSize', 9);
    end
end
exportgraphics(gcf, 'SymmetricBranchesTest_profiles.png', 'Resolution', 110);

function b = join_branches(b_down, b_up)
% reversed downward branch followed by the upward one (start point once)
f = {'m', 'E', 'dE', 'amp', 'n_neg', 'n_neg_full', 'fold', 'cross_uniform', 'u'};
for i = 1:numel(f)
    a = b_down.(f{i});
    c = b_up.(f{i});
    b.(f{i}) = [fliplr(a), c(2:end)];
end
b.i_start = numel(b_down.m);
b.stop_reason = {b_down.stop_reason, b_up.stop_reason};
end

function n_neg = n_neg_even(u, Z, Grid, Para)
% number of negative eigenvalues of the Hessian on the even, mean-zero
% perturbations of the cell (columns of Z): the unstable directions
sz = size(u);
n = numel(u);
ne = size(Z, 2);
JZ = zeros(n, ne);
for c0 = 1:64:ne
    cols = c0:min(c0 + 63, ne);
    JZ(:, cols) = reshape(SSAV_FCH_helpers.steady_Jv(u, reshape(Z(:, cols), [sz, numel(cols)]), ...
        Para, Grid), n, []);
end
Je = Z' * JZ;
n_neg = sum(eig((Je + Je') / 2) < -1e-6);
end

function report(name, b)
st = b.stable;
d = diff([0 st 0]);
a = find(d == 1);
z = find(d == -1) - 1;
win = arrayfun(@(i) sprintf('[%.3f, %.3f] ', min(b.m([a(i) z(i)])), max(b.m([a(i) z(i)]))), ...
    1:numel(a), 'UniformOutput', false);
fprintf('%s: %d points in %.1f s, m in [%.3f, %.3f], stop: %s | %s\n', name, numel(b.m), b.run_time, ...
    min(b.m), max(b.m), b.stop_reason{1}, b.stop_reason{2});
fprintf('   folds at m = %s, uniform crossings at m = %s\n', mat2str(b.m(b.fold), 4), ...
    mat2str(b.m(b.cross_uniform), 4));
fprintf('   stable for m in %s\n', [win{:}]);
end

function plot_by_stability(m, E, stable, color)
% solid where stable, dashed where not; segments overlap by one point
n = numel(m);
i = 1;
while i <= n
    j = i;
    while j < n && stable(j+1) == stable(i)
        j = j + 1;
    end
    k = min(j + 1, n);
    if stable(i), style = '-'; else, style = '--'; end
    plot(m(i:k), E(i:k), style, 'Color', color, 'LineWidth', 2, 'HandleVisibility', 'off');
    i = j + 1;
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
