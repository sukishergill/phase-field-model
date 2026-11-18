% HexSymmetryTest.m
%
% Hexagonal spot branch continued with the full hexagonal (D6) symmetry
% enforced, so it cannot turn into ovals, dumbbells or stripes. Starts from
% a dark hex lattice (u ~ -1 spots in u ~ +1) with spacing d at mass m0 on
% one exact hex cell (width d, height d*sqrt(3), grid spacing ~h), and runs
% continuation_FCH restricted to the symmetric subspace
% (SSAV_FCH_helpers.hex_symmetry, opts.restrict), through folds and through
% the uniform state, in both directions until m = 0. E(-u) = E(u), so the
% m < 0 half of the branch is the mirror image (u -> -u, m -> -m) of the
% m > 0 half; with the symmetry enforced the full branch is a closed loop
% over -1 < m < 1 (for d = pi: |m| < 0.957).
%
% Plots E against m (solid: stable to all even perturbations of the cell,
% dashed: unstable), with the uniform state, and profiles along the branch.
%
% Settings (set any of these before running the script to override):
%   hst_m0      starting mass (0.57)
%   hst_d       lattice spacing (pi)
%   hst_h       target grid spacing (0.08)
%   hst_replot  true: skip the run and replot from HexSymmetryTest_results.mat

if ~exist('hst_m0', 'var'), hst_m0 = 0.57; end
if ~exist('hst_d', 'var'), hst_d = pi; end
if ~exist('hst_h', 'var'), hst_h = 0.08; end
if ~exist('hst_replot', 'var'), hst_replot = false; end

out_file = 'HexSymmetryTest_results.mat';

if ~hst_replot
    Para.epsilon = 0.18;
    Para.eta1 = Para.epsilon^2;
    Para.eta2 = Para.eta1;
    Para.tau = 0;

    [u0, Grid] = SSAV_FCH_helpers.hex_lattice_state(hst_m0, hst_d, 'dark', Para, hst_h);
    [S, sym_info] = SSAV_FCH_helpers.hex_symmetry(Grid);
    fprintf('hex cell %.3f x %.3f, grid %d x %d, symmetric subspace dimension %d\n', ...
        Grid.L(1), Grid.L(2), Grid.N(1), Grid.N(2), sym_info.dim);

    opts = struct('even', true, 'symmetry', S, 'restrict', true, ...
        'stop_on_uniform', false, 'ds_max', 0.05, 'max_steps', 3000, ...
        'm_range', [0 0.99], 'verbose', false);
    t_run = tic;
    opts.direction = 1;
    b_up = continuation_FCH(u0, Grid, Para, opts);
    opts.direction = -1;
    b_down = continuation_FCH(u0, Grid, Para, opts);
    run_time = toc(t_run);
    fprintf('towards larger m: %d points, stop: %s\n', numel(b_up.m), b_up.stop_reason);
    fprintf('towards smaller m: %d points, stop: %s\n', numel(b_down.m), b_down.stop_reason);

    % one branch, from the m = 0 end of the downward run through the start
    % and along the upward run to its m = 0 end
    branch = join_branches(b_down, b_up);
    fprintf('%d points in %.1f s\n', numel(branch.m), run_time);
    fprintf('folds at m = %s\n', mat2str(branch.m(branch.fold), 4));
    fprintf('crosses the uniform state at m = %s\n', mat2str(branch.m(branch.cross_uniform), 4));
    save(out_file, 'branch', 'Grid', 'Para', 'hst_m0', 'hst_d', 'hst_h', 'run_time');
else
    load(out_file, 'branch', 'Grid', 'Para', 'hst_m0', 'hst_d', 'hst_h', 'run_time');
end

%% energy diagram
mm = linspace(-1, 1, 801);
Ec = arrayfun(@(m) SSAV_FCH_helpers.steady_energy(m*ones(Grid.N(2), Grid.N(1)), Para, Grid), mm);
stable = branch.n_neg_full == 0;

figure('Position', [100 100 1100 480], 'Color', 'w');
for panel = 1:2
    subplot(1, 2, panel);
    hold on;
    plot(mm, Ec, 'k-.', 'DisplayName', '$u = m$');
    plot_by_stability(branch.m, branch.E, stable, [0.85 0.33 0.10]);
    plot(NaN, NaN, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2, 'DisplayName', 'hex (stable)');
    plot(NaN, NaN, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2, 'DisplayName', 'hex (unstable)');
    plot(branch.m(branch.i_start), branch.E(branch.i_start), 'o', 'Color', [0.85 0.33 0.10], ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 7, 'HandleVisibility', 'off');
    xlabel('$m$', 'Interpreter', 'latex');
    ylabel('$E$', 'Interpreter', 'latex');
    xlim([0 1]);
    if panel == 1
        title(sprintf('hex symmetry enforced, $d = %.2f$, $\\epsilon = %.2f$', hst_d, Para.epsilon), ...
            'Interpreter', 'latex');
        lgd = legend('Location', 'north');
        set(lgd, 'Interpreter', 'latex');
    else
        ylim([min(branch.E) - 3e-4, min(branch.E) + 9e-3]);
        ax = gca;   ax.YAxis.Exponent = 0;   ytickformat('%.3f');
        title('zoom (dot: start)', 'Interpreter', 'latex');
    end
    set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
    box on;
    hold off;
end
exportgraphics(gcf, 'HexSymmetryTest_E_vs_m.png', 'Resolution', 150);

%% profiles along the branch (2 x 2 cells each), evenly spaced in arclength
s = [0, cumsum(sqrt(diff(branch.m).^2 + diff(branch.amp).^2))];   % approximate arclength
k_show = arrayfun(@(t) find(s >= t, 1), linspace(0, s(end), 13));
k_show = unique(k_show(1:end-1), 'stable');
figure('Position', [50 50 1500 640], 'Color', 'w');
for j = 1:numel(k_show)
    k = k_show(j);
    subplot(2, 6, j);
    imagesc([0 2*Grid.L(1)], [0 2*Grid.L(2)], repmat(branch.u{k}, 2, 2), [-1 1]);
    axis image off;
    title(sprintf('m = %+.3f\n%s', branch.m(k), ternary(stable(k), 'stable', 'unstable')), 'FontSize', 10);
end
exportgraphics(gcf, 'HexSymmetryTest_profiles.png', 'Resolution', 110);

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
