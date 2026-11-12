% SpotStripeTest.m
%
% 2D analogue of CosTest.m: starts the pseudo-arclength continuation
% (arclength_FCH_func_2D) from two initial profiles at the same mass m0,
%   stripe  one stripe of u ~ -1 per cell width d, in a u ~ +1 background
%   spots   hexagonally packed spots of u ~ -1, nearest-neighbour spacing d
% and plots both m-E branches together (solid = stable, dashed = unstable),
% with the energy of the uniform state u = m for reference.
%
% Both run on the same periodic cell, d wide and d*sqrt(3) high: the
% smallest rectangle that tiles a hexagonal lattice exactly (two spots per
% cell, one at the centre and one shared between the four corners). The
% stripe then has period d. arclength_FCH_func_2D needs the same number of
% grid points in each direction, so the grid spacing differs in x and y
% (dy = sqrt(3)*dx). Energies are per unit area.
%
% Both profiles are even about the cell centre, as arclength_FCH_func_2D
% requires, and are sized so the u ~ -1 phase covers the area fraction
% (1 - m0)/2 that a sharp-interface state of mass m0 would have. The
% function relaxes them to steady states itself before continuing. It does
% not force the two spots of the cell to stay identical.
%
% Settings (set any of these before running the script to override):
%   sst_m0      starting mass (0.57)
%   sst_d       spot spacing = stripe period = cell width (pi)
%   sst_N       grid points per direction (56: dx = 0.056, dy = 0.097 at d = pi)
%   sst_replot  true: skip the runs and replot from SpotStripeTest_results.mat

if ~exist('sst_m0', 'var'), sst_m0 = 0.57; end
if ~exist('sst_d', 'var'), sst_d = pi; end
if ~exist('sst_N', 'var'), sst_N = 56; end
if ~exist('sst_replot', 'var'), sst_replot = false; end

out_file = 'SpotStripeTest_results.mat';

if ~sst_replot
    m0 = sst_m0;
    N = sst_N;
    d = sst_d;
    L = [d, d*sqrt(3)];
    eps = 0.18;     eta = eps^2;
    Grid = SSAV_FCH_helpers.generate_Grid(L, [N, N], 2);

    layer = @(s) (1 - tanh(s / (sqrt(2)*eps))) / 2;    % ~1 inside (s < 0), ~0 outside
    frac = (1 - m0) / 2;                    % area fraction of the u ~ -1 phase
    w = frac * L(1);                        % stripe width
    R = sqrt(frac * prod(L) / (2*pi));      % spot radius (two spots per cell)

    % spots at the cell centre and at the corners (+-d/2, +-Ly/2)
    r1 = hypot(Grid.xx, Grid.yy);
    ddx = Grid.xx - L(1)/2;   ddx = ddx - L(1)*round(ddx/L(1));
    ddy = Grid.yy - L(2)/2;   ddy = ddy - L(2)*round(ddy/L(2));
    r2 = hypot(ddx, ddy);

    names = {'stripe', 'hex spots'};
    u0s = {1 - 2*layer(abs(Grid.xx) - w/2), ...
           1 - 2*layer(r1 - R) - 2*layer(r2 - R)};

    branches = cell(numel(names), 1);
    for q = 1:numel(names)
        fprintf('\n\n===================== %s =====================\n', names{q});
        u0 = u0s{q} - mean(u0s{q}(:)) + m0;     % force exact mass = m0
        t_run = tic;
        [u_vals, E_vals, m_vals, diagnostics] = arclength_FCH_func_2D(u0, m0, eps, eta, Grid);
        branches{q} = struct('name', names{q}, 'u0', u0, 'u_vals', {u_vals}, ...
            'E_vals', E_vals, 'm_vals', m_vals, 'diagnostics', diagnostics, ...
            'i_start', diagnostics.dec.steps + 1, 'run_time', toc(t_run));
    end

    % energy of the uniform state u = m, from the same energy routine
    mm = linspace(0, 1, 501);
    Ec = zeros(size(mm));
    for i = 1:numel(mm)
        uc = mm(i) * ones(N, N);
        Ec(i) = SSAV_FCH_helpers.compute_E_direct(uc, fftn(uc), eps, eta, Grid);
    end

    save(out_file, 'branches', 'names', 'm0', 'N', 'L', 'd', 'eps', 'eta', 'mm', 'Ec', 'Grid');
else
    load(out_file, 'branches', 'names', 'm0', 'N', 'L', 'd', 'eps', 'eta', 'mm', 'Ec', 'Grid');
end

%% energy diagram: E against m, full range and zoomed on the two branches
% Thick lines: the part of each run that is still the starting pattern
% (see main_piece below). Thin, pale lines: the rest of the same run, after
% it has passed a fold at low m or gone through the uniform state.
colors = [0.00 0.45 0.74; 0.85 0.33 0.10];
figure('Position', [100 100 1100 480], 'Color', 'w');
for panel = 1:2
    subplot(1, 2, panel);
    hold on;
    for q = 1:numel(branches)
        b = branches{q};
        st = b.diagnostics.stability;
        [ia, ib] = main_piece(b);
        pale = 0.55 + 0.45*colors(q, :);
        plot(NaN, NaN, '-', 'Color', colors(q, :), 'LineWidth', 2.5, 'DisplayName', b.name);   % legend entry
        plot_branch_by_stability(b.m_vals(1:ia), b.E_vals(1:ia), st(1:ia), pale, 1);
        plot_branch_by_stability(b.m_vals(ib:end), b.E_vals(ib:end), st(ib:end), pale, 1);
        plot_branch_by_stability(b.m_vals(ia:ib), b.E_vals(ia:ib), st(ia:ib), colors(q, :), 2.5);
        plot(b.m_vals(b.i_start), b.E_vals(b.i_start), 'o', 'Color', colors(q, :), ...
            'MarkerFaceColor', colors(q, :), 'MarkerSize', 7, 'HandleVisibility', 'off');
    end
    plot(mm, Ec, 'k-.', 'DisplayName', '$u = m$');
    if panel == 1
        title(sprintf('$\\epsilon = %.2f$, cell $%.2f \\times %.2f$, $%d^2$', eps, L(1), L(2), N), 'Interpreter', 'latex');
    else
        E_all = cell2mat(cellfun(@(b) b.E_vals(:)', branches(:)', 'UniformOutput', false));
        E_lo = min(E_all);
        ylim([E_lo - 3e-4, E_lo + 9e-3]);
        ax = gca;   ax.YAxis.Exponent = 0;   ytickformat('%.3f');
        title('zoom (solid: stable, dashed: unstable)', 'Interpreter', 'latex');
    end
    xlabel('$m$', 'Interpreter', 'latex');
    ylabel('$E$', 'Interpreter', 'latex');
    xlim([0 1]);
    locs = {'northwest', 'northwest'};
    lgd = legend('Location', locs{panel});
    set(lgd, 'Interpreter', 'latex');
    set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
    box on;
    hold off;
end
sgtitle(sprintf('Stripe and hex-spot runs started at $m_0 = %.2f$ (dots), spacing $d = %.2f$; pale: same run on other patterns', m0, d), ...
    'Interpreter', 'latex', 'FontSize', 16);
exportgraphics(gcf, 'SpotStripeTest_E_vs_m.png', 'Resolution', 150);

%% the two relaxed starting states (2 x 2 cells, true aspect ratio)
figure('Position', [100 100 800 560], 'Color', 'w');
for q = 1:numel(branches)
    b = branches{q};
    subplot(1, 2, q);
    imagesc([0 2*L(1)], [0 2*L(2)], repmat(b.u_vals{b.i_start}, 2, 2), [-1 1]);
    axis image off;
    colorbar;
    title(sprintf('%s, $m = %.2f$, $E = %.5f$', b.name, b.m_vals(b.i_start), b.E_vals(b.i_start)), ...
        'Interpreter', 'latex');
end
exportgraphics(gcf, 'SpotStripeTest_start_states.png', 'Resolution', 150);

function [ia, ib] = main_piece(b)
% Index range of the run that is still the starting pattern: from the
% starting point out to the nearest fold (turning point of m) on each
% side, or the end of the run. If the run reaches the uniform state
% (std(u) < 0.01) on the stretch just beyond that fold, the piece is
% extended to there: that stretch is the same pattern shrinking to zero
% amplitude. Beyond these limits the run is on a different pattern (or on
% a shifted copy of the same one).

n = numel(b.m_vals);
folds = find(diff(sign(diff(b.m_vals))) ~= 0) + 1;
is_uniform = cellfun(@(u) std(u(:)) < 0.01, b.u_vals);

lo = [1, folds(folds < b.i_start)];
hi = [folds(folds > b.i_start), n];
ia = lo(end);
ib = hi(1);

if numel(hi) > 1
    k = find(is_uniform(ib:hi(2)), 1, 'first');
    if ~isempty(k), ib = ib + k - 1; end
end
if numel(lo) > 1
    k = find(is_uniform(lo(end-1):ia), 1, 'last');
    if ~isempty(k), ia = lo(end-1) + k - 1; end
end
end

function plot_branch_by_stability(m_vals, E_vals, stability, color, width)
% Plot one branch, solid where stable (stability == +1) and dashed where
% unstable (-1). Consecutive same-sign runs are separate segments in one
% colour, each extended one point into its neighbour (as in CosTest.m).
% No legend entries: the caller adds one per branch.

n = numel(m_vals);
i = 1;
while i <= n
    j = i;
    while j < n && stability(j+1) == stability(i)
        j = j + 1;
    end
    idx_end = min(j+1, n);

    if stability(i) > 0
        style = '-';
    else
        style = '--';
    end

    plot(m_vals(i:idx_end), E_vals(i:idx_end), style, 'Color', color, ...
        'LineWidth', width, 'HandleVisibility', 'off');

    i = j + 1;
end
end
