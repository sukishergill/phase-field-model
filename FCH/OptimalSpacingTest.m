% OptimalSpacingTest.m
%
% Stripes, hex network and hex spots on cells of locally optimal size:
% continuation_FCH_cell with optimize_d = true (the cell size d is an
% unknown with dE/dd = 0), each started at one mass from the best spacing
% of SpacingScan, then continued in m in both directions (0 <= m < 1).
%   stripes      1D, period d                      start m = 0.40, d = 2.4
%   hex network  light hex lattice (dark walls)    start m = 0.25, d = 2.5
%   hex spots    dark hex lattice                  start m = 0.50, d = 4.4
% Grid point counts are fixed per family, chosen for the largest expected d.
%
% Plots E against m (solid: a local minimum in both the profile and d,
% n_neg_d = 0; dashed: not), the uniform state, the fixed-d scan envelope
% from SpacingScan_results.mat if present (thin dotted), and d against m.
% Stability for 'hex' is within the D6-symmetric subspace, and for the
% stripes within 1D; see the report printed at the end.
%
% Settings (set any of these before running the script to override):
%   ost_replot  true: skip the runs and replot from OptimalSpacingTest_results.mat

if ~exist('ost_replot', 'var'), ost_replot = false; end
out_file = 'OptimalSpacingTest_results.mat';

if ~ost_replot
    Para.epsilon = 0.18;
    Para.eta1 = Para.epsilon^2;
    Para.eta2 = Para.eta1;
    Para.tau = 0;
    layer = @(s) (1 - tanh(s / (sqrt(2)*Para.epsilon))) / 2;
    base = struct('optimize_d', true, 'verbose', false, 'max_steps', 3000, 'd_range', [0.8 12]);

    % starting states
    m0 = 0.40;  d0 = 2.4;
    G1 = SSAV_FCH_helpers.generate_Grid(d0, 128, 1);
    u = 1 - 2*layer(abs(G1.x) - (1 - m0)/2*d0/2);
    starts(1) = struct('name', 'stripes', 'cell', 'stripe', 'u', u - mean(u) + m0, 'd', d0);
    u = SSAV_FCH_helpers.hex_lattice_state(0.25, 2.5, 'light', Para, 2.5/48);
    starts(2) = struct('name', 'hex network', 'cell', 'hex', 'u', u, 'd', 2.5);
    u = SSAV_FCH_helpers.hex_lattice_state(0.50, 4.4, 'dark', Para, 4.4/72);
    starts(3) = struct('name', 'hex spots', 'cell', 'hex', 'u', u, 'd', 4.4);

    fam = cell(1, 3);
    for q = 1:3
        o = base;
        o.cell = starts(q).cell;
        t0 = tic;
        o.direction = 1;   b_up = continuation_FCH_cell(starts(q).u, starts(q).d, Para, o);
        o.direction = -1;  b_dn = continuation_FCH_cell(starts(q).u, starts(q).d, Para, o);
        f = {'m', 'd', 'E', 'amp', 'n_neg', 'n_neg_d', 'G', 'fold', 'u'};
        b = struct('name', starts(q).name);
        for i = 1:numel(f)
            b.(f{i}) = [fliplr(b_dn.(f{i})), b_up.(f{i})(2:end)];
        end
        b.stop_reason = {b_dn.stop_reason, b_up.stop_reason};
        b.run_time = toc(t0);
        fam{q} = b;
        fprintf('%-12s %4d points in %4.1f s, m in [%.3f, %.3f], d in [%.3f, %.3f], max |G| %.1e, stops: %s | %s\n', ...
            b.name, numel(b.m), b.run_time, min(b.m), max(b.m), min(b.d), max(b.d), max(abs(b.G)), b.stop_reason{:});
    end
    save(out_file, 'fam', 'Para');
else
    load(out_file, 'fam', 'Para');
end

%% lowest family at each m (points that are local minima in profile and d)
mq = 0:0.001:0.99;
Ebest = inf(3, numel(mq));
for q = 1:3
    b = fam{q};
    ok = b.n_neg_d(1:end-1) == 0 & b.n_neg_d(2:end) == 0 & b.amp(1:end-1) > 0.02 & b.amp(2:end) > 0.02;
    for k = find(ok)
        if b.m(k) == b.m(k+1), continue, end
        sel = mq >= min(b.m(k), b.m(k+1)) & mq <= max(b.m(k), b.m(k+1));
        E = b.E(k) + (b.E(k+1) - b.E(k)) * (mq(sel) - b.m(k)) / (b.m(k+1) - b.m(k));
        Ebest(q, sel) = min(Ebest(q, sel), E);
    end
end
Eu = 0.5*(mq.^3 - mq).^2 - Para.eta2*(1 - mq.^2).^2/4;
names = [cellfun(@(b) b.name, fam, 'UniformOutput', false), {'uniform'}];
[~, win] = min([Ebest; Eu], [], 1);
ch = [1, find(diff(win) ~= 0) + 1];
fprintf('\nlowest energy (local minima in profile and cell size):\n');
for i = 1:numel(ch)
    fprintf('  from m = %.3f: %s\n', mq(ch(i)), names{win(ch(i))});
end

%% plots
colors = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10];
have_scan = exist('SpacingScan_results.mat', 'file') == 2;
figure('Position', [100 100 1150 480], 'Color', 'w');
subplot(1, 2, 1);
hold on;
plot(mq, Eu, 'k-.', 'DisplayName', '$u = m$');
for q = 1:3
    b = fam{q};
    plot_by_stability(b.m, b.E, b.n_neg_d == 0, colors(q, :));
    plot(NaN, NaN, '-', 'Color', colors(q, :), 'LineWidth', 2, 'DisplayName', [b.name ', optimal $d$']);
end
if have_scan
    sc = scan_envelope('SpacingScan_results.mat');
    for q = 1:3
        plot(sc.m, sc.E(q, :), ':', 'Color', colors(q, :), 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    plot(NaN, NaN, 'k:', 'DisplayName', 'best of the $d$ scan (1.6--5)');
end
xlabel('$m$', 'Interpreter', 'latex');
ylabel('$E$', 'Interpreter', 'latex');
xlim([0 1]);
ylim([-0.0085, 0.002]);
ax = gca;   ax.YAxis.Exponent = 0;   ytickformat('%.3f');
lgd = legend('Location', 'northwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('locally optimal cell size, $\\epsilon = %.2f$ (dashed: not a minimum)', Para.epsilon), ...
    'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
box on;
hold off;

subplot(1, 2, 2);
hold on;
for q = 1:3
    b = fam{q};
    plot_by_stability(b.m, b.d, b.n_neg_d == 0, colors(q, :));
    plot(NaN, NaN, '-', 'Color', colors(q, :), 'LineWidth', 2, 'DisplayName', b.name);
end
xlabel('$m$', 'Interpreter', 'latex');
ylabel('optimal spacing $d$', 'Interpreter', 'latex');
xlim([0 1]);
lgd = legend('Location', 'northwest');
set(lgd, 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
box on;
hold off;
exportgraphics(gcf, 'OptimalSpacingTest_E_vs_m.png', 'Resolution', 150);

function sc = scan_envelope(file)
% lowest stable fixed-d energies per family from SpacingScan (same
% classification as SpacingScan.m, cell size as d)
S = load(file, 'H', 'St', 'ss_d');
sc.m = 0:0.0025:0.99;
sc.E = inf(3, numel(sc.m));
for a = 1:numel(S.ss_d)
    for q = 1:3
        if q == 1
            b = S.St{a};
            if ~b.ok, continue, end
            same = true(1, numel(b.m) - 1);
        else
            b = S.H{a};
            if ~b.ok, continue, end
            part = cellfun(@(u, m) dark_net(double(u) - m), b.u, num2cell(b.m));
            in = part == (q == 2);
            same = b.run(1:end-1) == b.run(2:end) & in(1:end-1) & in(2:end);
        end
        ok = same & b.stable(1:end-1) & b.stable(2:end) & b.amp(1:end-1) > 0.02 & b.amp(2:end) > 0.02;
        for k = find(ok)
            if b.m(k) == b.m(k+1), continue, end
            sel = sc.m >= min(b.m(k), b.m(k+1)) & sc.m <= max(b.m(k), b.m(k+1));
            E = b.E(k) + (b.E(k+1) - b.E(k)) * (sc.m(sel) - b.m(k)) / (b.m(k+1) - b.m(k));
            sc.E(q, sel) = min(sc.E(q, sel), E);
        end
    end
end
sc.E(isinf(sc.E)) = NaN;
end

function tf = dark_net(v)
% the dark region (v < 0) is connected across the cell (see SpacingScan.m)
[ny, nx] = size(v);
mask = repmat(v < 0, 3, 3);
idx = find(mask);
if isempty(idx), tf = false; return, end
[r, c] = ind2sub(size(mask), idx);
lab = zeros(size(mask));
lab(idx) = 1:numel(idx);
right = mask(:, [2:end, end]) & mask;  right(:, end) = false;
down = mask([2:end, end], :) & mask;   down(end, :) = false;
e1 = [lab(right), lab(circshift(right, 1, 2))];
e2 = [lab(down), lab(circshift(down, 1, 1))];
comp = conncomp(graph([e1(:, 1); e2(:, 1)], [e1(:, 2); e2(:, 2)], [], numel(idx)));
in_mid = r > ny & r <= 2*ny & c > nx & c <= 2*nx;
tf = false;
for cc = unique(comp(in_mid))
    sel = comp == cc;
    if max(r(sel)) - min(r(sel)) + 1 > ny || max(c(sel)) - min(c(sel)) + 1 > nx
        tf = true;
        return
    end
end
end

function plot_by_stability(x, y, stable, color)
% solid where stable, dashed where not; segments overlap by one point
n = numel(x);
i = 1;
while i <= n
    j = i;
    while j < n && stable(j+1) == stable(i)
        j = j + 1;
    end
    k = min(j + 1, n);
    if stable(i), style = '-'; else, style = '--'; end
    plot(x(i:k), y(i:k), style, 'Color', color, 'LineWidth', 2, 'HandleVisibility', 'off');
    i = j + 1;
end
end
