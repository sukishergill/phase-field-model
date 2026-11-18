% SpacingScan.m
%
% For each spacing d in ss_d, the stripe and hexagonal branches with their
% symmetry enforced (symmetric_branches: stripes of period d as a 1D
% problem; hex states on one exact hex cell restricted to D6 symmetry,
% started from dark and from light spots), 0 <= m < 1. Then, at each m,
% the lowest energy over d of the STABLE points of each family:
%   stripes, hex spots (isolated dark domains), hex network (dark region
%   connected across the cell: honeycomb walls, triangle or kagome
%   networks around light domains).
% A hex cell of size d can also hold patterns with a shorter spacing
% (d/sqrt(3), d/2, ...), so each state is also given its effective spacing
% d_eff from its dominant Fourier mode: d_eff = 4*pi/(sqrt(3)*|k|) for hex
% states, 2*pi/|k| for stripes.
% Stable = no negative Hessian eigenvalue for even perturbations of the
% hex cell (2D, so stripes may bend or break up). Energies are per unit
% area, so different d compare directly. E(-u) = E(u): m < 0 mirrors m > 0.
%
% Settings (set any of these before running the script to override):
%   ss_d       spacings (1.6:0.2:5.0)
%   ss_h       target grid spacing (0.08)
%   ss_replot  true: skip the runs and replot from SpacingScan_results.mat
% Runs the spacings in parallel when a parallel pool is available.

if ~exist('ss_d', 'var'), ss_d = 1.6:0.2:5.0; end
if ~exist('ss_h', 'var'), ss_h = 0.08; end
if ~exist('ss_replot', 'var'), ss_replot = false; end

out_file = 'SpacingScan_results.mat';

if ~ss_replot
    Para.epsilon = 0.18;
    Para.eta1 = Para.epsilon^2;
    Para.eta2 = Para.eta1;
    Para.tau = 0;
    nd = numel(ss_d);
    H = cell(1, nd);
    St = cell(1, nd);
    T = zeros(1, nd);
    bopts = struct('h', ss_h);
    t_all = tic;
    parfor a = 1:nd
        t0 = tic;
        [H{a}, St{a}] = symmetric_branches(ss_d(a), Para, bopts);
        T(a) = toc(t0);
    end
    fprintf('%d spacings in %.0f s (%.0f s of compute)\n', nd, toc(t_all), sum(T));
    save(out_file, 'H', 'St', 'T', 'ss_d', 'ss_h', 'Para', '-v7.3');
else
    load(out_file, 'H', 'St', 'T', 'ss_d', 'ss_h', 'Para');
end

%% classify the hex states, and the effective spacing of every state
for a = 1:numel(ss_d)
    L = [ss_d(a), ss_d(a)*sqrt(3)];
    if H{a}.ok
        n = numel(H{a}.m);
        H{a}.part2 = cell(1, n);
        H{a}.d_eff = zeros(1, n);
        for k = 1:n
            u = double(H{a}.u{k});
            H{a}.part2{k} = ternary(dark_percolates(u - H{a}.m(k)), 'network', 'spots');
            H{a}.d_eff(k) = 4*pi / (sqrt(3) * dominant_k(u, L));
        end
    end
    if St{a}.ok
        St{a}.d_eff = cellfun(@(u) 2*pi / dominant_k(double(u(:))', [ss_d(a), 1]), St{a}.u);
    end
end

%% lowest stable energy over d, per family
mq = 0:0.0025:0.99;
fam = {'stripes', 'hex spots', 'hex network'};
Emin = inf(numel(fam), numel(mq));
dmin = nan(numel(fam), numel(mq));
src = zeros(numel(fam), numel(mq), 2);     % (spacing index, point index) of the lower end
for a = 1:numel(ss_d)
    for q = 1:3
        if q == 1
            b = St{a};
            if ~b.ok, continue, end
            same = true(1, numel(b.m) - 1);
        else
            b = H{a};
            if ~b.ok, continue, end
            want = ternary(q == 2, 'spots', 'network');
            in = strcmp(b.part2, want);
            same = b.run(1:end-1) == b.run(2:end) & in(1:end-1) & in(2:end);
        end
        patterned = b.amp > 0.02;                  % leave out the (near-)uniform state
        ok = same & b.stable(1:end-1) & b.stable(2:end) & patterned(1:end-1) & patterned(2:end);
        for k = find(ok)
            m1 = b.m(k);  m2 = b.m(k+1);
            sel = mq >= min(m1, m2) & mq <= max(m1, m2);
            if ~any(sel) || m1 == m2, continue, end
            E = b.E(k) + (b.E(k+1) - b.E(k)) * (mq(sel) - m1) / (m2 - m1);
            idx = find(sel);
            better = E < Emin(q, idx);
            Emin(q, idx(better)) = E(better);
            dmin(q, idx(better)) = b.d_eff(k);
            src(q, idx(better), 1) = a;
            src(q, idx(better), 2) = k;
        end
    end
end
Ec = arrayfun(@(m) 0.5*(m^3 - m)^2 - Para.eta2*(1 - m^2)^2/4, mq);     % uniform state u = m
Emin(isinf(Emin)) = NaN;

%% table: lowest stable energy (and its d) per family
fprintf('\n   m     uniform   |  stripes (d_eff)    hex spots (d_eff)   hex network (d_eff)\n');
for mt = 0.05:0.05:0.95
    [~, i] = min(abs(mq - mt));
    fprintf('  %.2f  %+.5f  |', mq(i), Ec(i));
    for q = 1:3
        if isnan(Emin(q, i))
            fprintf('  %-18s', '-');
        else
            fprintf('  %+.5f (%.1f)  ', Emin(q, i), dmin(q, i));
        end
    end
    fprintf('\n');
end

%% plots: lowest stable energy over d, and the d that attains it
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
figure('Position', [100 100 1150 480], 'Color', 'w');
subplot(1, 2, 1);
hold on;
plot(mq, Ec, 'k-.', 'DisplayName', '$u = m$');
for q = 1:3
    plot(mq, Emin(q, :), '-', 'Color', colors(q, :), 'LineWidth', 2, 'DisplayName', fam{q});
end
xlabel('$m$', 'Interpreter', 'latex');
ylabel('lowest stable $E$ over $d$', 'Interpreter', 'latex');
xlim([0 1]);
E_lo = min(Emin(:), [], 'omitnan');
ylim([E_lo - 5e-4, 0.006]);
ax = gca;   ax.YAxis.Exponent = 0;   ytickformat('%.3f');
lgd = legend('Location', 'northwest');
set(lgd, 'Interpreter', 'latex');
title(sprintf('symmetry enforced, $\\epsilon = %.2f$, $d \\in [%.1f, %.1f]$', Para.epsilon, ss_d(1), ss_d(end)), ...
    'Interpreter', 'latex');
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
box on;
hold off;

subplot(1, 2, 2);
hold on;
for q = 1:3
    plot(mq, dmin(q, :), '.', 'Color', colors(q, :), 'MarkerSize', 8, 'DisplayName', fam{q});
end
xlabel('$m$', 'Interpreter', 'latex');
ylabel('effective spacing $d_{\mathrm{eff}}$ of the lowest stable state', 'Interpreter', 'latex');
xlim([0 1]);
ylim([1, ss_d(end) + 0.2]);
lgd = legend('Location', 'northwest');
set(lgd, 'Interpreter', 'latex');
title(sprintf('cells $d = %.1f$ to $%.1f$', ss_d(1), ss_d(end)), 'Interpreter', 'latex');
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
box on;
hold off;
exportgraphics(gcf, 'SpacingScan_E_vs_m.png', 'Resolution', 150);

%% the winning state of each family at a few masses (2 x 2 cells)
m_show = [0.15 0.3 0.45 0.6 0.75 0.9];
figure('Position', [50 50 1400 760], 'Color', 'w');
for q = 1:3
    for j = 1:numel(m_show)
        [~, i] = min(abs(mq - m_show(j)));
        subplot(3, numel(m_show), (q - 1)*numel(m_show) + j);
        if src(q, i, 1) == 0
            axis off;
            title(sprintf('%s, m = %.2f\nnone stable', fam{q}, m_show(j)), 'FontSize', 9);
            continue
        end
        a = src(q, i, 1);
        k = src(q, i, 2);
        if q == 1
            u = double(St{a}.u{k});
            ny = round(numel(u) * sqrt(3));
            u = repmat(u(:)', ny, 1);
            L = [ss_d(a), ss_d(a)*sqrt(3)];
        else
            u = double(H{a}.u{k});
            L = [ss_d(a), ss_d(a)*sqrt(3)];
        end
        imagesc([0 2*L(1)], [0 2*L(2)], repmat(u, 2, 2), [-1 1]);
        axis image off;
        title(sprintf('%s, m = %.2f\ncell %.1f, d_{eff} = %.2f, E = %.5f', fam{q}, m_show(j), ss_d(a), dmin(q, i), Emin(q, i)), 'FontSize', 9);
    end
end
exportgraphics(gcf, 'SpacingScan_states.png', 'Resolution', 110);

function tf = dark_percolates(v)
% true if the region v < 0 has a connected piece that reaches across a
% whole cell of the periodic pattern (a network); false if it is made of
% isolated domains (spots). Uses 4-neighbour connectivity on 3 x 3 copies
% of the cell, without wrapping.
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
g = graph([e1(:, 1); e2(:, 1)], [e1(:, 2); e2(:, 2)], [], numel(idx));
comp = conncomp(g);
% components with pixels in the central copy: does one span a full cell?
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

function kk = dominant_k(u, L)
% |k| of the largest nonzero Fourier mode of u on a cell of size L = [Lx Ly]
[ny, nx] = size(u);
A = abs(fft2(u - mean(u(:))));
A(1, 1) = 0;
[~, p] = max(A(:));
[r, c] = ind2sub([ny, nx], p);
i = mod(c - 1 + nx/2, nx) - nx/2;
j = mod(r - 1 + ny/2, ny) - ny/2;
kk = 2*pi*hypot(i/L(1), j/L(2));
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
