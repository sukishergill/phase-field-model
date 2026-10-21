% Plot dE = E - E_uniform(m) against m for every branch from
% continue_m_sweep.m: one colour per starting state (both continuation
% directions in the same colour), solid stable, dashed unstable, dots at
% folds. Saves a .png next to the branches file.

% newest branches file (set br_file by hand to use another)
br_files = dir('ex5_m_sweep_*_branches.mat');
[~, newest] = max([br_files.datenum]);
br_file = br_files(newest).name;
B = load(br_file, 'branches');
br_cells = B.branches;                     % nm x 2 (directions +1, -1)
[nm, nd] = size(br_cells);
m0 = cellfun(@(b) b.m(1), br_cells(:, 1));

f = figure('Position', [100 100 1100 750]);
ax = axes(f);
hold(ax, 'on');
colors = turbo(nm + 2);
colors = colors(2:end-1, :);
h = gobjects(nm, 1);
for i = 1:nm
    c = colors(i, :);
    for d = 1:nd
        b = br_cells{i, d};
        stable = b.n_neg == 0;
        edges = [0, find(diff(stable) ~= 0), numel(stable)];
        for r = 1:numel(edges) - 1
            idx = edges(r) + 1 : min(edges(r+1) + 1, numel(stable));
            if stable(edges(r) + 1), style = '-'; else, style = '--'; end
            hl = plot(ax, b.m(idx), b.dE(idx), style, 'Color', c, 'LineWidth', 1.8);
            if d == 1 && r == 1, h(i) = hl; end
        end
        plot(ax, b.m(b.fold), b.dE(b.fold), 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4);
    end
    plot(ax, m0(i), br_cells{i, 1}.dE(1), 'p', 'Color', 'k', 'MarkerFaceColor', c, 'MarkerSize', 12);
end
legend(ax, h, arrayfun(@(m) sprintf('start m = %.2f', m), m0, 'UniformOutput', false), ...
    'Location', 'southwest');
xlabel(ax, '$m$', 'Interpreter', 'latex');
ylabel(ax, '$E - E_{\mathrm{uniform}}(m)$', 'Interpreter', 'latex');
G = br_cells{1, 1}.Grid;
title(ax, sprintf('All branches, %dx%d, L = %.2f: solid stable, dashed unstable, dots folds, stars starting states', ...
    G.N(1), G.N(end), G.L(1)), 'Interpreter', 'none');
box(ax, 'on');
grid(ax, 'on');

png = strrep(br_file, '.mat', '_dE_vs_m.png');
exportgraphics(f, png, 'Resolution', 110);
fprintf('Saved %s\n', png);
