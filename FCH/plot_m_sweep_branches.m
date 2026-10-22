% Plot dE = E - E_uniform(m) against m for every branch from
% continue_m_sweep.m: one colour per starting state (both continuation
% directions in the same colour), solid stable, dashed unstable, dots at
% folds. If reference_patterns.mat exists (reference_patterns_sweep.m),
% the identical-stripe and hex-spot-lattice branches are overlaid as thin
% lines (solid stable, dotted unstable within their own symmetry class).
% Saves a .png next to the branches file.

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
labels = arrayfun(@(m) sprintf('start m = %.2f', m), m0, 'UniformOutput', false);

% reference patterns (reference_patterns_sweep.m), if computed: n identical
% stripes (gray) and ideal hex spot lattices (magenta), thin lines
xl = [min(cellfun(@(b) min(b.m), br_cells(:))) - 0.05, max(cellfun(@(b) max(b.m), br_cells(:))) + 0.05];
if ~isempty(dir('reference_patterns*.mat'))
    Rf = load_reference_patterns();     % all reference_patterns*.mat, merged
    fam = {Rf.stripes, Rf.spots};
    fam_col = {[0.45 0.45 0.45], [0.75 0 0.75]};
    fam_lab = {@(a) sprintf('%d', Rf.n_stripes(a)), @(a) sprintf('d=%.2f', Rf.spot_d(a))};
    fam_name = {'n identical stripes (label n)', 'hex spot lattice (label d)'};
    labelled = {{}, {}};
    for F = 1:2
        for a = 1:size(fam{F}, 1)
            best = [inf, NaN];
            for d = 1:size(fam{F}, 2)
                b = fam{F}{a, d};
                stable = b.n_neg == 0;
                edges = [0, find(diff(stable) ~= 0), numel(stable)];
                for r = 1:numel(edges) - 1
                    idx = edges(r) + 1 : min(edges(r+1) + 1, numel(stable));
                    if stable(edges(r) + 1), style = '-'; else, style = ':'; end
                    hr = plot(ax, b.m(idx), b.dE(idx), style, 'Color', fam_col{F}, 'LineWidth', 0.9);
                end
                dEin = b.dE;
                dEin(b.m < xl(1) | b.m > xl(2)) = inf;   % label within the plotted range
                [emin, kmin] = min(dEin);
                if emin < best(1), best = [emin, b.m(kmin)]; end
            end
            if isfinite(best(1)) && ~ismember(fam_lab{F}(a), labelled{F})
                text(ax, best(2), best(1), fam_lab{F}(a), 'Color', fam_col{F}, 'FontSize', 8, ...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
                labelled{F}{end+1} = fam_lab{F}(a);
            end
        end
        h(end+1) = hr; %#ok<AGROW>
        labels{end+1} = fam_name{F}; %#ok<AGROW>
    end
    xlim(ax, xl);
end
legend(ax, h, labels, 'Location', 'southwest');
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
