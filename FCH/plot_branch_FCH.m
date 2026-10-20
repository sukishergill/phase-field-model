function plot_branch_FCH(branches, ax, yfield)
% Plot branches from continuation_FCH against m: solid where stable
% (n_neg = 0), dashed where unstable, circles at folds, squares at possible
% branch points.
%
%   plot_branch_FCH(branch)               % dE = E - E_uniform(m) vs m
%   plot_branch_FCH({b1, b2, ...})        % several, one colour each
%   plot_branch_FCH(branch, ax)           % into existing axes ([] = new figure)
%   plot_branch_FCH(branch, ax, 'lambda') % chemical potential vs m
%
% yfield: 'dE' (default), 'lambda' or 'E'. Sheets of a branch can be
% nearly on top of each other in dE; lambda = dE_total/dm separates them
% and shows each fold as a turning point.

if ~iscell(branches), branches = {branches}; end
if nargin < 2 || isempty(ax)
    figure;
    ax = gca;
end
if nargin < 3, yfield = 'dE'; end
hold(ax, 'on');
colors = lines(numel(branches));
for b = 1:numel(branches)
    br = branches{b};
    c = colors(b, :);
    stable = br.n_neg == 0;
    % split into runs of equal stability; overlap by one point so the
    % pieces join up
    edges = [0, find(diff(stable) ~= 0), numel(stable)];
    for r = 1:numel(edges) - 1
        idx = edges(r) + 1 : min(edges(r+1) + 1, numel(stable));
        if stable(edges(r) + 1), style = '-'; else, style = '--'; end
        plot(ax, br.m(idx), br.(yfield)(idx), style, 'Color', c, 'LineWidth', 2);
    end
    plot(ax, br.m(br.fold), br.(yfield)(br.fold), 'o', 'Color', c, 'MarkerFaceColor', c);
    plot(ax, br.m(br.bif), br.(yfield)(br.bif), 's', 'Color', c, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
end
xlabel(ax, '$m$', 'Interpreter', 'latex');
ylabels = struct('dE', '$E - E_{\mathrm{uniform}}(m)$', 'lambda', '$\lambda$ (chemical potential)', 'E', '$E$');
ylabel(ax, ylabels.(yfield), 'Interpreter', 'latex');
title(ax, 'solid: stable, dashed: unstable, $\bullet$ fold, $\square$ possible branch point', ...
    'Interpreter', 'latex');
box(ax, 'on');
end
