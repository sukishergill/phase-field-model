% Which pattern family has the lowest dE = E - E_uniform(m) at each mass:
%   hex spots  hex lattices of dark spots (any d) and of light domains with
%              d < 4                     (reference_patterns*.mat)
%   honeycomb  hex lattices of light domains with d >= 4 (membrane cells),
%              best stable spacing at each scanned m, cells that fit the
%              box                        (honeycomb_scan.mat, interpolated)
%   stripes    n identical stripes        (reference_patterns*.mat)
%   other      the m-sweep branches       (all *_tofold_branches.mat, merged)
% The spots / honeycomb split at d = 4 is a convention: light lattices turn
% continuously into honeycombs as d grows.
% Prints the m-intervals each family wins, with the margin to the runner-up,
% and plots the lower envelopes. "other" only exists where the sweep
% branches reach, so outside those m it cannot win by construction; a sweep
% state within tie_tol of a reference counts as that reference pattern.

stable_only = true;     % false: also count unstable stretches of every family
tie_tol = 5e-6;         % "other" within this of a reference is that reference pattern (covers interpolation error between branch points)
dm = 0.001;
mq = 0:dm:0.95;         % tau = 0: E is symmetric under u -> -u, m -> -m

B.branches = load_sweep_branches();     % all continue_m_sweep outputs, merged
Rf = load_reference_patterns();         % all reference_patterns*.mat, merged

names = {'hex spots', 'honeycomb', 'stripes', 'other'};
nf = numel(names);
env = inf(nf, numel(mq));
who = nan(nf, numel(mq));
[env(1, :), who(1, :)] = lowest_dE_at_m(Rf.spots, mq, stable_only);
[env(3, :), who(3, :)] = lowest_dE_at_m(Rf.stripes, mq, stable_only);
[env(4, :), who(4, :)] = lowest_dE_at_m(B.branches, mq, stable_only);
% honeycomb: best stable spacing from honeycomb_scan, linear in m between
% scanned masses (only where both neighbouring scan points have one)
HC = struct('hc_m', [], 'E_honey', [], 'd_honey', []);
if exist('honeycomb_scan.mat', 'file')
    HC = load('honeycomb_scan.mat', 'hc_m', 'E_honey', 'd_honey');
    for q = 1:numel(mq)
        a = find(HC.hc_m <= mq(q), 1, 'last');
        if isempty(a) || a == numel(HC.hc_m)
            if ~isempty(a) && abs(HC.hc_m(a) - mq(q)) < 1e-12, env(2, q) = HC.E_honey(a); who(2, q) = a; end
            continue
        end
        t = (mq(q) - HC.hc_m(a)) / (HC.hc_m(a+1) - HC.hc_m(a));
        env(2, q) = (1 - t)*HC.E_honey(a) + t*HC.E_honey(a+1);     % inf if either is inf
        who(2, q) = a + (t > 0.5);
    end
    env(2, isnan(env(2, :))) = inf;
end
% label of the winning member: d for spots / honeycomb, n for stripes, start m for sweep
member = @(F, w) member_label(F, w, Rf, B, HC);

% a sweep state within tie_tol of a reference is the same pattern (e.g.
% the m = 0.85 sweep state is a single stripe): credit the reference
same = abs(env(4, :) - min(env(1:3, :), [], 1)) < tie_tol;
env(4, same) = inf;
who(4, same) = NaN;
[best, win] = min(env, [], 1);
win(~isfinite(best)) = 0;
sorted = sort(env, 1);
margin = sorted(2, :) - sorted(1, :);            % gap to the runner-up (inf if alone)

fprintf('Lowest-energy family by m (stable_only = %d, dm = %.3f)\n', stable_only, dm);
fprintf('  %-10s  %-17s  %-24s  %s\n', 'family', 'm range', 'winning member(s)', 'margin to runner-up (min..max)');
k = 1;
while k <= numel(mq)
    j = k;
    while j < numel(mq) && win(j+1) == win(k)
        j = j + 1;
    end
    if win(k) == 0
        fprintf('  %-10s  [%+.3f, %+.3f]\n', '(none)', mq(k), mq(j));
    else
        mem = unique(arrayfun(@(q) string(member(win(k), who(win(k), q))), k:j), 'stable');
        mg = margin(k:j);
        mg = mg(isfinite(mg));
        if isempty(mg), mtxt = 'only family with data';
        else, mtxt = sprintf('%.1e .. %.1e', min(mg), max(mg));
        end
        fprintf('  %-10s  [%+.3f, %+.3f]  %-24s  %s\n', names{win(k)}, mq(k), mq(j), strjoin(mem, ', '), mtxt);
    end
    k = j + 1;
end
fprintf('"other" (sweep) branches cover m in: %s\n', intervals_str(mq, ~isinf(env(4, :)) | same));
if exist('HC', 'var')
    fprintf('honeycomb scanned at m = %.2f..%.2f (no stable honeycomb at m = %s)\n', HC.hc_m(1), HC.hc_m(end), ...
        mat2str(HC.hc_m(~isfinite(HC.E_honey))));
end

%% plot: the lower envelopes, the overall minimum marked by family
cols = [0.75 0 0.75; 0.1 0.6 0.6; 0.45 0.45 0.45; 0.85 0.33 0.1];
f = figure('Position', [100 100 1100 650]);
ax = axes(f); hold(ax, 'on');
h = gobjects(nf, 1);
for F = 1:nf
    y = env(F, :); y(~isfinite(y)) = NaN;
    h(F) = plot(ax, mq, y, '-', 'Color', cols(F, :), 'LineWidth', 1.2);
end
for F = 1:nf
    y = best; y(win ~= F) = NaN;
    plot(ax, mq, y, '-', 'Color', cols(F, :), 'LineWidth', 5);
end
legend(ax, h, strcat(names, ' (lowest member)'), 'Location', 'southwest');
xlabel(ax, '$m$', 'Interpreter', 'latex');
ylabel(ax, 'lowest $E - E_{\mathrm{uniform}}(m)$', 'Interpreter', 'latex');
title(ax, 'Lowest energy per family; thick = overall lowest at that m', 'Interpreter', 'none');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [mq(1) mq(end)]);
exportgraphics(f, 'lowest_family_by_m.png', 'Resolution', 110);
fprintf('Saved lowest_family_by_m.png\n');

function s = member_label(F, w, Rf, B, HC)
if isnan(w), s = '-'; return, end
% w is a linear index into the family's (member x direction) cell array
% (for the honeycomb: an index into the scanned masses)
switch F
    case 1
        i = mod(w - 1, size(Rf.spots, 1)) + 1;     % row = member
        s = Rf.spot_label{i};
    case 2
        s = sprintf('d=%.2f', HC.d_honey(w));
    case 3
        i = mod(w - 1, size(Rf.stripes, 1)) + 1;
        s = sprintf('n=%d', Rf.n_stripes(i));
    case 4
        i = mod(w - 1, size(B.branches, 1)) + 1;
        s = sprintf('start %.2f', B.branches{i, 1}.m(1));
end
end

function s = intervals_str(x, mask)
d = diff([0 mask 0]);
st = find(d == 1); en = find(d == -1) - 1;
s = strjoin(arrayfun(@(a, b) sprintf('[%.3f, %.3f]', x(a), x(b)), st, en, 'UniformOutput', false), ' ');
end
