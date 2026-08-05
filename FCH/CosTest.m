% CosTest.m
%
% Builds initial profiles
%   u0 = m0 + small random noise + 0.1*cos(k*2*pi*x/L)
% for k = 1, ..., 7 at mass m0 = 0.3, runs the FCH pseudo-arclength
% continuation (arclength_FCH_func) starting from each, and plots all the
% resulting m-E branches together for comparison.
%
% Each k gets its own reproducible random seed (rng(k)) so re-running this
% script reproduces the same starting profiles every time. Mass is forced
% to exactly m0 after adding the noise + cosine (their combination won't
% generally average to exactly zero over a finite grid, and any mismatch
% between the intended and actual starting mass gets silently carried
% through the whole run -- see the mean-centering discussion earlier).

L = 4*pi;
N = 128;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, 1);
m0 = 0.3;
eps = 0.1;      eta = eps^2;

K = 1:11;
branches = cell(numel(K), 1);

figure;
hold on;
colors = get(gca, 'ColorOrder');

for idx = 1:numel(K)
    k = K(idx);
    fprintf('\n\n===================== k = %d =====================\n', k);

    rng(k);     % reproducible per-k random seed
    u0 = m0 + 0.00*(rand(N, 1) - 0.5) + .95*(1-m0)*cos(k*2*pi*Grid.x/L);
    u0 = u0 - mean(u0) + m0;    % force exact mass = m0

    try
        [u_vals, E_vals, m_vals, diagnostics] = arclength_FCH_func(u0, m0);

        branches{idx} = struct('k', k, 'u0', u0, 'u_vals', {u_vals}, ...
            'E_vals', E_vals, 'm_vals', m_vals, 'diagnostics', diagnostics);

        color = colors(mod(idx-1, size(colors,1)) + 1, :);
        plot_branch_by_stability(m_vals, E_vals, diagnostics.stability, ...
            color, sprintf('k=%d', k));
    catch err
        fprintf(2, 'k=%d FAILED: %s\n', k, err.message);
        branches{idx} = struct('k', k, 'u0', u0, 'error', err.message);
    end
end

xlabel('$m$', 'Interpreter', 'latex');
ylabel('$E$', 'Interpreter', 'latex');
lgd = legend('Location', 'best');
set(lgd, 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
set(gca, 'TickLabelInterpreter', 'latex');
title('FCH branches from $u_0 = m_0 + \mathrm{noise} + 0.1\cos(k \cdot 2\pi x/L)$', ...
    'Interpreter', 'latex');
mm = linspace(0,1,1001);
Ec = (mm-mm.^3).^2 - eta*(1-mm.^2).^2/4;
hold on;
plot(mm,Ec,'-.', 'DisplayName', 'u=m');
hold off;

save('CosTest_results.mat', 'branches', 'K', 'm0', 'N', 'L','eta','eps');

function plot_branch_by_stability(m_vals, E_vals, stability, color, dispname)
% Plot a single branch, solid where stable (stability==+1) and dashed
% where unstable (stability==-1). Consecutive same-sign runs are drawn as
% separate line segments (all in the same color), each extended one point
% into its neighbour so the curve reads as continuous. Only the first
% segment gets a legend entry.

n = numel(m_vals);
first_seg = true;
i = 1;
while i <= n
    j = i;
    while j < n && stability(j+1) == stability(i)
        j = j + 1;
    end
    idx_end = min(j+1, n);    % overlap one point with the next segment

    if stability(i) > 0
        style = '-';
    else
        style = '--';
    end

    if first_seg
        plot(m_vals(i:idx_end), E_vals(i:idx_end), style, 'Color', color, ...
            'LineWidth', 2, 'DisplayName', dispname);
        first_seg = false;
    else
        plot(m_vals(i:idx_end), E_vals(i:idx_end), style, 'Color', color, ...
            'LineWidth', 2, 'HandleVisibility', 'off');
    end

    i = j + 1;
end
end
