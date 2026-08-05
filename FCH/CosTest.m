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

L = 2*pi;
N = 128;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, 1);
m0 = 0.3;

K = 1:7;
branches = cell(numel(K), 1);

figure;
hold on;

for idx = 1:numel(K)
    k = K(idx);
    fprintf('\n\n===================== k = %d =====================\n', k);

    rng(k);     % reproducible per-k random seed
    u0 = m0 + 0.02*(rand(N, 1) - 0.5) + 0.1*cos(k*2*pi*Grid.x/L);
    u0 = u0 - mean(u0) + m0;    % force exact mass = m0

    try
        [u_vals, E_vals, m_vals, diagnostics] = arclength_FCH_func(u0, m0);

        branches{idx} = struct('k', k, 'u0', u0, 'u_vals', {u_vals}, ...
            'E_vals', E_vals, 'm_vals', m_vals, 'diagnostics', diagnostics);

        plot(m_vals, E_vals, 'LineWidth', 2, 'DisplayName', sprintf('k=%d', k));
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
plot(mm,Ec,'-.');
hold off;

save('CosTest_results.mat', 'branches', 'K', 'm0', 'N', 'L');
