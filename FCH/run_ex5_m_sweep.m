% Run FCH_BDF2_SAV_ex5 for each mass m in m_vals (in parallel) and plot the
% final profile of every run in one figure.

m_vals = 0.05:0.1:0.85;    % m = 0 avoided: its unstable band is a narrow sliver
nm = numel(m_vals);

final_u = cell(nm, 1);
Grids = cell(nm, 1);
energy = cell(nm, 1);       % per run: t_vals, Eu (SAV); t_direct, Eu_direct (true energy, sampled)
Times = cell(nm, 1);

parfor i = 1:nm
    [final_u{i}, Grids{i}, energy{i}, Times{i}] = run_ex5(m_vals(i));
end

% Save before plotting so a plotting error can't lose the runs. Timestamped
% so a new sweep never overwrites an old one.
out_file = sprintf('ex5_m_sweep_%s.mat', datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
save(out_file, 'm_vals', 'final_u', 'Grids', 'energy', 'Times');
fprintf('Saved results to %s\n', out_file);

%% Plot all final profiles
figure;
tiledlayout('flow', 'TileSpacing', 'compact');
for i = 1:nm
    nexttile;
    Grid = Grids{i};
    if isvector(final_u{i})
        plot(Grid.x, final_u{i}, 'LineWidth', 2)
        xlim([Grid.x(1) Grid.x(end)])
    else
        surf(Grid.xx, Grid.yy, final_u{i}), shading interp, view(2);
        axis equal tight off;
        colorbar;
    end
    title(sprintf('$m = %.2f$', m_vals(i)), 'Interpreter', 'latex')
end

function [u_final, Grid, energy, Time] = run_ex5(m_override) %#ok<INUSD> read by the script via exist()
% parfor can't call a script directly, so run it inside this function's
% workspace with m_override set.
FCH_BDF2_SAV_ex5;
u_final = Results.uu{end};
energy.t_vals = Results.t_vals;
energy.Eu = Results.Eu;
energy.Eu_direct = Results.Eu_direct;
energy.t_direct = Results.t_direct;
end
