% Polish every final state from an m-sweep, check stability, and descend
% from saddles to stable states (descend_FCH); save the results and plot
% them in one figure titled by m.

% newest output of run_ex5_m_sweep.m (set in_file by hand to use another)
sweep_files = dir('ex5_m_sweep_*.mat');
sweep_files = sweep_files(~endsWith({sweep_files.name}, '_polished.mat'));   % skip our own output
[~, newest] = max([sweep_files.datenum]);
in_file = sweep_files(newest).name;
fprintf('Polishing %s\n', in_file);

% Para isn't saved in the sweep file; these match FCH_BDF2_SAV_ex5
Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;

S = load(in_file, 'final_u', 'Grids', 'm_vals');
m_vals = S.m_vals;
nm = numel(m_vals);

u_polished = cell(nm, 1);
infos = cell(nm, 1);
hists = cell(nm, 1);

% Polish, check stability, and descend from saddles to a stable state.
% even_only: stay with even states (stable = stable to even perturbations)
even_only = true;
fprintf('   m    dE first polish  neg   dE final         neg  rounds  time%s\n', '');
for i = 1:nm
    tic;
    [u_polished{i}, hists{i}, infos{i}] = descend_FCH(S.final_u{i}, S.Grids{i}, Para, ...
        struct('verbose', false, 'even_only', even_only));
    h = hists{i};
    flag = '';
    if ~infos{i}.converged, flag = '  NOT CONVERGED'; end
    if h(end).n_neg > 0, flag = [flag '  STILL UNSTABLE']; end
    if infos{i}.collapsed, flag = [flag '  COLLAPSED']; end
    fprintf('  %.2f  %+.10f  %2d    %+.10f  %2d   %3d   %5.1fs%s\n', m_vals(i), ...
        h(1).dE, h(1).n_neg, h(end).dE, h(end).n_neg, numel(h) - 1, toc, flag);
end

% drop the function handles before saving (they capture the whole workspace)
infos_save = cellfun(@(s) rmfield(s, {'mu', 'Jv'}), infos, 'UniformOutput', false);
[~, base] = fileparts(in_file);
out_file = [base '_polished.mat'];
save(out_file, 'm_vals', 'u_polished', 'infos_save', 'hists', 'Para', 'in_file', 'even_only');
fprintf('Saved polished states to %s\n', out_file);

%% Plot all polished states
figure;
tiledlayout('flow', 'TileSpacing', 'compact');
for i = 1:nm
    nexttile;
    Grid = S.Grids{i};
    if isvector(u_polished{i})
        plot(Grid.x, u_polished{i}, 'LineWidth', 2)
        xlim([Grid.x(1) Grid.x(end)])
    else
        surf(Grid.xx, Grid.yy, u_polished{i}), shading interp, view(2);
        axis equal tight off;
        colorbar;
    end
    title(sprintf('$m = %.2f$, $\\Delta E = %.6f$', m_vals(i), infos{i}.dE), 'Interpreter', 'latex')
end
