% Polish every final state from an m-sweep with newton_polish_FCH, save the
% polished states, and plot them in one figure titled by m.

in_file = 'Ex5_test.mat';

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

fprintf('   m    Newton  |R| before   |R| after    lambda      E             time\n');
for i = 1:nm
    tic;
    [u_polished{i}, infos{i}] = newton_polish_FCH(S.final_u{i}, S.Grids{i}, Para, ...
        struct('verbose', false, 'maxit', 3000));   % slow coarsening can take ~1000 PTC steps
    info = infos{i};
    flag = '';
    if ~info.converged, flag = '  NOT CONVERGED'; end
    if info.collapsed, flag = [flag '  COLLAPSED']; end
    fprintf('  %.1f   %3d     %.2e     %.2e    %+.5f   %+.8f   %5.1fs%s\n', m_vals(i), ...
        numel(info.res_hist) - 1, info.res_hist(1), info.res_hist(end), info.lambda, ...
        info.E, toc, flag);
end

% drop the function handles before saving (they capture the whole workspace)
infos_save = cellfun(@(s) rmfield(s, {'mu', 'Jv'}), infos, 'UniformOutput', false);
[~, base] = fileparts(in_file);
out_file = [base '_polished.mat'];
save(out_file, 'm_vals', 'u_polished', 'infos_save', 'Para', 'in_file');
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
    title(sprintf('$m = %.1f$', m_vals(i)), 'Interpreter', 'latex')
end
