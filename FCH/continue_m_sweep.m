% Continue every polished state of an m-sweep in m (both directions) with
% continuation_FCH, save the branches, and plot dE = E - E_uniform(m)
% against m for all of them (solid stable, dashed unstable).
%
% Cost per point (even mode, direct solver): ~0.3 s at 64^2, ~6 s at
% 128^2 near soft (nearly marginal) states. Branches run in parallel.

% newest polished sweep (set in_file by hand to use another)
pol_files = dir('ex5_m_sweep_*_polished.mat');
[~, newest] = max([pol_files.datenum]);
in_file = pol_files(newest).name;
fprintf('Continuing states from %s\n', in_file);

PL = load(in_file, 'u_polished', 'm_vals', 'Para', 'even_only');
SW = load(strrep(in_file, '_polished', ''), 'Grids');
Para = PL.Para;
nm = numel(PL.m_vals);

copts = struct('even', PL.even_only, 'max_steps', 150, 'verbose', false, ...
    'store_u', true, 'polish_first', false);

dirs = [1 -1];
branches = cell(nm, numel(dirs));
parfor job = 1:nm*numel(dirs)
    [i, d] = ind2sub([nm, numel(dirs)], job);
    o = copts;
    o.direction = dirs(d);
    branches{job} = continuation_FCH(PL.u_polished{i}, SW.Grids{i}, Para, o);
end
branches = reshape(branches, nm, numel(dirs));

for i = 1:nm
    for d = 1:numel(dirs)
        b = branches{i, d};
        fprintf('start m = %.2f, dir %+d: %3d points, m in [%.3f, %.3f], %d fold(s): %s\n', ...
            PL.m_vals(i), dirs(d), numel(b.m), min(b.m), max(b.m), sum(b.fold), b.stop_reason);
    end
end

[~, base] = fileparts(in_file);
out_file = [base '_branches.mat'];
save(out_file, 'branches', 'in_file', 'copts', 'dirs', '-v7.3');
fprintf('Saved branches to %s\n', out_file);

%% Plot: one colour per starting state, both directions
figure;
ax = gca;
plot_branch_FCH(reshape(branches', 1, []), ax);
