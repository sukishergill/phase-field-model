% Continue every polished state of an m-sweep in m (both directions) with
% continuation_FCH, up to the first fold or until m leaves [0, 1], and
% assemble the branches. Plot / compare with plot_m_sweep_branches.m,
% compare_to_references.m, lowest_family_by_m.m.
%
% Each branch is saved to its own file in out_dir as soon as it finishes,
% and branches that already have a file are skipped, so an interrupted run
% resumes where it stopped. Progress prints in the command window.
%
% Cost per point (even mode, direct solver, 128^2): ~4-7 s on the client,
% roughly 3x that on a single-threaded parallel worker.

% newest polished sweep (set in_file by hand to use another)
pol_files = dir('ex5_m_sweep_*_polished.mat');
[~, newest] = max([pol_files.datenum]);
in_file = pol_files(newest).name;
[~, base] = fileparts(in_file);
out_dir = 'branches_tofold';
fprintf('Continuing states from %s into %s/\n', in_file, out_dir);

PL = load(in_file, 'u_polished', 'm_vals', 'Para', 'even_only');
SW = load(strrep(in_file, '_polished', ''), 'Grids');
Para = PL.Para;
nm = numel(PL.m_vals);

copts = struct('even', PL.even_only, 'max_steps', 5000, 'stop_at_fold', true, ...
    'm_range', [0 1], 'verbose', false, 'store_u', true, 'polish_first', false);

dirs = [1 -1];
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
branch_file = @(i, d) fullfile(out_dir, sprintf('%s_m%.2f_dir%+d.mat', base, PL.m_vals(i), dirs(d)));

njobs = nm * numel(dirs);
done = false(1, njobs);
for job = 1:njobs
    [i, d] = ind2sub([nm, numel(dirs)], job);
    done(job) = exist(branch_file(i, d), 'file') == 2;
end
todo = find(~done);
fprintf('%d of %d branches already done; running %d\n', sum(done), njobs, numel(todo));

% progress in the command window as each branch finishes
q = parallel.pool.DataQueue;
t_start = tic;
afterEach(q, @(msg) fprintf('[%6.1f min] %s\n', toc(t_start)/60, msg));
parfor jj = 1:numel(todo)
    [i, d] = ind2sub([nm, numel(dirs)], todo(jj));
    o = copts;
    o.direction = dirs(d);
    b = continuation_FCH(PL.u_polished{i}, SW.Grids{i}, Para, o);
    save_branch(branch_file(i, d), b);
    send(q, sprintf('done: start m = %.2f, dir %+d, %d points, m -> %.4f (%s)', PL.m_vals(i), ...
        dirs(d), numel(b.m), b.m(end), b.stop_reason));
end

% assemble
branches = cell(nm, numel(dirs));
for job = 1:njobs
    [i, d] = ind2sub([nm, numel(dirs)], job);
    S = load(branch_file(i, d), 'b');
    branches{i, d} = S.b;
end
for i = 1:nm
    for d = 1:numel(dirs)
        b = branches{i, d};
        fprintf('start m = %.2f, dir %+d: %4d points, m in [%.4f, %.4f]: %s\n', ...
            PL.m_vals(i), dirs(d), numel(b.m), min(b.m), max(b.m), b.stop_reason);
    end
end
out_file = [base '_tofold_branches.mat'];
save(out_file, 'branches', 'in_file', 'copts', 'dirs', '-v7.3');
fprintf('Saved branches to %s\n', out_file);

function save_branch(file, b)
% save inside parfor (save cannot be called directly in a parfor body)
save(file, 'b', '-v7.3');
end
