% Checks for continuation_FCH. Prints PASS/FAIL per check.
%   1. every branch point is a converged steady state with mean(u) = m
%   2. re-polishing branch points independently (plain Newton) gives the
%      same energy
%   3. n_neg changes by exactly 1 at each fold (exchange of stability)
%   4. direct (even) and GMRES solvers trace the same branch
%   5. 2D (64^2, even): fold detected with n_neg changing by 1

clear
pass = @(ok, name) fprintf('%s  %s\n', ternary(ok, 'PASS', 'FAIL'), name);

P.epsilon = 0.18;
P.eta1 = P.epsilon^2;
P.eta2 = P.eta1;
P.tau = 0;

% 1D even bilayer pair, mean ~ 0.36
G1 = SSAV_FCH_helpers.generate_Grid(2*pi, 128, 1);
u1 = -(tanh((G1.x + 1)/0.18) - tanh((G1.x - 1)/0.18) - 1);

b = continuation_FCH(u1, G1, P, struct('even', true, 'max_steps', 120, 'verbose', false));

%% 1. converged, mass correct
res = cellfun(@(u) norm(reshape(SSAV_FCH_helpers.steady_mu(u, P, G1) - mean(SSAV_FCH_helpers.steady_mu(u, P, G1)), [], 1)) / sqrt(numel(u)), b.u);
mass_err = max(abs(cellfun(@mean, b.u) - b.m));
pass(max(res) < 1e-9 && mass_err < 1e-12, sprintf('1D: %d points, max |R| %.1e, max |mean(u) - m| %.1e', numel(b.m), max(res), mass_err));

%% 2. independent re-polish at fixed m
idx = round(linspace(2, numel(b.m), 6));
dE = zeros(size(idx));
for j = 1:numel(idx)
    [~, ij] = newton_polish_FCH(b.u{idx(j)}, G1, P, struct('ptc', false, 'even', true, 'verbose', false));
    dE(j) = abs(ij.E - b.E(idx(j)));
end
pass(max(dE) < 1e-12, sprintf('1D: re-polished branch points, max |E diff| %.1e', max(dE)));

%% 3. stability exchange at folds
kf = find(b.fold);
jumps = abs(b.n_neg(kf) - b.n_neg(kf - 1));
pass(~isempty(kf) && all(jumps == 1), sprintf('1D: %d fold(s), n_neg jumps %s', numel(kf), mat2str(jumps)));

%% 4. direct vs GMRES: same branch (compare dE(m) before the first fold)
o = struct('even', true, 'max_steps', 30, 'verbose', false, 'ds_max', 0.01);
bd = continuation_FCH(u1, G1, P, setfield(o, 'linear_solver', 'direct')); %#ok<SFLD>
bg = continuation_FCH(u1, G1, P, setfield(o, 'linear_solver', 'gmres'));  %#ok<SFLD>
fd = find(bd.fold, 1);  if isempty(fd), fd = numel(bd.m) + 1; end
fg = find(bg.fold, 1);  if isempty(fg), fg = numel(bg.m) + 1; end
md = bd.m(1:fd-1);  Ed = bd.dE(1:fd-1);
mg = bg.m(1:fg-1);  Eg = bg.dE(1:fg-1);
in = mg >= min(md) & mg <= max(md);
diffE = max(abs(interp1(md, Ed, mg(in), 'pchip') - Eg(in)));
pass(nnz(in) >= 5 && diffE < 1e-6, sprintf('direct vs GMRES: %d shared points, max |dE diff| %.1e', nnz(in), diffE));

%% 5. 2D, 64^2, even: an m = 0.3 sweep state folds near m = 0.347
G2 = SSAV_FCH_helpers.generate_Grid([2*pi 2*pi], [64 64], 2);
S = load('Ex5_test_polished.mat', 'u_polished');
b2 = continuation_FCH(S.u_polished{4}, G2, P, struct('even', true, 'max_steps', 40, 'verbose', false));
kf2 = find(b2.fold);
pass(~isempty(kf2) && all(abs(b2.n_neg(kf2) - b2.n_neg(kf2 - 1)) == 1) && max(b2.res) < 1e-9, ...
    sprintf('2D 64^2: fold(s) at m = %s, n_neg %s -> %s, max |R| %.1e', mat2str(b2.m(kf2), 4), ...
    mat2str(b2.n_neg(kf2 - 1)), mat2str(b2.n_neg(kf2)), max(b2.res)));

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
