% For each starting state of the m-sweep, compare its dE = E - E_uniform(m)
% with the lowest stable reference pattern at the same m: n identical
% stripes and ideal hex spot lattices (reference_patterns_sweep.m).
% Reference stability is within their own symmetry class (1D one-period
% for stripes, cell-periodic for spots); interpolated linearly in m
% between neighbouring branch points. Set stable_only = false to compare
% against every reference state, stable or not (any simple-pattern
% steady state with lower energy counts against a 'complex state is
% lower' claim).

if ~exist('stable_only', 'var'), stable_only = true; end
br_files = dir('ex5_m_sweep_*_branches.mat');
[~, newest] = max([br_files.datenum]);
B = load(br_files(newest).name, 'branches');
Rf = load('reference_patterns.mat', 'stripes', 'n_stripes', 'spots', 'spot_d');

fprintf('%s vs reference_patterns.mat (stable_only = %d)\n', br_files(newest).name, stable_only);
fprintf(' start m   sweep dE     best stripes dE (n)   best hex spots dE (d)   sweep - best ref\n');
for i = 1:size(B.branches, 1)
    b0 = B.branches{i, 1};
    mq = b0.m(1);
    e0 = b0.dE(1);
    [es, ks] = best_ref(Rf.stripes, mq, stable_only);
    [eh, kh] = best_ref(Rf.spots, mq, stable_only);
    lab_s = '-';  if ~isnan(ks), lab_s = sprintf('%d', Rf.n_stripes(ks)); end
    lab_h = '-';  if ~isnan(kh), lab_h = sprintf('%.2f', Rf.spot_d(kh)); end
    if e0 < min(es, eh), verdict = 'sweep state lower'; else, verdict = 'reference lower'; end
    fprintf('  %.2f   %+.6f   %+.6f (n=%s)       %+.6f (d=%s)        %+.6f  %s\n', ...
        mq, e0, es, lab_s, eh, lab_h, e0 - min(es, eh), verdict);
end

function [e, which] = best_ref(fam, mq, stable_only)
% lowest dE at m = mq over all branches of a family (stable stretches only
% if stable_only)
e = inf;
which = NaN;
for a = 1:size(fam, 1)
    for d = 1:size(fam, 2)
        b = fam{a, d};
        for k = 1:numel(b.m) - 1
            ok = ~stable_only || (b.n_neg(k) == 0 && b.n_neg(k+1) == 0);
            if ok && b.m(k) ~= b.m(k+1) && ...
                    (b.m(k) - mq)*(b.m(k+1) - mq) <= 0
                t = (mq - b.m(k)) / (b.m(k+1) - b.m(k));
                ek = b.dE(k) + t*(b.dE(k+1) - b.dE(k));
                if ek < e
                    e = ek;
                    which = a;
                end
            end
        end
    end
end
end
