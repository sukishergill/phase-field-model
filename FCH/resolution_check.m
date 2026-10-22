% Grid-resolution check of the energy comparison at a fixed mass m_target:
% the lowest "other" (m-sweep branch) state against the lowest n-stripe and
% hex-lattice states at the same m. Each is taken from the nearest stable
% point on its branch, shifted to exactly m_target, polished there, then
% Fourier-interpolated to finer grids and re-polished. If each dE converges
% and the gaps keep their sign and size, the comparison is resolved.
%
% Settings (override before running):
%   m_target        (0.35)
%   other_factors   refinement factors for the sweep state ([1 2]: 128^2, 256^2)
%   stripe_factors  for the one-period stripe cell   ([1 2 4 8])
%   hex_factors     for the two-spot hex cell        ([1 2 4])

if ~exist('m_target', 'var'), m_target = 0.35; end
if ~exist('other_factors', 'var'), other_factors = [1 2]; end
if ~exist('stripe_factors', 'var'), stripe_factors = [1 2 4 8]; end
if ~exist('hex_factors', 'var'), hex_factors = [1 2 4]; end

[B.branches, br_names] = load_sweep_branches();   % all continue_m_sweep outputs, merged
Rf = load_reference_patterns();
Para = B.branches{1}.Para;

fams = {'other', 'stripes', 'hex'};
sets = {B.branches, Rf.stripes, Rf.spots};
factors = {other_factors, stripe_factors, hex_factors};
labels = {@(w) sprintf('start m %.2f', B.branches{mod(w-1, size(B.branches, 1)) + 1, 1}.m(1)), ...
          @(w) sprintf('n = %d', Rf.n_stripes(mod(w-1, size(Rf.stripes, 1)) + 1)), ...
          @(w) sprintf('d = %.2f', Rf.spot_d(mod(w-1, size(Rf.spots, 1)) + 1))};
fprintf('resolution check at m = %.4f (%s)\n', m_target, strjoin(br_names, ' + '));

dE = cell(1, 3);
for F = 1:3
    [e0, w] = lowest_dE_at_m(sets{F}, m_target, true);
    if ~isfinite(e0)
        fprintf('%-8s no stable state at this m\n', fams{F});
        continue
    end
    b = sets{F}{w};
    % nearest stable point on that branch, shifted to exactly m_target
    [~, k] = min(abs(b.m - m_target) + 1e3*(b.n_neg ~= 0 | isnan(b.dE)));
    u0 = b.u{k} + (m_target - b.m(k));
    G0 = b.Grid;
    popts = struct('even', true, 'verbose', false, 'maxit', 3000);
    if isfield(b.opts, 'symmetry') && ~isempty(b.opts.symmetry)
        % The stored symmetry has the coarse grid's shift baked in (the hex
        % cell's half-cell shift). Use the size-generic half-shift instead,
        % after checking it is the same operation on the coarse grid.
        half_shift = @(v) (v + circshift(v, size(v)/2)) / 2;
        r = randn(size(u0));
        assert(max(abs(half_shift(r) - b.opts.symmetry(r)), [], 'all') < 1e-14, ...
            'stored symmetry is not the half-cell shift; cannot refine it');
        popts.symmetry = half_shift;
    end
    fprintf('%-8s %s, branch point m = %.4f (dE %+.6f) shifted by %+.1e\n', fams{F}, ...
        labels{F}(w), b.m(k), b.dE(k), m_target - b.m(k));
    dE{F} = nan(size(factors{F}));
    for j = 1:numel(factors{F})
        N = factors{F}(j) * G0.N;
        dim = 1 + (size(u0, 2) > 1);             % column vector: 1D, else 2D
        G = SSAV_FCH_helpers.generate_Grid(G0.L, N, dim);
        u = SSAV_FCH_helpers.fourier_resample(u0, N);
        t0 = tic;
        [~, inf_p] = newton_polish_FCH(u, G, Para, popts);
        dE{F}(j) = inf_p.dE;
        fprintf('           x%d  %s: dE = %+.10f  |R| %.1e  converged %d  (%.0f s)\n', ...
            factors{F}(j), mat2str(N), inf_p.dE, inf_p.res_hist(end), inf_p.converged, toc(t0));
    end
end

fprintf('gaps at m = %.4f (negative: the other state is lower)\n', m_target);
for F = 2:3
    if isempty(dE{1}) || isempty(dE{F}), continue, end
    fprintf('  other - %-7s  coarsest %+.3e   finest %+.3e   change %.1e\n', fams{F}, ...
        dE{1}(1) - dE{F}(1), dE{1}(end) - dE{F}(end), abs((dE{1}(end) - dE{F}(end)) - (dE{1}(1) - dE{F}(1))));
end
