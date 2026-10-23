function branch = continuation_FCH(u0, Grid, Para, opts)
% Pseudo-arclength continuation of FCH steady states in the mass m.
%
% Unknowns (v, m) with u = m + v, mean(v) = 0. Equations
%   F(v, m) = mu(u) - mean(mu(u)) = 0          (steady state at mass m)
%   <v - v0, tv> + (m - m0) tm = ds             (arclength, RMS inner product)
%   T'(u - u0) = 0                              (translation phase conditions)
% where (tv, tm) is the unit tangent at the previous point (v0, m0) and the
% columns of T are the translation directions of u0 (none in even mode).
% Each corrector step solves the bordered system
%   [ J      F_m  T ] [dv   ]     [ F ]
%   [ tv'/n  tm   0 ] [dm   ] = - [ g ]
%   [ T'     0    0 ] [sigma]     [ h ]
% with F_m = dF/dm = J*1 - mean. The bordered matrix stays nonsingular at
% folds in m, so the branch is followed around them. The corrector is plain
% (chord) Newton, not PTC, so unstable stretches are followed too.
%
% Linear solver
%   even mode ('direct'): J is assembled in the even, mean-zero subspace
%     (65^2 - 1 = 4224 unknowns at 128^2, ~1 s) once per accepted point and
%     used for the tangent, the stability count (full eig, ~1 s) and as the
%     chord matrix of the next corrector (re-assembled if contraction is
%     poor). Large domains have many soft interface modes that make
%     Fourier-preconditioned GMRES need hundreds of iterations per solve.
%   full space ('gmres'): matrix-free unrestarted GMRES with the Fourier
%     preconditioner of newton_polish_FCH, stability via stability_FCH.
%
% Inputs
%   u0    a (near-)steady state; m0 = mean(u0)
%   Grid  from SSAV_FCH_helpers.generate_Grid (1D or 2D)
%   Para  epsilon, eta1, eta2, tau
%   opts  (all optional)
%     direction      +1: start towards larger m, -1: smaller     (+1)
%     ds             initial arclength step (RMS units)         (0.01)
%     ds_min, ds_max step-size limits                           (1e-5, 0.2)
%     its_grow       grow ds x1.5 when the corrector needs <= this many
%                    iterations                                 (5)
%     its_shrink     shrink ds x0.7 when it needs >= this many  (9)
%                    (the even-mode chord corrector typically needs 4-7,
%                    so a threshold of 3 kept ds from ever growing back
%                    after a failure halved it)
%     stop_at_fold   stop at the first fold in m (true), or at the
%                    k-th fold (a number k)                     (false)
%     max_steps      max continuation steps                     (200)
%     m_range        stop when m leaves [m_min, m_max]          ([-0.95 0.95])
%     even           stay in the even subspace                  (false)
%     linear_solver  'auto' (direct if even), 'direct', 'gmres' ('auto')
%     tol            corrector RMS residual tolerance           (1e-9)
%     newton_maxit   corrector iterations before ds is halved   (15)
%     gmres_maxit    unrestarted GMRES iterations per solve     (1000)
%     beta           GMRES preconditioner shift                 (1)
%     polish_first   polish u0 with newton_polish_FCH first     (true)
%     stability      compute n_neg at every point               (true)
%     neg_tol        eigenvalues below -neg_tol count negative  (1e-6)
%     stab_k         eigenvalues kept per point                 (4)
%     store_u        keep every profile in branch.u             (true)
%     stop_on_uniform  stop when the pattern dies out (amp < uniform_tol):
%                    the branch has reached the uniform branch u = m,
%                    which solves the equations for every m       (true)
%     uniform_tol    amplitude threshold, RMS of u - m          (1e-4)
%     uniform_cross_tol  a local minimum of the amplitude below this
%                    means the branch passed through the uniform state
%                    between steps (it turns back in m there with the
%                    pattern sign-flipped; not a fold): treated as
%                    reaching the uniform branch                (0.02)
%     symmetry       extra symmetry projection v -> S(v), applied after
%                    the even projection (e.g. averaging with a half-cell
%                    shift to keep two spots identical); S must be linear,
%                    idempotent and commute with J                ([])
%     fold_tm_tol    a fold is counted when dm/ds changes sign relative
%                    to the last value with |dm/ds| > fold_tm_tol, so
%                    jitter on near-vertical stretches is ignored (0.01)
%     verbose        print one line per step                    (true)
%
% Output: branch struct with per-point arrays
%   m, E, dE (E - E_uniform(m)), lambda (chemical potential), amp
%   (RMS of u - m), n_neg, lam (stab_k smallest eigenvalues, one row per
%   point), tm (dm/ds), fold (fold in m between this point and the
%   previous), bif (n_neg changed without a fold: possible branch point),
%   cross_uniform (passed through the uniform state u = m),
%   ds, newton_its, res; u (cell of profiles if store_u); plus Grid, Para,
%   opts, solver, stop_reason.

if nargin < 4, opts = struct(); end
opts = set_defaults(opts, struct('direction', 1, 'ds', 0.01, 'ds_min', 1e-5, ...
    'ds_max', 0.2, 'its_grow', 5, 'its_shrink', 9, 'stop_at_fold', false, 'max_steps', 200, 'm_range', [-0.95 0.95], 'even', false, ...
    'linear_solver', 'auto', 'tol', 1e-9, 'newton_maxit', 15, 'gmres_maxit', 1000, ...
    'beta', 1, 'polish_first', true, 'stability', true, 'neg_tol', 1e-6, 'stab_k', 4, ...
    'store_u', true, 'stop_on_uniform', true, 'uniform_tol', 1e-4, 'uniform_cross_tol', 0.02, ...
    'verbose', true, ...
    'symmetry', [], 'fold_tm_tol', 0.01));
solver = opts.linear_solver;
if strcmp(solver, 'auto')
    if opts.even, solver = 'direct'; else, solver = 'gmres'; end
end
if strcmp(solver, 'direct') && ~opts.even
    error('continuation_FCH: the direct solver works in the even subspace (opts.even = true)');
end

sz = size(u0);
n = numel(u0);
rms = @(a) norm(a(:)) / sqrt(numel(a));
proj0 = @(v) v - mean(v(:));
if opts.even
    sym = @(v) SSAV_FCH_helpers.symmetrize_even(v);
else
    sym = @(v) v;
end
if ~isempty(opts.symmetry)
    sym = @(v) opts.symmetry(sym(v));
end
mu_fun = @(u) SSAV_FCH_helpers.steady_mu(u, Para, Grid);
Jv_fun = @(u, v) SSAV_FCH_helpers.steady_Jv(u, v, Para, Grid);

if strcmp(solver, 'direct')
    Z = SSAV_FCH_helpers.even_basis(sz);    % orthonormal, even, mean zero
    ne = size(Z, 2);
else
    Pinv_sym = 1 ./ (Para.epsilon^2*Grid.k2 + opts.beta).^2;
    Pinv_sym(1) = 0;
    apply_Pinv = @(v) real(ifftn(Pinv_sym .* fftn(v)));
end

% ---- starting point ----
if opts.polish_first
    u0 = newton_polish_FCH(sym(u0), Grid, Para, struct('even', opts.even, ...
        'symmetry', opts.symmetry, 'verbose', false, 'maxit', 3000));
end
u0 = sym(u0);
m = mean(u0(:));
v = u0 - m;
T = phase_dirs(m + v);
lin = linearize(m + v);         % J (and F_m) at the current point

% initial tangent: fix dm = direction through the arclength row
[tv, tm] = tangent(lin, zeros(sz), opts.direction);

branch = init_branch();
record(m + v, m, tm, 0, 0, rms(proj0(mu_fun(m + v))), lin);
if opts.verbose
    fprintf(' step     m          dE              n_neg   ds        its  note\n');
    print_point(1, '');
end

ds = opts.ds;
stop_reason = 'max_steps';
left_start = false;
tm_sign = sign(tm);             % last significant direction in m (fold hysteresis)
k_sig = 1;                      % index of the last point with |tm| > fold_tm_tol
for step = 1:opts.max_steps
    % ---- predictor / corrector, halving ds on failure ----
    while true
        [v_new, m_new, its, ok, res] = corrector(v, m, tv, tm, ds, lin);
        if ok, break, end
        ds = ds / 2;
        if ds < opts.ds_min, break, end
    end
    if ~ok
        stop_reason = 'ds below ds_min (corrector failed)';
        break
    end

    % ---- new point: linearize once, new tangent oriented by the old one ----
    v = v_new;  m = m_new;
    T = phase_dirs(m + v);
    lin = linearize(m + v);
    [tv, tm] = tangent(lin, tv, tm);

    record(m + v, m, tm, ds, its, res, lin);
    k = numel(branch.m);
    if abs(tm) > opts.fold_tm_tol
        if tm_sign ~= 0 && sign(tm) ~= tm_sign
            % turned around since the last significant point: mark the
            % fold where |dm/ds| is smallest in between
            [~, j] = min(abs(branch.tm(k_sig:k)));
            branch.fold(k_sig + j - 1) = true;
        end
        tm_sign = sign(tm);
        k_sig = k;
    end
    branch.bif(k) = ~branch.fold(k) && branch.n_neg(k) ~= branch.n_neg(k-1);
    if opts.verbose
        note = '';
        if any(branch.fold(k_sig:k)) && k_sig == k, note = 'FOLD (marked at the turning point)'; end
        if branch.bif(k), note = 'n_neg changed: possible branch point'; end
        print_point(k, note);
    end

    % ---- step-size control ----
    if its <= opts.its_grow
        ds = min(opts.ds_max, 1.5*ds);
    elseif its >= opts.its_shrink
        ds = max(opts.ds_min, 0.7*ds);
    end


    if m < opts.m_range(1) || m > opts.m_range(2)
        stop_reason = 'left m_range';
        break
    end
    if opts.stop_on_uniform && branch.amp(k) < opts.uniform_tol
        stop_reason = sprintf('reached the uniform branch u = m at m = %.4f', m);
        break
    end
    if k >= 3 && branch.amp(k-1) < opts.uniform_cross_tol && ...
            branch.amp(k-1) < branch.amp(k-2) && branch.amp(k-1) < branch.amp(k)
        % stepped across the uniform state: a turning point in m but not a fold
        branch.fold(k-2:k) = false;
        branch.cross_uniform(k-1) = true;
        if opts.stop_on_uniform
            stop_reason = sprintf('passed through the uniform branch u = m near m = %.4f', branch.m(k-1));
            break
        end
    end
    if opts.stop_at_fold && sum(branch.fold) >= double(opts.stop_at_fold)
        kf_all = find(branch.fold);
        stop_reason = sprintf('fold %d at m = %.4f', numel(kf_all), branch.m(kf_all(end)));
        break
    end
    % closed loop (isola): back within one step of the start after having
    % moved well away from it
    dist_start = sqrt(rms(m + v - branch.u_start)^2 + (m - branch.m(1))^2);
    left_start = left_start || dist_start > 4*ds;
    if left_start && dist_start < ds
        stop_reason = 'returned to the starting point (closed branch)';
        break
    end
end

branch.stop_reason = stop_reason;
branch.solver = solver;
branch.Grid = Grid;
branch.Para = Para;
branch.opts = opts;
branch = rmfield(branch, 'u_start');
if opts.verbose
    fprintf('stopped: %s (%d points, %d fold(s), %d possible branch point(s))\n', ...
        stop_reason, numel(branch.m), sum(branch.fold), sum(branch.bif));
end

% ===================== nested functions =====================

    function T = phase_dirs(u)
        if opts.even
            T = zeros(n, 0);        % translations are odd: not in the even subspace
        else
            T = SSAV_FCH_helpers.translation_directions(u, Grid);
        end
    end

    function lin = linearize(u)
        % everything the linear solves need at u
        lin.u = u;
        Fm = proj0(Jv_fun(u, ones(sz)));
        if strcmp(solver, 'direct')
            JZ = zeros(n, ne);
            for c0 = 1:64:ne                        % J applied to 64 basis vectors at a time
                cols = c0:min(c0 + 63, ne);
                JZ(:, cols) = reshape(Jv_fun(u, reshape(Z(:, cols), [sz, numel(cols)])), n, []);
            end
            Je = Z' * JZ;
            lin.Je = (Je + Je') / 2;
            lin.Fm = Z' * Fm(:);                     % F_m is even and mean zero
        else
            lin.Fm = Fm(:);
        end
    end

    function [dv, dm] = bordered_solve(lin, rv, ra, rp, tv_row, tm_row, lin_tol)
        % solve [J F_m T; tv_row'/n tm_row 0; T' 0 0] [dv; dm; s] = [rv; ra; rp]
        if strcmp(solver, 'direct')
            % in coordinates v = Z c (Z orthonormal; rv is even, mean zero)
            B = [lin.Je, lin.Fm; (Z' * tv_row(:))' / n, tm_row];
            x = B \ [Z' * rv(:); ra];
            dv = reshape(Z * x(1:ne), sz);
            dm = x(ne + 1);
        else
            p = size(T, 2);
            tvr = tv_row(:);
            u = lin.u;
            K = @(y) [reshape(proj0(Jv_fun(u, apply_Pinv(reshape(y(1:n), sz)))), n, 1) + lin.Fm*y(n+1) + T*y(n+2:end); ...
                      tvr' * reshape(apply_Pinv(reshape(y(1:n), sz)), n, 1) / n + tm_row*y(n+1); ...
                      T' * reshape(apply_Pinv(reshape(y(1:n), sz)), n, 1)];
            [y, ~] = gmres(K, [rv(:); ra; rp], [], lin_tol, min(opts.gmres_maxit, n + 1 + p));  % 2 outputs: no console message
            dv = sym(proj0(apply_Pinv(reshape(y(1:n), sz))));
            dm = y(n+1);
        end
    end

    function [tv_new, tm_new] = tangent(lin, tv_row, tm_row)
        % (dv, dm) with J dv + F_m dm = 0, T'dv = 0 and <(dv,dm),(tv_row,tm_row)> = 1,
        % normalized; the last condition orients it along the previous tangent
        [dv, dm] = bordered_solve(lin, zeros(sz), 1, zeros(size(T, 2), 1), tv_row, tm_row, 1e-8);
        nrm = sqrt(rms(dv)^2 + dm^2);
        tv_new = dv / nrm;
        tm_new = dm / nrm;
    end

    function [v1, m1, its, ok, res] = corrector(v0, m0, tv0, tm0, ds, lin0)
        % (chord) Newton on the extended system from the predictor; the
        % direct solver reuses lin0 and re-linearizes if contraction is poor
        v1 = sym(v0 + ds*tv0);
        m1 = m0 + ds*tm0;
        lin_c = lin0;
        ok = false;
        res = inf;
        for it_n = 1:opts.newton_maxit
            its = it_n - 1;                         % corrector updates so far
            u1 = m1 + v1;
            Fv = proj0(mu_fun(u1));
            g = sum((v1 - v0).*tv0, 'all')/n + (m1 - m0)*tm0 - ds;
            h = T' * reshape(v1 - v0, n, 1);
            res_new = rms([Fv(:); g; h / sqrt(n)]);
            if res_new < opts.tol
                ok = true;
                res = res_new;
                return
            end
            if res_new > 10*res && it_n > 2
                return                              % diverging
            end
            if strcmp(solver, 'direct')
                if it_n > 1 && res_new > 0.3*res
                    lin_c = linearize(u1);          % chord contracting poorly: true Newton
                end
            else
                lin_c = linearize(u1);              % cheap for GMRES: only F_m
            end
            res = res_new;
            lin_tol = min(1e-3, max(1e-8, 0.1*opts.tol/res));
            [dv, dm] = bordered_solve(lin_c, -Fv, -g, -h, tv0, tm0, lin_tol);
            v1 = sym(v1 + dv);
            m1 = m1 + dm;
        end
        its = opts.newton_maxit;
        res = rms(proj0(mu_fun(m1 + v1)));
        ok = res < opts.tol;
    end

    function b = init_branch()
        b = struct('m', [], 'E', [], 'dE', [], 'lambda', [], 'amp', [], 'n_neg', [], ...
            'lam', zeros(0, opts.stab_k), 'tm', [], 'fold', false(0), 'bif', false(0), ...
            'cross_uniform', false(0), ...
            'ds', [], 'newton_its', [], 'res', [], 'u', {{}}, 'u_start', u0);
    end

    function record(u, m, tm, ds, its, res, lin)
        E = SSAV_FCH_helpers.steady_energy(u, Para, Grid);
        mu = mu_fun(u);
        branch.m(end+1) = m;
        branch.E(end+1) = E;
        branch.dE(end+1) = E - SSAV_FCH_helpers.steady_energy(m*ones(sz), Para, Grid);
        branch.lambda(end+1) = mean(mu(:));
        branch.amp(end+1) = rms(u - m);
        branch.tm(end+1) = tm;
        branch.ds(end+1) = ds;
        branch.newton_its(end+1) = its;
        branch.res(end+1) = res;
        branch.fold(end+1) = false;
        branch.bif(end+1) = false;
        branch.cross_uniform(end+1) = false;
        lam = nan(1, opts.stab_k);
        n_neg = NaN;
        if opts.stability
            if strcmp(solver, 'direct')
                ev = sort(eig(lin.Je));             % exact, ~1 s at 128^2 even
                n_neg = sum(ev < -opts.neg_tol);
                kk = min(opts.stab_k, numel(ev));
                lam(1:kk) = ev(1:kk);
            else
                [ls, ~, sinfo] = stability_FCH(u, Grid, Para, struct('even_only', opts.even, ...
                    'k', opts.stab_k, 'neg_tol', opts.neg_tol, 'verbose', false));
                n_neg = sinfo.n_neg;
                lam(1:numel(ls)) = ls;
            end
        end
        branch.n_neg(end+1) = n_neg;
        branch.lam(end+1, :) = lam;
        if opts.store_u
            branch.u{end+1} = u;
        end
    end

    function print_point(k, note)
        fprintf(' %4d  %+.5f  %+.10f   %3g   %.2e  %2d   %s\n', k - 1, branch.m(k), ...
            branch.dE(k), branch.n_neg(k), branch.ds(k), branch.newton_its(k), note);
    end
end

function opts = set_defaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i})
        opts.(names{i}) = defaults.(names{i});
    end
end
end
