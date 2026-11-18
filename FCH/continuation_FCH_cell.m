function br = continuation_FCH_cell(u0, d0, Para, opts)
% Pseudo-arclength continuation in the mass m of FCH steady states on one
% periodic cell of size d, with the cell's symmetry enforced, and (switch
% opts.optimize_d) optionally with d as an extra unknown fixed by
% dE/dd = 0: the branch of locally optimal cell sizes.
%
% Cells ('cell' option):
%   'stripe'  1D, period d; u a column (stripes, u independent of y)
%   'hex'     one exact hex cell, d x d*sqrt(3); u Ny x Nx, restricted to
%             the full hexagonal (D6) symmetry (SSAV_FCH_helpers.hex_symmetry)
% States are even; E(-u) = E(u), so m < 0 mirrors m > 0.
%
% Unknowns: the profile v = u - m in an orthonormal basis Z of the even
% (and for 'hex' D6-symmetric) mean-zero subspace, v = Z a; the mass m;
% and, with optimize_d, the cell size d. Equations:
%   Z'(mu - mean(mu)) = 0                         (steady state)
%   G(u, d) = 0                                   (optimize_d only)
%   <a - a0, ta>/n + (m - m0) tm + w^2 (d - d0) td = ds   (arclength)
% with w = opts.d_weight. Chord Newton with a direct bordered solve.
%
% Optimal cell size: with the grid values fixed, only the wavenumbers
% depend on d (k^2 ~ 1/d^2). With g = eps^2 Lap u - W'(u) and E the
% discrete energy per unit area (steady_energy),
%   dE/dd = -G/d,   G = mean(2 eps^2 g Lap u + eta1 eps^2 u Lap u),
% and at a steady state this is also the derivative of E along the family
% of steady states at fixed m. So G = 0 is the extra equation (zero stress
% on the cell). The cell is scaled as a whole (one size d), which is
% enough for stripes and hex states: their symmetry makes the stress
% isotropic. Other patterns need separate cell lengths.
%
% The number of grid points is fixed, so the grid spacing changes with d;
% a warning is printed if it exceeds h_max.
%
% Stability: n_neg = negative eigenvalues of the Hessian at fixed d on the
% subspace; with optimize_d also n_neg_d = negative eigenvalues of the
% Hessian in (a, d) together (0: a local minimum also with respect to the
% cell size; n_neg + 1: d is a local maximum of E). For 'hex' both count
% D6-symmetric perturbations only.
%
% Inputs
%   u0   starting state (near a steady state) on the cell of size d0
%   d0   cell size (with optimize_d: the starting guess, e.g. the best d
%        of a spacing scan)
%   Para epsilon, eta1, eta2, tau
%   opts (all optional)
%     cell          'stripe' or 'hex'                          ('stripe')
%     optimize_d    d is an unknown, fixed by dE/dd = 0        (false)
%     d_weight      weight w of d in the arclength norm        (0.1)
%     direction     +1: start towards larger m, -1: smaller    (+1)
%     ds, ds_min, ds_max  arclength step and limits           (0.01, 1e-5, 0.05)
%     max_steps     continuation steps                         (500)
%     m_range       stop when m leaves it                      ([0 0.99])
%     d_range       stop when d leaves it                      ([0.5 20])
%     amp_min       optimize_d: stop when the pattern amplitude (RMS of
%                   u - m) drops below this; the optimal d is not defined
%                   at the uniform state                       (0.02)
%     tol           residual tolerance                         (1e-9)
%     newton_maxit  corrector iterations before ds is halved   (12)
%     h_max         warn when the grid spacing exceeds this    (0.12)
%     polish_first  first solve at the starting m (with optimize_d: for
%                   u and d)                                   (true)
%     verbose       print one line per step                    (true)
%
% Output: struct with per-point arrays m, d, E, amp, n_neg, n_neg_d (NaN
% unless optimize_d), G (the stress; ~0 with optimize_d), res, its, ds,
% fold (fold in m), u (cell array of profiles); plus N, stop_reason, opts.

if nargin < 4, opts = struct(); end
opts = set_defaults(opts, struct('cell', 'stripe', 'optimize_d', false, 'd_weight', 0.1, ...
    'direction', 1, 'ds', 0.01, 'ds_min', 1e-5, 'ds_max', 0.05, 'max_steps', 500, ...
    'm_range', [0 0.99], 'd_range', [0.5 20], 'amp_min', 0.02, 'tol', 1e-9, ...
    'newton_maxit', 12, 'h_max', 0.12, 'polish_first', true, 'verbose', true));

eps2 = Para.epsilon^2;
fd = opts.optimize_d;           % d is an unknown
w2 = opts.d_weight^2;
switch opts.cell
    case 'stripe', Lhat = 1;           dim = 1;
    case 'hex',    Lhat = [1 sqrt(3)]; dim = 2;
    otherwise, error('continuation_FCH_cell: cell must be ''stripe'' or ''hex''');
end
proj0 = @(w) w - mean(w(:));
sz = size(u0);
if dim == 1
    Nsz = numel(u0);
else
    Nsz = [sz(2), sz(1)];
end
n = numel(u0);

% ---- basis of the symmetric subspace (independent of d) ----
Z = SSAV_FCH_helpers.even_basis(sz);
if dim == 2
    S = SSAV_FCH_helpers.hex_symmetry(grid_at(d0));
    ZS = zeros(size(Z));
    for j = 1:size(Z, 2)
        ZS(:, j) = reshape(S(reshape(Z(:, j), sz)), [], 1);
    end
    [U, sv] = svd(ZS, 'econ');
    sv = diag(sv);
    Z = U(:, sv > 1e-8 * max(sv));
    clear ZS U
end

% ---- starting point ----
m = mean(u0(:));
d = d0;
a = Z' * (u0(:) - m);
start_ok = true;
if opts.polish_first
    [a, d, ok_p, res_p] = polish(a, m, d);
    if ~ok_p
        warning('continuation_FCH_cell: starting polish did not converge (residual %.2e); not continuing', res_p);
        start_ok = false;
    end
end
u = to_u(a, m);
lin = linearize(u, d);
[ta, tm, td] = tangent(lin, [], [], [], true);
if sign(tm) ~= opts.direction, ta = -ta; tm = -tm; td = -td; end

br = struct('m', [], 'd', [], 'E', [], 'amp', [], 'n_neg', [], 'n_neg_d', [], 'G', [], ...
    'res', [], 'its', [], 'ds', [], 'fold', false(0), 'u', {{}});
record(u, m, d, lin, 0, 0, residual_norm(u, d));
if opts.verbose
    fprintf(' step     m        d        E             n_neg n_neg_d  ds        its\n');
    print_point(1);
end
warned = false;

ds = opts.ds;
stop_reason = 'max_steps';
if ~start_ok, stop_reason = 'starting polish did not converge'; end
tm_prev = tm;
for step = 1:start_ok*opts.max_steps
    while true
        [a1, m1, d1, its, ok, res] = corrector(a, m, d, ta, tm, td, ds, lin);
        if ok, break, end
        ds = ds / 2;
        if ds < opts.ds_min, break, end
    end
    if ~ok
        stop_reason = 'ds below ds_min (corrector failed)';
        break
    end
    a = a1;  m = m1;  d = d1;
    u = to_u(a, m);
    lin = linearize(u, d);
    [ta, tm, td] = tangent(lin, ta, tm, td, false);
    record(u, m, d, lin, ds, its, res);
    k = numel(br.m);
    br.fold(k) = sign(tm) ~= sign(tm_prev) && tm_prev ~= 0;
    tm_prev = tm;
    if opts.verbose, print_point(k); end
    if ~warned && d / Nsz(1) > opts.h_max
        warning('continuation_FCH_cell: grid spacing %.3f > h_max = %.3f at d = %.3f (fixed %s grid)', ...
            d / Nsz(1), opts.h_max, d, mat2str(Nsz));
        warned = true;
    end

    if its <= 4
        ds = min(opts.ds_max, 1.5*ds);
    elseif its >= 8
        ds = max(opts.ds_min, 0.7*ds);
    end

    if m < opts.m_range(1) || m > opts.m_range(2)
        stop_reason = 'left m_range';
        break
    end
    if d < opts.d_range(1) || d > opts.d_range(2)
        stop_reason = 'left d_range';
        break
    end
    if fd && br.amp(k) < opts.amp_min
        stop_reason = 'pattern amplitude below amp_min (near the uniform state)';
        break
    end
end
br.N = Nsz;
br.stop_reason = stop_reason;
br.opts = opts;
if opts.verbose
    fprintf('stopped: %s (%d points, %d fold(s) in m)\n', stop_reason, numel(br.m), sum(br.fold));
end

% ===================== nested functions =====================

    function Gr = grid_at(d)
        Gr = SSAV_FCH_helpers.generate_Grid(d*Lhat, Nsz, dim);
    end

    function u = to_u(a, m)
        u = m + reshape(Z * a, sz);
    end

    function [F, G, gradG, lu] = evaluate(u, d)
        % projected residual Z'F, the stress G and its gradient in the mean
        % inner product (dG = mean(gradG .* du))
        Gr = grid_at(d);
        lap = @(w) SSAV_FCH_helpers.spectral_lap(w, Gr);
        [~, W1, W2] = SSAV_FCH_helpers.compute_F(u, Para.tau);
        lu = lap(u);
        g = eps2*lu - W1;
        mu = eps2*lap(g) - W2.*g + Para.eta1*g + (Para.eta1 - Para.eta2)*W1;
        F = Z' * reshape(proj0(mu), [], 1);
        G = mean(2*eps2*g.*lu + Para.eta1*eps2*u.*lu, 'all');
        gradG = 2*eps2*(eps2*lap(lu) - W2.*lu + lap(g)) + 2*Para.eta1*eps2*lu;
    end

    function r = residual_norm(u, d)
        [F, G] = evaluate(u, d);
        r = sqrt(sum(F.^2)/n + fd*G^2);
    end

    function lin = linearize(u, d)
        % blocks of the Jacobian of (Z'F, G) with respect to (a, m, d)
        Gr = grid_at(d);
        [~, G, gradG, lu] = evaluate(u, d);
        ne = size(Z, 2);
        JZ = zeros(n, ne);
        for c0 = 1:64:ne
            cols = c0:min(c0 + 63, ne);
            JZ(:, cols) = reshape(SSAV_FCH_helpers.steady_Jv(u, reshape(Z(:, cols), ...
                [sz, numel(cols)]), Para, Gr), n, []);
        end
        Je = Z' * JZ;
        lin.Je = (Je + Je') / 2;
        lin.Fm = Z' * reshape(proj0(SSAV_FCH_helpers.steady_Jv(u, ones(sz), Para, Gr)), [], 1);
        lin.Fd = -(Z' * gradG(:)) / d;              % d(Z'F)/dd
        lin.Gv = (Z' * gradG(:))' / n;              % dG/da
        lin.Gm = mean(gradG(:));                    % dG/dm
        lin.Gd = -(2/d) * (2*eps2^2*mean(lu(:).^2) + G);
        lin.G = G;
    end

    function A = jac(lin)
        % Jacobian of the equations (without the arclength row) in the unknowns
        if fd
            A = [lin.Je, lin.Fm, lin.Fd; lin.Gv, lin.Gm, lin.Gd];
        else
            A = [lin.Je, lin.Fm];
        end
    end

    function r = arc_row(ta, tm, td)
        if fd
            r = [ta' / n, tm, w2*td];
        else
            r = [ta' / n, tm];
        end
    end

    function [a, d, ok, res] = polish(a, m, d)
        % Newton at fixed m, with backtracking. With optimize_d:
        %   1. solve for a at the given d;
        %   2. walk d downhill in energy (towards G = 0) through steady
        %      states at fixed d, each warm-started from the last, until G
        %      changes sign, and take the secant estimate of the root;
        %   3. Newton for a and d together, each step in d limited to 10%
        %      of d.
        % (Newton for both from a d far from the optimum can fail or run
        % off to the uniform state, where G = 0 for every d.)
        if fd
            fd = false;
            [a, ~, ok0, res] = polish(a, m, d);
            if ~ok0
                fd = true;
                ok = false;
                return
            end
            [~, G] = evaluate(to_u(a, m), d);
            s = sign(G);                % dE/dd = -G/d: energy decreases towards sign(G)
            hstep = 0.02*d;
            for walk = 1:300*(s ~= 0)
                d_try = d + s*hstep;
                [a_try, ~, ok1] = polish(a, m, d_try);
                if ~ok1
                    hstep = hstep / 2;
                    if hstep < 1e-4*d, break, end
                    continue
                end
                [~, G_try] = evaluate(to_u(a_try, m), d_try);
                if sign(G_try) ~= s
                    t = G / (G - G_try);                % secant: G = 0 at d + t*(d_try - d)
                    a = a + t*(a_try - a);
                    d = d + t*(d_try - d);
                    break
                end
                a = a_try;  d = d_try;  G = G_try;
                hstep = min(1.5*hstep, 0.1*d);
            end
            fd = true;
        end
        ok = false;
        for it = 1:50
            uu = to_u(a, m);
            [F, G] = evaluate(uu, d);
            res = sqrt(sum(F.^2)/n + fd*G^2);
            if res < opts.tol, ok = true; return, end
            ln = linearize(uu, d);
            if fd
                x = -[ln.Je, ln.Fd; ln.Gv, ln.Gd] \ [F; G];
                x = x * min(1, 0.1*d / max(abs(x(end)), realmin));
            else
                x = [-(ln.Je \ F); 0];
            end
            s = 1;
            while s > 2^-10
                at = a + s*x(1:end-1);
                dt = d + s*x(end);
                if residual_norm(to_u(at, m), dt) < res, break, end
                s = s / 2;
            end
            a = at;  d = dt;
        end
    end

    function [ta, tm, td] = tangent(lin, ta0, tm0, td0, initial)
        % null direction of the Jacobian, oriented by the old tangent
        % (initial: by dm = 1)
        A = jac(lin);
        if initial
            row = zeros(1, size(A, 2));
            row(size(Z, 2) + 1) = 1;
        else
            row = arc_row(ta0, tm0, td0);
        end
        x = [A; row] \ [zeros(size(A, 1), 1); 1];
        ne = size(Z, 2);
        ta = x(1:ne);
        tm = x(ne + 1);
        if fd, td = x(ne + 2); else, td = 0; end
        nrm = sqrt(sum(ta.^2)/n + tm^2 + w2*td^2);
        ta = ta / nrm;  tm = tm / nrm;  td = td / nrm;
    end

    function [a1, m1, d1, its, ok, res] = corrector(a0, m0, d0, ta, tm, td, ds, lin0)
        % chord Newton from the predictor; re-linearizes if contraction is poor
        a1 = a0 + ds*ta;  m1 = m0 + ds*tm;  d1 = d0 + ds*td;
        ln = lin0;
        B = [];
        ok = false;
        res = inf;
        ne = size(Z, 2);
        for it = 1:opts.newton_maxit
            its = it - 1;
            uu = to_u(a1, m1);
            [F, G] = evaluate(uu, d1);
            arc = (a1 - a0)' * ta / n + (m1 - m0)*tm + w2*(d1 - d0)*td - ds;
            res_new = sqrt(sum(F.^2)/n + fd*G^2 + arc^2);
            if res_new < opts.tol
                ok = true;
                res = res_new;
                return
            end
            if res_new > 10*res && it > 2
                return                              % diverging
            end
            if isempty(B) || (it > 1 && res_new > 0.3*res)
                if it > 1, ln = linearize(uu, d1); end
                B = decomposition([jac(ln); arc_row(ta, tm, td)]);
            end
            res = res_new;
            if fd
                x = -(B \ [F; G; arc]);
                d1 = d1 + x(ne + 2);
            else
                x = -(B \ [F; arc]);
            end
            a1 = a1 + x(1:ne);
            m1 = m1 + x(ne + 1);
        end
        its = opts.newton_maxit;
    end

    function record(u, m, d, lin, ds, its, res)
        br.m(end+1) = m;
        br.d(end+1) = d;
        br.E(end+1) = SSAV_FCH_helpers.steady_energy(u, Para, grid_at(d));
        br.amp(end+1) = sqrt(mean((u(:) - m).^2));
        br.n_neg(end+1) = sum(eig(lin.Je) < -1e-6);
        if fd
            % Hessian of E in (a, d), times n: [Je, Fd; Fd', n*(G/d^2 - Gd/d)]
            H = [lin.Je, lin.Fd; lin.Fd', n*(lin.G/d^2 - lin.Gd/d)];
            br.n_neg_d(end+1) = sum(eig((H + H')/2) < -1e-6);
        else
            br.n_neg_d(end+1) = NaN;
        end
        br.G(end+1) = lin.G;
        br.res(end+1) = res;
        br.its(end+1) = its;
        br.ds(end+1) = ds;
        br.fold(end+1) = false;
        br.u{end+1} = u;
    end

    function print_point(k)
        fprintf(' %4d  %+.5f  %.4f  %+.9f   %3d   %3g    %.2e  %2d\n', k - 1, br.m(k), ...
            br.d(k), br.E(k), br.n_neg(k), br.n_neg_d(k), br.ds(k), br.its(k));
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
