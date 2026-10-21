function [u, info] = newton_polish_FCH(u0, Grid, Para, opts)
% Polish an approximate FCH steady state with Newton-GMRES.
%
% Solves  mu(u) = lambda,  mean(u) = mean(u0)  for u, where
%   E(u)  = int 1/2 (eps^2 Lap u - W'(u))^2 - eta1 eps^2/2 |grad u|^2 - eta2 W(u)
%   mu(u) = (eps^2 Lap - W''(u) + eta1)(eps^2 Lap u - W'(u)) + (eta1 - eta2) W'(u)
% is the energy used by FCH_BDF2_SAV_adap_dt (compute_E + compute_G) and
% W = compute_F(u, Para.tau). lambda = mean(mu) is the Lagrange multiplier
% for the mass constraint, so the residual is R(u) = mu(u) - mean(mu(u)).
%
% Translation invariance: shifted copies of a steady state are steady
% states, so grad u spans (nearly) null directions of the Jacobian. Each
% Newton step solves the bordered system
%   [ J   T ] [du   ]   [ -R ]
%   [ T'  0 ] [sigma] = [ -g ],   g = T'(u - u_ref),
% where the columns of T are the independent translation directions of
% u_ref (0, 1 or 2 of them -- see translation_directions). u_ref is reset
% to the current iterate after every accepted step, so T always matches the
% current pattern and g = 0 at the start of each step. sigma -> 0
% at a solution; if it doesn't, the state is at a singular point that is
% not a translation (e.g. a symmetry-breaking bifurcation).
%
% Globalization (opts.ptc, default on): pseudo-transient continuation.
% Each step solves (J + (1/dtau) M) du = -R with M = (-Lap)^{-1} on
% mean-zero functions, i.e. a linearized backward-Euler step of the PDE
% u_t = Lap mu with pseudo-time step dtau. Steps are accepted when the
% energy decreases (or the residual halves); dtau grows on acceptance and
% shrinks on rejection, so far from a solution this follows the gradient
% flow with large implicit steps and near one it becomes Newton. Because it
% follows the flow it converges to stable states; it will move away from a
% saddle rather than converge to it (use ptc = false, or continuation).
%
% The linear solve is matrix-free GMRES with an exact Jacobian-vector
% product and a Fourier-diagonal right preconditioner (eps^2 k^2 + beta)^2,
% which matches J at high wavenumber. Works for 1D (column) or 2D input.
%
% Inputs
%   u0    initial guess (N x 1 or Ny x Nx, same layout as the time-stepper)
%   Grid  from SSAV_FCH_helpers.generate_Grid
%   Para  needs epsilon, eta1, eta2, tau
%   opts  (all optional)
%     tol           RMS residual tolerance               (1e-9)
%     maxit         max Newton iterations                (100)
%     gmres_maxit   max GMRES iterations per solve       (60)
%                   (capping inexact early PTC solves: 2-5x faster at
%                   128^2 than 300, same converged states)
%     eta_max       max GMRES relative tolerance         (0.1)
%     beta          preconditioner shift                 (1)
%     phase         use translation phase conditions     (true)
%     rank_tol      relative rank threshold for T        (1e-8)
%     verbose       print iteration table                (true)
%     ptc           pseudo-transient continuation        (true)
%                   false = plain Newton + residual line search
%     dtau0         initial pseudo-time step             (1)
%     dtau_min      give up below this pseudo-time step  (1e-8)
%     dtau_max      cap on the pseudo-time step          (1e12)
%     dtau_grow     max growth of dtau per accepted step (3)
%     stall_tol     accept a stalled residual below this     (1e-8)
%     stall_window  iterations without halving = stalled     (50)
%     even          keep u even in x and y (symmetrize u0 and
%                   every accepted iterate)                  (false)
%     symmetry      extra symmetry projection u -> S(u) applied with
%                   the even one (e.g. a half-cell shift average) ([])
%
% Outputs
%   u     polished state
%   info  struct: converged (res < tol, or stalled below stall_tol),
%         stalled, res_hist, E_hist, gmres_its, step_hist,
%         dtau_hist, accepted, sigma_rms (phase-condition multipliers on the
%         RMS scale of the residual; -> 0 at a solution),
%         lambda, E (energy consistent with mu), E_direct (compute_E_direct),
%         E_uniform (energy of u = m), dE = E - E_uniform,
%         T (translation directions of the final state), gram_eigs, collapsed,
%         plus function handles mu(u) and Jv(u, v) for reuse (e.g. eigs).

if nargin < 4, opts = struct(); end
opts = set_defaults(opts, struct('tol', 1e-9, 'maxit', 100, 'gmres_maxit', 60, ...
    'eta_max', 0.1, 'beta', 1, 'phase', true, 'rank_tol', 1e-8, 'verbose', true, ...
    'ptc', true, 'dtau0', 1, 'dtau_min', 1e-8, 'dtau_max', 1e12, 'dtau_grow', 3, ...
    'stall_tol', 1e-8, 'stall_window', 50, 'even', false, 'symmetry', []));

eps2 = Para.epsilon^2;
sz = size(u0);
n = numel(u0);
m = mean(u0(:));
rms = @(v) norm(v(:)) / sqrt(numel(v));

% even mode: keep u even in x and y. A Newton step from an even state is
% even in exact arithmetic, but PTC follows the flow and amplifies roundoff
% asymmetry near saddles whose unstable modes are odd (~1e-3 in the Ex5
% sweep), so re-symmetrize every accepted iterate.
if opts.even
    [u0, asym0] = SSAV_FCH_helpers.symmetrize_even(u0);
    if opts.verbose && asym0 > 1e-6
        fprintf('even mode: removed relative asymmetry %.1e from u0\n', asym0);
    end
end
if ~isempty(opts.symmetry)
    u0 = opts.symmetry(u0);
end

mu_fun = @(u) SSAV_FCH_helpers.steady_mu(u, Para, Grid);
Jv_fun = @(u, v) SSAV_FCH_helpers.steady_Jv(u, v, Para, Grid);

% (-Lap)^{-1} on mean-zero functions: the H^{-1} metric of u_t = Lap mu
Minv_sym = 1 ./ Grid.k2;
Minv_sym(1) = 0;
apply_Minv = @(v) real(ifftn(Minv_sym .* fftn(v)));

u = u0;
u_ref = u0;
if opts.phase
    [T, gram_eigs] = SSAV_FCH_helpers.translation_directions(u_ref, Grid, opts.rank_tol);
else
    T = zeros(n, 0);
    gram_eigs = [];
end
p = size(T, 2);

[R, g] = residual(u);
res = rms([R(:); g]);
E = SSAV_FCH_helpers.steady_energy(u, Para, Grid);
dtau = opts.dtau0;

info.res_hist = res;
info.E_hist = E;
info.gmres_its = [];
info.step_hist = [];
info.dtau_hist = [];
info.accepted = logical([]);
info.sigma_rms = [];

if opts.verbose
    fprintf('Newton polish (%s): n = %d, %d translation direction(s), gram eigs = %s\n', ...
        ternary(opts.ptc, 'PTC', 'line search'), n, p, mat2str(gram_eigs', 3));
    if opts.ptc
        fprintf('  it   |R|_rms      GMRES  dtau       acc  |sigma|_rms  E\n');
    else
        fprintf('  it   |R|_rms      GMRES  step    |sigma|_rms |u-m|_rms\n');
    end
    fprintf('  %2d   %.3e\n', 0, res);
end

eta_k = opts.eta_max;
stalled = false;
for it = 1:opts.maxit
    if res < opts.tol, break, end

    % PTC adds (1/dtau)(-Lap)^{-1}: a linearized backward-Euler step of
    % u_t = Lap mu with step dtau. dtau -> inf recovers Newton.
    if opts.ptc
        shift = 1/dtau;
    else
        shift = 0;
    end

    % right preconditioner: inverse of (eps^2 k^2 + beta)^2 + shift/k^2,
    % which matches J + shift*Minv at high and (for PTC) low wavenumber
    Pinv = 1 ./ ((eps2*Grid.k2 + opts.beta).^2 + shift*Minv_sym);
    Pinv(1) = 0;
    apply_Pinv = @(v) real(ifftn(Pinv .* fftn(v)));

    % bordered system, right-preconditioned: unknowns y = [y_u; sigma],
    % du = Pinv*y_u
    A = @(v) proj0(Jv_fun(u, v) + shift*apply_Minv(v));
    K = @(y) [reshape(A(apply_Pinv(reshape(y(1:n), sz))), n, 1) + T*y(n+1:end); ...
              T' * reshape(apply_Pinv(reshape(y(1:n), sz)), n, 1)];
    rhs = -[R(:); g];

    % don't solve the linear system more accurately than the target needs
    eta_use = max(eta_k, 0.5*opts.tol/res);
    [y, flag, relres, iters] = gmres(K, rhs, [], eta_use, min(opts.gmres_maxit, n + p));
    if flag ~= 0 && opts.verbose
        fprintf('  GMRES flag %d (relres %.2e after %d its)\n', flag, relres, iters(2));
    end
    du = apply_Pinv(reshape(y(1:n), sz));
    sigma = y(n+1:end);

    if opts.ptc
        % Accept if the energy goes down (we are following the gradient
        % flow, so the residual may rise, e.g. while leaving a saddle) or the
        % residual drops a lot (near convergence the energy change is below
        % roundoff). Grow dtau on acceptance, shrink it on rejection.
        step = 1;
        u_try = u + du;
        [R_try, g_try] = residual(u_try);
        res_try = rms([R_try(:); g_try]);
        E_try = SSAV_FCH_helpers.steady_energy(u_try, Para, Grid);
        accept = E_try <= E + 1e-14 || res_try < 0.5*res;
        dtau_used = dtau;
        if accept
            dtau = min(opts.dtau_max, dtau * min(opts.dtau_grow, max(1.5, res/res_try)));
        else
            dtau = dtau / 4;
        end
    else
        % backtracking line search on the residual norm
        step = 1;
        while true
            u_try = u + step*du;
            [R_try, g_try] = residual(u_try);
            res_try = rms([R_try(:); g_try]);
            if res_try <= (1 - 1e-4*step)*res || step < 2^-10
                break
            end
            step = step / 2;
        end
        accept = true;
        dtau_used = inf;
        E_try = SSAV_FCH_helpers.steady_energy(u_try, Para, Grid);
    end

    if accept
        % Eisenstat-Walker (choice 2) forcing term for the next GMRES solve
        eta_new = 0.9*(res_try/res)^2;
        if 0.9*eta_k^2 > 0.1, eta_new = max(eta_new, 0.9*eta_k^2); end
        eta_k = min(opts.eta_max, eta_new);

        u = u_try;  R = R_try;  g = g_try;  res = res_try;  E = E_try;
        if opts.even || ~isempty(opts.symmetry)
            if opts.even, u = SSAV_FCH_helpers.symmetrize_even(u); end
            if ~isempty(opts.symmetry), u = opts.symmetry(u); end
            [R, g] = residual(u);
            res = rms([R(:); g]);
            E = SSAV_FCH_helpers.steady_energy(u, Para, Grid);
        end

        % Re-anchor the phase conditions at the current state. PTC can
        % reshape the pattern so much that the starting state's translation
        % directions end up nearly orthogonal to the current ones (87 deg
        % in one m = 0.5 run), which leaves the bordered system nearly
        % singular and stalls Newton.
        if opts.phase
            u_ref = u;
            [T, gram_eigs] = SSAV_FCH_helpers.translation_directions(u_ref, Grid, opts.rank_tol);
            p = size(T, 2);
            [R, g] = residual(u);
        end
    end

    info.res_hist(end+1) = res;
    info.E_hist(end+1) = E;
    info.gmres_its(end+1) = iters(2);
    info.step_hist(end+1) = step;
    info.dtau_hist(end+1) = dtau_used;
    info.accepted(end+1) = accept;
    info.sigma_rms(end+1) = norm(sigma)/sqrt(n);

    if opts.verbose
        if opts.ptc
            fprintf('  %2d   %.3e   %4d   %.2e   %d    %.3e    %+.10f\n', it, res, iters(2), ...
                dtau_used, accept, norm(sigma)/sqrt(n), E);
        else
            fprintf('  %2d   %.3e   %4d   %-6.4g  %.3e   %.3e\n', it, res, iters(2), ...
                step, norm(sigma)/sqrt(n), rms(u - m));
        end
    end

    if ~opts.ptc && step < 2^-10
        if opts.verbose, fprintf('  line search failed -- stopping\n'); end
        break
    end
    if opts.ptc && dtau < opts.dtau_min
        if opts.verbose, fprintf('  dtau fell below dtau_min -- stopping\n'); end
        break
    end

    % Roundoff floor: some states can't reach tol (one m = 0.7 state sits at
    % ~2e-9 for thousands of steps). Stop once the residual is below
    % stall_tol and hasn't halved over the last stall_window iterations.
    w = opts.stall_window;
    if it > w && res < opts.stall_tol && ...
            min(info.res_hist(end-w+1:end)) > 0.5*min(info.res_hist(1:end-w))
        stalled = true;
        if opts.verbose
            fprintf('  residual stalled at %.1e (roundoff floor) -- stopping\n', res);
        end
        break
    end
end

info.converged = res < opts.tol || stalled;
info.stalled = stalled;
mu = mu_fun(u);
info.lambda = mean(mu(:));
% E: the energy the solver uses (u is an exact critical point of it).
% E_direct: compute_E_direct, which zeroes the Nyquist mode of its
% derivatives, so it differs at the ~1e-9 level.
info.E = SSAV_FCH_helpers.steady_energy(u, Para, Grid);
info.E_direct = SSAV_FCH_helpers.compute_E_direct(u, fftn(u), Para.epsilon, Para.eta1, Grid);
% energy relative to the uniform state u = m (same energy functional)
info.E_uniform = SSAV_FCH_helpers.steady_energy(m*ones(sz), Para, Grid);
info.dE = info.E - info.E_uniform;
info.T = reshape(T, [sz, p]);
info.gram_eigs = gram_eigs;
% the constant state u = m solves mu = const for every m; flag a collapse
info.collapsed = rms(u - m) < 1e-3;
info.mu = mu_fun;
info.Jv = Jv_fun;

if opts.verbose
    if info.converged
        fprintf('Converged: |R|_rms = %.2e, lambda = %.6f, E = %.10f, dE = %.10f\n', ...
            res, info.lambda, info.E, info.dE);
    else
        fprintf('NOT converged: |R|_rms = %.2e after %d iterations\n', res, numel(info.res_hist) - 1);
    end
    if info.collapsed
        fprintf('WARNING: collapsed to the constant state u = %.3f\n', m);
    end
end

    function [R, g] = residual(u)
        R = proj0(mu_fun(u));
        g = T' * (u(:) - u_ref(:));
    end

end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function v = proj0(v)
% remove the mean (k = 0 mode)
v = v - mean(v(:));
end

function opts = set_defaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i})
        opts.(names{i}) = defaults.(names{i});
    end
end
end
