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
% u_ref = u0 (0, 1 or 2 of them -- see translation_directions). sigma -> 0
% at a solution; if it doesn't, the state is at a singular point that is
% not a translation (e.g. a symmetry-breaking bifurcation).
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
%     tol           RMS residual tolerance               (1e-10)
%     maxit         max Newton iterations                (30)
%     gmres_maxit   max GMRES iterations per solve       (300)
%     eta_max       max GMRES relative tolerance         (0.1)
%     beta          preconditioner shift                 (1)
%     phase         use translation phase conditions     (true)
%     rank_tol      relative rank threshold for T        (1e-8)
%     verbose       print iteration table                (true)
%
% Outputs
%   u     polished state
%   info  struct: converged, res_hist, gmres_its, step_hist, sigma (columns
%         multiply unit-norm T, so norm(sigma)/sqrt(numel(u)) is on the RMS
%         scale of the residual),
%         lambda, E, T (translation directions), gram_eigs, collapsed,
%         plus function handles mu(u) and Jv(u, v) for reuse (e.g. eigs).

if nargin < 4, opts = struct(); end
opts = set_defaults(opts, struct('tol', 1e-10, 'maxit', 30, 'gmres_maxit', 300, ...
    'eta_max', 0.1, 'beta', 1, 'phase', true, 'rank_tol', 1e-8, 'verbose', true));

eps2 = Para.epsilon^2;
sz = size(u0);
n = numel(u0);
m = mean(u0(:));
rms = @(v) norm(v(:)) / sqrt(numel(v));

lap = @(v) real(ifftn(-Grid.k2 .* fftn(v)));
mu_fun = @(u) compute_mu(u, Para, lap);
Jv_fun = @(u, v) apply_J(u, v, Para, lap);

% right preconditioner: inverse of (eps^2 k^2 + beta)^2, mean mode removed
Pinv = 1 ./ (eps2*Grid.k2 + opts.beta).^2;
Pinv(1) = 0;
apply_Pinv = @(v) real(ifftn(Pinv .* fftn(v)));

u = u0;
u_ref = u0;
if opts.phase
    [T, gram_eigs] = translation_directions(u_ref, Grid, opts.rank_tol);
else
    T = zeros(n, 0);
    gram_eigs = [];
end
p = size(T, 2);

[R, g] = residual(u);
res = rms([R(:); g]);

info.res_hist = res;
info.gmres_its = [];
info.step_hist = [];
info.sigma = zeros(p, 0);

if opts.verbose
    fprintf('Newton polish: n = %d, %d translation direction(s), gram eigs = %s\n', ...
        n, p, mat2str(gram_eigs', 3));
    fprintf('  it   |R|_rms      GMRES  step    |sigma|_rms |u-m|_rms\n');
    fprintf('  %2d   %.3e\n', 0, res);
end

eta_k = opts.eta_max;
for it = 1:opts.maxit
    if res < opts.tol, break, end

    % bordered system, right-preconditioned: unknowns y = [y_u; sigma],
    % du = Pinv*y_u
    Ju = @(v) J_fun(u, v);
    K = @(y) [reshape(proj0(Ju(apply_Pinv(reshape(y(1:n), sz)))), n, 1) + T*y(n+1:end); ...
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

    % Eisenstat-Walker (choice 2) forcing term for the next GMRES solve
    eta_new = 0.9*(res_try/res)^2;
    if 0.9*eta_k^2 > 0.1, eta_new = max(eta_new, 0.9*eta_k^2); end
    eta_k = min(opts.eta_max, eta_new);

    u = u_try;  R = R_try;  g = g_try;  res = res_try;

    info.res_hist(end+1) = res;
    info.gmres_its(end+1) = iters(2);
    info.step_hist(end+1) = step;
    info.sigma(:, end+1) = sigma;

    if opts.verbose
        fprintf('  %2d   %.3e   %4d   %-6.4g  %.3e   %.3e\n', it, res, iters(2), ...
            step, norm(sigma)/sqrt(n), rms(u - m));
    end

    if step < 2^-10
        if opts.verbose, fprintf('  line search failed -- stopping\n'); end
        break
    end
end

info.converged = res < opts.tol;
mu = mu_fun(u);
info.lambda = mean(mu(:));
info.E = SSAV_FCH_helpers.compute_E_direct(u, fftn(u), Para.epsilon, Para.eta1, Grid);
info.T = reshape(T, [sz, p]);
info.gram_eigs = gram_eigs;
% the constant state u = m solves mu = const for every m; flag a collapse
info.collapsed = rms(u - m) < 1e-3;
info.mu = mu_fun;
info.Jv = Jv_fun;

if opts.verbose
    if info.converged
        fprintf('Converged: |R|_rms = %.2e, lambda = %.6f, E = %.10f\n', res, info.lambda, info.E);
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

    function w = J_fun(u, v)
        w = Jv_fun(u, v);
    end

end

function v = proj0(v)
% remove the mean (k = 0 mode)
v = v - mean(v(:));
end

function mu = compute_mu(u, Para, lap)
[~, W1, W2] = SSAV_FCH_helpers.compute_F(u, Para.tau);
eps2 = Para.epsilon^2;
g = eps2*lap(u) - W1;
mu = eps2*lap(g) - W2.*g + Para.eta1*g + (Para.eta1 - Para.eta2)*W1;
end

function w = apply_J(u, v, Para, lap)
% exact Jacobian of compute_mu at u applied to v (W''' = 6u for any tau)
[~, W1, W2] = SSAV_FCH_helpers.compute_F(u, Para.tau);
eps2 = Para.epsilon^2;
g = eps2*lap(u) - W1;
dg = eps2*lap(v) - W2.*v;
w = eps2*lap(dg) - W2.*dg + Para.eta1*dg - 6*u.*v.*g + (Para.eta1 - Para.eta2)*W2.*v;
end

function [T, gram_eigs] = translation_directions(u, Grid, rank_tol)
% Orthonormal basis for span{du/dx, du/dy}, dropping directions in which u
% is (numerically) invariant: stripes have one, a constant state has none.
% Using the eigenvectors of the Gram matrix also handles tilted stripes,
% whose invariant direction is a combination of x and y.
U = fftn(u);
if size(u, 2) == 1
    D = real(ifft(1i*Grid.kxx .* U));
else
    ux = real(ifft2(1i*Grid.kxx .* U));
    uy = real(ifft2(1i*Grid.kyy .* U));
    D = [ux(:), uy(:)];
end
n = size(D, 1);
[V, L] = eig((D'*D) / n);
[gram_eigs, order] = sort(diag(L), 'descend');
V = V(:, order);
% relative threshold, with an absolute floor so a constant state
% (derivatives at roundoff level) gets no directions at all
keep = gram_eigs > max(rank_tol*gram_eigs(1), 1e-12);
T = D * V(:, keep);
T = T ./ vecnorm(T);
end

function opts = set_defaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i})
        opts.(names{i}) = defaults.(names{i});
    end
end
end
