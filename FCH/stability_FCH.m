function [lam, modes, info] = stability_FCH(u, Grid, Para, opts)
% Linear stability of an FCH steady state (e.g. from newton_polish_FCH).
%
% The state is a local energy minimum at fixed mass iff the Jacobian J of
% mu (= Hessian of E) is positive on mean-zero perturbations orthogonal to
% the translation directions (shifting the pattern costs no energy, so those
% eigenvalues are 0 and are excluded, as is the mean/mass mode).
%
% Methods
%   'dense'  assemble J and diagonalize on the constrained subspace; exact
%            eigenvalues. Default for numel(u) <= 5000 (64^2: ~20 s).
%   'eigs'   Lanczos on P^{-1/2} J P^{-1/2}, P = (eps^2 k^2 + beta)^2 (the
%            Newton preconditioner). By Sylvester's law of inertia it has the
%            same number of negative/zero eigenvalues as J, so the stability
%            verdict is exact, but the eigenvalues are rescaled; info.rayleigh
%            gives the true curvature d'Jd/d'd of each returned mode.
%
% Inputs
%   u, Grid, Para  as for newton_polish_FCH (u should be polished)
%   opts (optional)
%     k        number of smallest eigenvalues to return    (6)
%     method   'auto', 'dense' or 'eigs'                   ('auto')
%     neg_tol  eigenvalues below -neg_tol count as negative (1e-6)
%     beta     preconditioner shift (eigs)                 (1)
%     verbose  print a summary                             (true)
%
% Outputs
%   lam    k smallest eigenvalues (ascending), translations/mean excluded
%   modes  matching modes, size [size(u), k], unit RMS, mean zero,
%          orthogonal to the translations; negative ones are energy-
%          decreasing directions at fixed mass
%   info   n_neg, n_near_zero (|lam| < neg_tol, marginal), stable, method,
%          rayleigh, n_translations, converged (false: eigs fell back to
%          Ritz values)

if nargin < 4, opts = struct(); end
defaults = struct('k', 6, 'method', 'auto', 'neg_tol', 1e-6, 'beta', 1, 'verbose', true);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}), opts.(names{i}) = defaults.(names{i}); end
end

sz = size(u);
n = numel(u);
method = opts.method;
if strcmp(method, 'auto')
    if n <= 5000, method = 'dense'; else, method = 'eigs'; end
end

% Jacobian and translation directions at u, from the Newton solver
[~, ninfo] = newton_polish_FCH(u, Grid, Para, struct('maxit', 0, 'verbose', false));
Jv = @(v) reshape(ninfo.Jv(u, reshape(v, sz)), n, 1);
T = reshape(ninfo.T, n, []);
Q = orth([ones(n, 1), T]);          % excluded directions: mean + translations
proj = @(v) v - Q*(Q'*v);
k = min(opts.k, n - size(Q, 2));

switch method
    case 'dense'
        J = zeros(n);
        e = zeros(n, 1);
        for j = 1:n
            e(j) = 1;
            J(:, j) = Jv(e);
            e(j) = 0;
        end
        Z = null(Q');                       % orthonormal basis of the allowed subspace
        Jr = Z' * J * Z;
        [W, L] = eig((Jr + Jr')/2);
        [lam, ord] = sort(diag(L));
        lam = lam(1:k);
        D = Z * W(:, ord(1:k));
        eigs_converged = true;

    case 'eigs'
        % P^{-1/2}: Fourier multiplier 1/(eps^2 k^2 + beta)
        Pmh_sym = 1 ./ (Para.epsilon^2*Grid.k2 + opts.beta);
        Pmh = @(v) reshape(real(ifftn(Pmh_sym .* fftn(reshape(v, sz)))), n, 1);
        % constraint Q'd = 0 on d = P^{-1/2} w  <=>  (P^{-1/2} Q)' w = 0
        Qt = orth(cell2mat(arrayfun(@(c) Pmh(Q(:, c)), 1:size(Q, 2), 'UniformOutput', false)));
        projt = @(w) w - Qt*(Qt'*w);
        % shift the excluded directions above the wanted (small) end of the
        % spectrum; keep it modest -- a large shift stretches the spectrum
        % and stalls Lanczos (1e3 failed to converge at 64^2)
        shift = 1;
        B = @(w) projt(Pmh(Jv(Pmh(projt(w))))) + shift*(Qt*(Qt'*w));
        % name-value pairs: a struct here is read as eigs' old-style
        % options (tol, issym, ...) and these settings would be ignored
        [W, L, flag] = eigs(B, n, k, 'smallestreal', 'IsFunctionSymmetric', true, ...
            'Tolerance', 1e-8, 'SubspaceDimension', max(40, 4*k), 'MaxIterations', 3000);
        lam = diag(L);
        eigs_converged = flag == 0;
        if ~eigs_converged
            % unconverged eigenvalues come back NaN; use Ritz values from the
            % returned vectors instead so the verdict never silently reads
            % NaN as "not negative"
            warning('stability_FCH:eigs', ['eigs did not converge all %d eigenvalues; ' ...
                'using Ritz values (verdict less certain)'], k);
            for c = 1:k
                lam(c) = (W(:, c)' * B(W(:, c))) / (W(:, c)' * W(:, c));
            end
        end
        [lam, ord] = sort(lam);
        D = zeros(n, k);
        for c = 1:k
            D(:, c) = Pmh(W(:, ord(c)));
        end

    otherwise
        error('stability_FCH: unknown method ''%s''', method);
end

% unit RMS, mean zero, orthogonal to translations
D = proj(D);
D = D ./ (vecnorm(D) / sqrt(n));
rayleigh = zeros(k, 1);
for c = 1:k
    rayleigh(c) = (D(:, c)' * Jv(D(:, c))) / (D(:, c)' * D(:, c));
end
modes = reshape(D, [sz, k]);

info.method = method;
info.n_translations = size(T, 2);
info.n_neg = sum(lam < -opts.neg_tol);
info.n_near_zero = sum(abs(lam) <= opts.neg_tol);
info.stable = info.n_neg == 0;
info.rayleigh = rayleigh;
info.converged = eigs_converged;

if opts.verbose
    if info.stable, verdict = 'STABLE'; else, verdict = 'UNSTABLE'; end
    fprintf('stability (%s, %d translation dir(s) excluded): %s, %d negative', ...
        method, info.n_translations, verdict, info.n_neg);
    if info.n_near_zero > 0
        fprintf(', %d near zero (marginal)', info.n_near_zero);
    end
    if ~info.converged
        fprintf(' [eigs NOT converged: Ritz values]');
    end
    fprintf('\n  smallest eigenvalues: %s\n', mat2str(lam', 3));
end
end
