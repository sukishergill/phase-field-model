function [u_vals, E_vals, m_vals, diagnostics] = arclength_FCH_func(u0, m0)
% ARCLENGTH_FCH_FUNC  Pseudo-arclength continuation of FCH steady states.
%
% [u_vals, E_vals, m_vals, diagnostics] = arclength_FCH_func(u0, m0)
%
% Uses a cosine-only (even-symmetric) representation: u is assumed even
% about the domain center, so its Fourier coefficients ck are real and
% palindromic (ck(N+2-j) = ck(j)). The state is stored as the reduced,
% real vector c = ck(1:N/2+1) of independent cosine amplitudes rather than
% the full complex ck. This eliminates the translation mode by
% construction -- du/dx is odd, so translation simply isn't a direction
% available within the even/cosine subspace -- instead of relying on
% bordering to paper over it after the fact.
%
% Inputs:
%   u0 - initial profile (physical space, column vector of length N),
%        typically the SAV time-stepper's output (e.g. Results.uu{end}).
%        u0 does not need to already be a converged steady state: it is
%        pre-relaxed here with a semi-implicit gradient descent before
%        continuation starts (see gradient_descent_cos below), since ANY
%        spatially constant profile is an *exact* root of compute_mu, and
%        starting the corrector too far from the true non-trivial branch
%        can converge it onto that trivial constant branch instead.
%   m0 - initial mass (target mean value). The branch is continued in
%        both increasing and decreasing m from this starting point.
%
% Outputs:
%   u_vals - 1xK cell array of profiles along the combined branch, both
%            directions spliced together into one continuous sequence
%            (ordered from the decreasing-m end through to the
%            increasing-m end).
%   E_vals - 1xK vector of energies, aligned with u_vals.
%   m_vals - 1xK vector of masses, aligned with u_vals.
%   diagnostics - struct with fields:
%       gd_resid        residual history from the initial gradient descent
%       polish_resid    residual history from the fixed-mass Newton polish
%                       that follows gradient descent (see
%                       newton_polish_cos below)
%       t0              initial tangent vector used (increasing-m sign)
%       stability       1xK vector aligned with m_vals/E_vals: +1 where
%                       linearly stable under u_t = mu (every eigenvalue
%                       of the Jacobian, restricted to the
%                       mass-conserving subspace, has negative real
%                       part), -1 otherwise. Only the sign is kept, not
%                       the full spectrum (see extend_branch below).
%       inc, dec        per-direction diagnostic structs from
%                       extend_branch (stop_reason, steps, final_ds)

eps = 0.1;      eta = eps^2;

% generate grid
dim = 1;
L = 2*pi;
N = numel(u0);      % must match whatever grid size u0 was generated with
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);
Nc = N/2 + 1;       % number of independent cosine (even) Fourier coefficients

ck = fft(u0) / N;
c = reduce_cos(ck, N);      % restrict to the cosine (even) subspace

% The SAV time-stepper's output may not be relaxed all the way to a
% genuine non-trivial steady state. That matters here more than it might
% seem: ANY spatially constant profile is an *exact* root of compute_mu,
% for every value of the mean -- so if the corrector starts too far from
% the true non-trivial branch, Newton's iteration can converge onto that
% trivial constant branch instead (this is the "energy locked at -eta/4"
% failure mode). Pre-relax with a robust (if slower) gradient descent
% before handing the profile to Newton.
gd_maxit = 4000;
gd_tol = 1e-3;
[c, gd_resid] = gradient_descent_cos(c, eps, eta, Grid, gd_maxit, gd_tol);

% Gradient descent is robust far from a solution but can plateau well
% short of full convergence -- a genuine limitation of this simple
% first-order scheme, not a step-count or step-size-floor issue (more
% iterations or a much smaller minimum step size don't move the plateau
% at all). Follow up with a plain, fixed-mass Newton polish: it converges
% much faster once reasonably close, which is exactly the regime
% gradient descent's plateau should leave us in.
polish_maxit = 50;
polish_tol = 1e-8;
[c, polish_resid] = newton_polish_cos(c, m0, eps, eta, Grid, polish_maxit, polish_tol);

ck = expand_cos(c, N);
u = real(ifft(ck * N));
m = m0;

E = SSAV_FCH_helpers.compute_E_direct(u, ck*N, eps, eta, Grid);       %#ok<NASGU> initial E(u), unused after this (extend_branch recomputes it)

% ------------------------------------------------------------------------
% Initial tangent. A "correct" null-vector from the SVD of the bordered
% Jacobian is only well-defined when there's a single, clearly-isolated
% smallest singular value; near clustered near-degenerate points (several
% modes hovering close to marginal at once) which direction gets picked is
% numerically fragile, and can carry a built-in bias toward jumping onto
% the wrong branch on the very first predictor step. A random tangent
% avoids that particular bias, but being random it can occasionally point
% toward the wrong nearby branch anyway -- initialize_t instead computes a
% deterministic tangent from the actual Jacobian at the actual current
% profile (see its comments), which has neither problem.
t0 = initialize_t(u, c, eps, eta, Grid);

% ------------------------------------------------------------------------
% Extend the branch in both directions from the same relaxed starting
% point, then splice the two histories into one continuous array.

ds_start = 1E-6;
ds_max = 5E-4;
max_steps = 10000;
loop_tol = 20 * ds_max;     % closed-loop detection threshold, in (c,m)-space

% Two-sided bound: folds can send either direction's branch past either
% boundary (an "increasing m" branch can fold over and start decreasing,
% and vice versa), so both directions use the same absolute [0, 1] bound
% rather than a single direction-specific target.
m_bounds = [0, 1];

fprintf('=== extending in increasing m ===\n');
[data_inc, uvals_inc, stab_inc, diag_inc] = extend_branch(c, m, u, t0, eps, eta, Grid, N, Nc, ...
    ds_start, ds_max, max_steps, m_bounds, loop_tol);

fprintf('=== extending in decreasing m ===\n');
[data_dec, uvals_dec, stab_dec, diag_dec] = extend_branch(c, m, u, -t0, eps, eta, Grid, N, Nc, ...
    ds_start, ds_max, max_steps, m_bounds, loop_tol);

% Splice: reversed decreasing branch (from its far end back to the shared
% starting point) followed by the increasing branch forward from the
% start (the starting point itself is dropped from the second half since
% it's already the last column of the first half).
data = [fliplr(data_dec), data_inc(:, 2:end)];
u_vals = [fliplr(uvals_dec), uvals_inc(2:end)];
stability = [fliplr(stab_dec), stab_inc(2:end)];

m_vals = data(1, :);
E_vals = data(2, :);

diagnostics = struct( ...
    'gd_resid', gd_resid, ...
    'polish_resid', polish_resid, ...
    't0', t0, ...
    'stability', stability, ...
    'inc', diag_inc, ...
    'dec', diag_dec);

end

function u = compute_u(Grid, m, eps)

u = tanh((Grid.x + pi*(m + 1)/2) / eps) - ...
    tanh((Grid.x- pi*(m + 1)/2) / eps) - 1;

% u = m*ones(size(Grid.x));         % test case

end


function Hv = fch_jac(u, v, eps, eta, Grid)
% Action of Jacobian on one vector

[~, dF, d2F] = SSAV_FCH_helpers.compute_F(u, 0);
d3F = 6*u;

what = -eps^2*Grid.k2.*fft(u)/Grid.N - fft(dF)/Grid.N;
w = real(ifft(what*Grid.N));

dv = real(ifft(Grid.N * v));

dwhat = -eps^2*Grid.k2.*v - fft(d2F.*dv)/Grid.N;
dw = real(ifft(Grid.N * dwhat));


% Hv = -eps^2*Grid.k2.*w - w.*fft(d2F.*dv)/Grid.N + eta*w;
Hv = -eps^2*Grid.k2.*dwhat ...
     - fft(d2F.*dw)/Grid.N ...
     - fft(d3F.*dv.*w)/Grid.N ...
     + eta*dwhat;
Hv = -Grid.k2.*Hv;

end

% Build the full Jacobian
function [J, J_FD] = fch_jac_matrix(u, ck, eps, eta, Grid)
%FCH_JAC_MATRIX  Dense FCH Hessian/Jacobian matrix assembled by matvecs.
%   Intended for small 1D or small 2D problems.

J = zeros(Grid.N);
J_FD = J;

e = zeros(size(u));
h = 1e-6;

for i = 1:Grid.N
    e(i) = 1;
    col = fch_jac(u, e, eps, eta, Grid);
    J_FD(:,i) = (compute_mu(ck+e*h, eps, eta, Grid) - ...
        compute_mu(ck-e*h, eps, eta, Grid))/(2*h);
    J(:,i) = col(:);
    e(i) = 0;

end
end
%------------------

function mu = compute_mu(ck, eps, eta, Grid)

u = real(ifft(Grid.N*ck));
[~, dF, d2F] = SSAV_FCH_helpers.compute_F(u, 0);

what = -eps^2.*Grid.k2.*ck - fft(dF)/Grid.N;
w = real(ifft(what*Grid.N));
mu = -eps^2*Grid.k2.*what - fft(d2F.*w)/Grid.N + eta*what;
mu = -Grid.k2.*mu;
end

function t = initialize_t(u, c, eps, eta, Grid)
% Deterministic initial tangent via implicit differentiation of the
% constrained system at the CURRENT profile (u, c) -- not a random
% direction, and not a finite difference of the idealized analytic
% tanh-front formula (which only reflects the simplest possible branch,
% not whichever branch (u, c) actually sits on).
%
% The mass parameter m enters only through the constraint c(1) = m, never
% through compute_mu directly, so differentiating the whole constrained
% system -- c(1) = m, compute_mu_cos(c) = 0 on rows 2:end -- with respect
% to m gives:
%   J(2:end,:) * dc/dm = 0      (interior/physics rows)
%   dc(1)/dm = 1                (from differentiating c(1) = m)
% This is a well-posed (Nc x Nc) linear system for dc/dm at a regular
% (non-fold) point; solve it directly, then set dm/ds = 1 as the
% reference direction before normalizing. Being deterministic, this has
% no risk of an unlucky draw pointing toward the wrong nearby branch, the
% way a random initial tangent occasionally can.

Nc = Grid.N/2 + 1;
J = fch_jac_matrix_cos(u, c, eps, eta, Grid);

e1_row = zeros(1, Nc);   e1_row(1) = 1;
A0 = [J(2:end,:); e1_row];
rhs = [zeros(Nc-1, 1); 1];

dc_dm = A0 \ rhs;

t = [dc_dm; 1];
t = t / norm(t, 2);

end

%------------------------------------------------------------------------
% cosine (even-symmetric) reduction: since the AC/FCH equations have no
% explicit x-dependence and only involve even-order derivatives and even
% nonlinearities, they are equivariant under x -> -x, i.e. they map the
% subspace of even u (real, palindromic ck) to itself. That lets us work
% entirely with the reduced, real vector of independent cosine amplitudes
% c = ck(1:N/2+1), instead of the full redundant complex ck.

function ck_full = expand_cos(c, N)
% Build the full length-N ck from the reduced cosine coefficients c
% (length N/2+1). Grid.x is centered at x=0 but is a 1-indexed grid
% x_i = i*dx - L/2, so the true reflection point sits half a period off
% from the standard (0-indexed) FFT phase origin. That means an even u
% does NOT simply give a real, palindromic ck here -- ck picks up a known
% linear phase exp(-2*pi*i*k*(N/2-1)/N). We build the real/palindromic
% ("de-rotated") array first, then undo that phase to get the actual ck.
ck_derot = zeros(N, 1);
ck_derot(1:N/2+1) = c;
ck_derot(N/2+2:N) = c(N/2:-1:2);
k = (0:N-1)';
phase = exp(2i*pi*k*(N/2-1)/N);
ck_full = ck_derot .* conj(phase);
end

function c = reduce_cos(ck_full, N)
% Restrict a full length-N ck (assumed to come from an even u) down to the
% reduced cosine coefficients c (length N/2+1). First undoes the grid's
% phase offset (see expand_cos) so the result is real & palindromic, then
% averages the two mirrored entries for robustness (e.g. against roundoff).
k = (0:N-1)';
phase = exp(2i*pi*k*(N/2-1)/N);
ck_derot = ck_full .* phase;
c = zeros(N/2+1, 1);
c(1) = real(ck_derot(1));
c(N/2+1) = real(ck_derot(N/2+1));
c(2:N/2) = real(0.5*(ck_derot(2:N/2) + ck_derot(N:-1:N/2+2)));
end

function mu_c = compute_mu_cos(c, eps, eta, Grid)
ck_full = expand_cos(c, Grid.N);
mu_full = compute_mu(ck_full, eps, eta, Grid);
mu_c = reduce_cos(mu_full, Grid.N);
end

function Hv = fch_jac_cos(u, vc, eps, eta, Grid)
v_full = expand_cos(vc, Grid.N);
Hv_full = fch_jac(u, v_full, eps, eta, Grid);
Hv = reduce_cos(Hv_full, Grid.N);
end

function [J, J_FD] = fch_jac_matrix_cos(u, c, eps, eta, Grid)
% Dense (N/2+1)x(N/2+1) real Jacobian in the reduced cosine basis.

Nc = Grid.N/2 + 1;
J = zeros(Nc);
J_FD = J;

e = zeros(Nc, 1);
h = 1e-6;

for i = 1:Nc
    e(i) = 1;
    col = fch_jac_cos(u, e, eps, eta, Grid);
    J_FD(:,i) = (compute_mu_cos(c+e*h, eps, eta, Grid) - ...
        compute_mu_cos(c-e*h, eps, eta, Grid))/(2*h);
    J(:,i) = col(:);
    e(i) = 0;
end
end

function [c, resid_hist] = gradient_descent_cos(c, eps, eta, Grid, maxit, tol)
% Semi-implicit ("stabilized") gradient descent toward a genuine FCH
% steady state at fixed mass, starting from c.
%
% compute_mu_cos(c) is (proportional to) the mass-conserving L2 gradient
% of the energy. Its linear part is 6th-order in derivatives and very
% stiff (eps^4*k^6), so a naive explicit step would need an impractically
% tiny dtau for stability. Instead the linear part is treated implicitly
% -- cheap here since it's diagonal in this cosine/Fourier basis -- while
% the nonlinear remainder is explicit, with adaptive step-size control:
% grow dtau on an accepted (energy-decreasing) step, halve and retry on
% any step that would increase the energy. The mean/mass mode needs no
% special handling: both the linear operator and compute_mu are
% structurally zero there (Grid.k2(1)=0), so c(1) never moves.

Nc = Grid.N/2 + 1;
k2r = Grid.k2(1:Nc);
Lop = -eps^4*k2r.^3 + eta*eps^2*k2r.^2;

ck = expand_cos(c, Grid.N);
u = real(ifft(ck * Grid.N));
E = SSAV_FCH_helpers.compute_E_direct(u, ck*Grid.N, eps, eta, Grid);

dtau = 1e-3;
dtau_min = 1e-8;
resid_hist = zeros(maxit, 1);

for it = 1:maxit
    mu = compute_mu_cos(c, eps, eta, Grid);
    resid_hist(it) = norm(mu);

    if mod(it, 100) == 0
        fprintf("  gradient descent it %d: norm(mu)=%e  dtau=%e\n", it, resid_hist(it), dtau);
    end

    if resid_hist(it) < tol
        resid_hist = resid_hist(1:it);
        return
    end

    % Near a critical point the true energy decrease per step shrinks with
    % the (already small) residual, and can end up smaller than the
    % floating-point noise in the energy evaluation itself -- the
    % step-halving would otherwise shrink dtau toward underflow forever
    % without making progress. Once dtau bottoms out, stop early and let
    % Newton (which uses proper second-order information) finish the job
    % from here rather than wasting iterations.
    if dtau < dtau_min
        fprintf("  gradient descent stalled at it %d (dtau underflowed): norm(mu)=%e -- stopping early, handing off to Newton\n", ...
            it, resid_hist(it));
        resid_hist = resid_hist(1:it);
        return
    end

    c_trial = (c + dtau*(mu - Lop.*c)) ./ (1 - dtau*Lop);
    ck_trial = expand_cos(c_trial, Grid.N);
    u_trial = real(ifft(ck_trial * Grid.N));
    E_trial = SSAV_FCH_helpers.compute_E_direct(u_trial, ck_trial*Grid.N, eps, eta, Grid);

    if E_trial > E
        dtau = dtau * 0.5;
        continue
    end

    c = c_trial;
    E = E_trial;
    dtau = dtau * 1.1;
end

end

function [c, resid_hist] = newton_polish_cos(c, m, eps, eta, Grid, maxit, tol)
% Plain, fixed-mass Newton polish: solve compute_mu_cos(c) = 0 on rows
% 2:end (row 1 is structurally always zero -- the mean-mode redundancy
% from Grid.k2(1)=0 -- so there's nothing to solve for there), holding
% c(1) = m fixed throughout. There's no arclength/tangent row here,
% unlike extend_branch's corrector: this just refines a single point
% before continuation starts, it doesn't step along the branch.
%
% This exists because gradient_descent_cos, despite being globally
% robust, can plateau well short of full convergence -- a genuine
% limitation of that simple first-order scheme (more iterations or a much
% smaller step-size floor don't move the plateau at all, since it isn't a
% budget or step-size problem). Newton converges much faster once
% reasonably close to a root, which is exactly the regime gradient
% descent's plateau should leave us in.

Nc = Grid.N/2 + 1;
c(1) = m;
resid_hist = zeros(maxit, 1);

for it = 1:maxit
    mu = compute_mu_cos(c, eps, eta, Grid);
    resid_hist(it) = norm(mu(2:end));

    if resid_hist(it) < tol
        resid_hist = resid_hist(1:it);
        return
    end

    ck = expand_cos(c, Grid.N);
    u = real(ifft(ck * Grid.N));
    J = fch_jac_matrix_cos(u, c, eps, eta, Grid);

    dc = zeros(Nc, 1);
    dc(2:end) = J(2:end, 2:end) \ (-mu(2:end));
    c = c + dc;     % c(1) unchanged: dc(1) = 0, mass stays exactly m

end

end

function [data_out, u_vals_out, stability_out, diag] = extend_branch(c, m, u, t, eps, eta, Grid, ...
    N, Nc, ds_start, ds_max, max_steps, m_bounds, loop_tol)
% Extends the FCH continuation branch from (c, m) along tangent t.
%
% Step size is adaptive, keyed off how many Newton iterations the
% corrector needed:
%   - if the corrector fails to converge within its iteration cap, the
%     step is REJECTED: ds is quartered and the same predictor/corrector
%     is retried from the last accepted point (never advances on a
%     failed step);
%   - on an accepted step, ds grows after a fast (<=3 iteration)
%     convergence, shrinks after a slow (>=10 iteration) one, and is left
%     alone in between.
%
% A converged Newton step is also rejected (same quarter-ds retry as an
% outright failure) if the actual displacement is wildly disproportionate
% to what was requested/recently seen -- severe local ill-conditioning can
% let Newton converge to a technically-valid but physically wrong root
% (i.e. the wrong branch), and this can happen anywhere along the branch,
% not just where the history so far has been badly behaved.
%
% Stops when:
%   - m leaves [m_bounds(1), m_bounds(2)] -- folds can send either
%     direction's branch past either boundary (an "increasing m" branch
%     can fold over and start decreasing, and vice versa), so this is
%     checked as an absolute two-sided bound, not relative to which
%     direction we started extending in, or
%   - the branch closes into a loop: the current (c, m) state comes back
%     within loop_tol of the starting point, after having moved at least
%     a little way from it first (so this can't fire immediately), or
%   - max_steps is reached, or
%   - ds underflows while retrying a rejected step (a robustness stop).
%
% stability_out is a 1xK vector aligned with data_out's columns: +1 where
% the point is linearly stable under the mass-conserving flow u_t = mu
% (every eigenvalue of the Jacobian, restricted to the mass-conserving
% subspace J(2:end,2:end), has negative real part), -1 otherwise. Just the
% sign is kept, not the full spectrum, to keep this cheap to store.
%
% diag is a struct with fields stop_reason (string), steps (number of
% accepted steps taken), and final_ds (ds at the point extension stopped).

ck = expand_cos(c, N);
E = SSAV_FCH_helpers.compute_E_direct(u, ck*N, eps, eta, Grid);
data_out = [m; E];
u_vals_out = cell(1, 1);
u_vals_out{1} = u;

J0 = fch_jac_matrix_cos(u, c, eps, eta, Grid);
stability_out = 2*all(real(eig(J0(2:end,2:end))) < 0) - 1;

c0 = c;
m0 = m;
ds = ds_start;
ds_min = ds_start * 1e-2;
s_accum = 0;
min_loop_s = 30 * ds_max;
prev_disp = NaN;    % no reference displacement yet for the very first step

diag.stop_reason = 'max_steps_reached';    % default; overridden below if
diag.steps = max_steps;                    % something else stops us first
diag.final_ds = ds;

for step = 1:max_steps

    accepted = false;
    while ~accepted

        % predictor
        c_curr = c + ds * t(1:end-1);
        m_curr = m + ds * t(end);
        ck_curr = expand_cos(c_curr, N);
        u_curr = real(ifft(ck_curr * N));

        % corrector (Newton)
        its = 0;
        ndx = 1;
        converged_tight = false;
        J_final = [];   % Jacobian at the converged point, once we have one --
                        % reused below for the tangent update and stability
                        % check instead of rebuilding it from scratch.

        while its < 100 && ndx > 1e-6

            F_curr = compute_mu_cos(c_curr, eps, eta, Grid);

            [J, ~] = fch_jac_matrix_cos(u_curr, c_curr, eps, eta, Grid);

            Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
            Jm = zeros(Nc-1, 1);

            G = [J(2:end,:), Jm; Fm; [(1/Nc)*t(1:end-1); t(end)]'];

            b = [F_curr(2:end); c_curr(1) - m_curr;...
                (1/Nc)*(c_curr - c)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

            dx = G \ (-b);

            c_new = c_curr + dx(1:end-1);
            ck_new = expand_cos(c_new, N);
            u_new = real(ifft(ck_new * N));

            m_new = m_curr + dx(end);
            m_curr = m_new;
            c_curr = c_new;
            u_curr = u_new;

            if norm(dx) < 1e-8 && norm(b) < 1e-8
                converged_tight = true;
                % One extra Jacobian build here, at the now-converged point
                % -- this replaces the separate rebuild that used to happen
                % in the "update tangent" step below and the stability
                % check, rather than adding a third build on top of them.
                J_final = fch_jac_matrix_cos(u_curr, c_curr, eps, eta, Grid);
                break
            end
            ndx = norm(dx);
            its = its + 1;

        end

        success = converged_tight || (ndx <= 1e-6);

        if success && isempty(J_final)
            % Exited via the loose ndx<=1e-6 threshold rather than the
            % tight break above -- still need a Jacobian at the converged
            % point for reuse below.
            J_final = fch_jac_matrix_cos(u_curr, c_curr, eps, eta, Grid);
        end

        % Even when Newton nominally converges, severe local
        % ill-conditioning can let it converge onto a completely
        % different (wrong) branch instead of the intended nearby point.
        % Warning sign: the actual displacement is wildly disproportionate
        % to the requested ds / the last accepted step's displacement
        % (whichever is larger, so this doesn't false-positive right after
        % ds has legitimately grown). Iteration count alone isn't used
        % here -- near a fold, genuinely correct convergence can just take
        % many iterations, so that on its own isn't a reliable sign of
        % having landed on the wrong branch. A large jump is treated as a
        % rejected step, exactly like an outright convergence failure.
        this_disp = sqrt((1/Nc)*sum((c_curr - c).^2) + (m_curr - m)^2);
        jump_ref = max(ds, prev_disp);
        suspicious = success && ~isnan(prev_disp) && this_disp > 10*jump_ref;

        if success && ~suspicious
            accepted = true;
            prev_disp = this_disp;
        else
            if suspicious
                fprintf(['  extend_branch: rejecting suspicious step at step %d ' ...
                    '(its=%d, disp=%e vs ref=%e) -- likely branch switch\n'], ...
                    step, its, this_disp, jump_ref);
            end
            % Newton failed to converge within the iteration cap, or the
            % step looked suspicious above -- reject the step (don't
            % advance), shrink ds, and retry from the same last-accepted
            % point.
            ds = ds / 4;
            if ds < ds_min
                fprintf('  extend_branch: ds underflowed at step %d -- stopping\n', step);
                diag.stop_reason = 'ds_underflow';
                diag.steps = step - 1;
                diag.final_ds = ds;
                return
            end
        end
    end

    % commit the accepted step
    c = c_curr;
    ck = expand_cos(c, N);
    u = u_curr;
    m = m_curr;
    s_accum = s_accum + ds;

    E = SSAV_FCH_helpers.compute_E_direct(u, ck * N, eps, eta, Grid);
    data_out = [data_out, [m; E]];
    u_vals_out{end+1} = u; %#ok<AGROW>

    % Stability at this point: reuse J_final (already built at this exact
    % (u, c) by the corrector above) instead of rebuilding it.
    is_stable = all(real(eig(J_final(2:end,2:end))) < 0);
    stability_out(end+1) = 2*is_stable - 1; %#ok<AGROW>

    % Adaptive step-size control based on how many Newton iterations it
    % took, graduated rather than a hard binary grow/shrink -- a step that
    % converges in 9 iterations is "successful but slower than ideal" and
    % should ease off gently, not get treated the same as one that took 4,
    % nor get the same sharp cut as one that took 15. The harshest cuts
    % are reserved for genuine non-convergence, handled separately above
    % (reject-and-quarter-ds retry).
    if its <= 3
        ds = min(ds * 1.5, ds_max);     % fast convergence -> grow
    elseif its <= 6
        % comfortably within the normal range -- leave ds alone
    elseif its <= 9
        ds = max(ds * 0.75, ds_min);    % converging, but slower than ideal
    else
        ds = max(ds * 0.5, ds_min);     % slow convergence -> shrink more
    end

    if mod(step, 100) == 0
        fprintf('  step %5d: m=%.6f  E=%.6f  its=%3d  ds=%e\n', step, m, E, its, ds);
    end

    % stop if m has left [m_bounds(1), m_bounds(2)] -- checked as an
    % absolute two-sided bound (not relative to m0/which direction we
    % started extending in), since a fold can send either direction's
    % branch past either boundary.
    if m <= m_bounds(1) || m >= m_bounds(2)
        fprintf('  extend_branch: m left [%.4f, %.4f] (m=%.6f) at step %d\n', ...
            m_bounds(1), m_bounds(2), m, step);
        diag.stop_reason = 'm_bound_reached';
        diag.steps = step;
        diag.final_ds = ds;
        break
    end

    % stop if the branch has closed into a loop back to the start
    if s_accum > min_loop_s
        dist_to_start = sqrt(sum((c - c0).^2) + (m - m0)^2);
        if dist_to_start < loop_tol
            fprintf('  extend_branch: closed loop detected at step %d (dist to start=%e, s=%.4f)\n', ...
                step, dist_to_start, s_accum);
            diag.stop_reason = 'closed_loop';
            diag.steps = step;
            diag.final_ds = ds;
            break
        end
    end

    % update tangent -- reuse J_final (built at this exact (u, c) by the
    % corrector above) instead of rebuilding the Jacobian a third time.
    Jm = zeros(Nc-1, 1);
    Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
    A = [J_final(2:end,:), Jm; Fm; [(1/Nc)*t(1:end-1); t(end)]'];
    t = A \ [zeros(Nc, 1); 1];
    t = t / sqrt((1/Nc)*sum(t(1:end-1).^2) + t(end)^2);

end

if strcmp(diag.stop_reason, 'max_steps_reached')
    diag.final_ds = ds;
end

end
