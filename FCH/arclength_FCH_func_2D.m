function [u_vals, E_vals, m_vals, diag1nostics] = arclength_FCH_func_2D(u0, m0, eps, eta, Grid)
% ARCLENGTH_FCH_FUNC  Pseudo-arclength continuation of FCH steady states.
%
% [u_vals, E_vals, m_vals, diag1nostics] = arclength_FCH_func(u0, m0)
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
%   diag1nostics - struct with fields:
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
%       inc, dec        per-direction diag1nostic structs from
%                       extend_branch (stop_reason, steps, final_ds)

% eps = 0.1;      eta = eps^2;

tic;

% generate grid
% dim = 2;
% Lx = 2*pi;      Ly = 2*pi;
% L = [Lx, Ly];
% 
% Nx = size(u0, 2);       Ny = size(u0, 1);
% % N = [Nx, Ny];
% 
% Grid = SSAV_FCH_helpers.generate_Grid(L, [Nx, Ny], dim);
N = prod(Grid.N);

Ncx = Grid.N(1)/2 + 1;
Ncy = Grid.N(2)/2 + 1;

Grid.Nc = [Ncx, Ncy];       % number of independent cosine (even) Fourier coefficients

% Weights of the arclength inner product on the reduced coefficients.
% Each reduced coefficient c(i,j) stands for 1, 2 or 4 entries of the full
% ck (interior indices are mirrored in each direction), so by Parseval
%   sum(Grid.wc .* c.^2) = mean((u - mean(u)).^2)
% with the mean mode c(1) given weight 0: it equals m, which the arclength
% norm counts separately. One unit of arclength is then an RMS change of 1
% in the profile (or a change of 1 in m), independent of the grid size.
% (The previous (1/Nc)*c'*c weighting made the profile part ~Nc times too
% light, so the arclength constraint was in effect "dm = ds" and did not
% limit how far the profile could move in one step.)
w1 = [1; 2*ones(Ncx - 2, 1); 1];
w2 = [1; 2*ones(Ncy - 2, 1); 1];
Grid.wc = reshape(w1 * w2', [], 1);
Grid.wc(1) = 0;

ck = fftn(u0) / prod(N);
c = reduce_cos(ck, Grid.N, Grid.Nc);      % restrict to the cosine (even) subspace

% The SAV time-stepper's output may not be relaxed all the way to a
% genuine non-trivial steady state. That matters here more than it might
% seem: ANY spatially constant profile is an *exact* root of compute_mu,
% for every value of the mean -- so if the corrector starts too far from
% the true non-trivial branch, Newton's iteration can converge onto that
% trivial constant branch instead (this is the "energy locked at -eta/4"
% failure mode). Pre-relax with a robust (if slower) gradient descent
% before handing the profile to Newton.
gd_maxit = 5000;
gd_tol = 1e-3;
[c, gd_resid] = gradient_descent_cos(c, eps, eta, Grid, gd_maxit, gd_tol);

% Gradient descent is robust far from a solution but can plateau well
% short of full convergence -- a genuine limitation of this simple
% first-order scheme, not a step-count or step-size-floor issue (more
% iterations or a much smaller minimum step size don't move the plateau
% at all). Follow up with a plain, fixed-mass Newton polish: it converges
% much faster once reasonably close, which is exactly the regime
% gradient descent's plateau should leave us in.
polish_maxit = 100;
polish_tol = 1e-8;
[c, polish_resid] = newton_polish_cos(c, m0, eps, eta, Grid, polish_maxit, polish_tol);

ck = expand_cos(c, Grid.N, Grid.Nc);
u = real(ifftn(ck * N));
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
ds_max = 0.05;              % in the weighted norm above: RMS profile change per step
max_steps = 5000;
loop_tol = ds_max;          % closed-loop detection threshold (weighted norm): within one step of the start

% Two-sided bound: folds can send either direction's branch past either
% boundary (an "increasing m" branch can fold over and start decreasing,
% and vice versa), so both directions use the same absolute [0, 1] bound
% rather than a single direction-specific target.
m_bounds = [0, 1];

fprintf('=== extending in increasing m ===\n');
[data_inc, uvals_inc, stab_inc, diag1_inc] = extend_branch(c, m, u, t0, eps, eta, Grid, ...
    ds_start, ds_max, max_steps, m_bounds, loop_tol);

fprintf('=== extending in decreasing m ===\n');
[data_dec, uvals_dec, stab_dec, diag1_dec] = extend_branch(c, m, u, -t0, eps, eta, Grid, ...
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

diag1nostics = struct( ...
    'gd_resid', gd_resid, ...
    'polish_resid', polish_resid, ...
    't0', t0, ...
    'stability', stability, ...
    'inc', diag1_inc, ...
    'dec', diag1_dec);

elapsedTime = toc;
disp(['Continuation run time:', num2str(elapsedTime), ' seconds.']);

end



function Hv = fch_jac(u, v, eps, eta, Grid)
% Action of Jacobian on one vector

N = prod(Grid.N);

[~, dF, d2F] = SSAV_FCH_helpers.compute_F(u, 0);
d3F = 6*u;

what = -eps^2*Grid.k2.*fftn(u)/N - fftn(dF)/N;
w = real(ifftn(what*N));

dv = real(ifftn(N * v));

dwhat = -eps^2*Grid.k2.*v - fftn(d2F.*dv)/N;
dw = real(ifftn(N * dwhat));

Hv = -eps^2*Grid.k2.*dwhat ...
     - fftn(d2F.*dw)/N ...
     - fftn(d3F.*dv.*w)/N ...
     + eta*dwhat;
Hv = -Grid.k2.*Hv;

end

% Build the full Jacobian
function J = fch_jac_matrix(u, ck, eps, eta, Grid)
%FCH_JAC_MATRIX  Dense FCH Hessian/Jacobian matrix assembled by matvecs.
%   Intended for small 1D or small 2D problems.

J = zeros(Grid.N);

e = zeros(size(u));

for i = 1:Grid.N
    e(i) = 1;
    col = fch_jac(u, e, eps, eta, Grid);
    J(:,i) = col(:);
    e(i) = 0;

end
end
%------------------

function mu = compute_mu(ck, eps, eta, Grid)
N = prod(Grid.N);
u = real(ifftn(N*ck));
[~, dF, d2F] = SSAV_FCH_helpers.compute_F(u, 0);

what = -eps^2.*Grid.k2.*ck - fftn(dF)/N;
w = real(ifftn(what*N));
mu = -eps^2*Grid.k2.*what - fftn(d2F.*w)/N + eta*what;
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

Nc = prod(Grid.Nc);
J = fch_jac_matrix_cos(u, c, eps, eta, Grid);

e1_row = zeros(1, Nc);   e1_row(1) = 1;
A0 = [J(2:end,:); e1_row];
rhs = [zeros(Nc-1, 1); 1];

dc_dm = A0 \ rhs;

t = [dc_dm; 1];
t = t / sqrt(sum(Grid.wc .* t(1:end-1).^2) + t(end)^2);     % same norm as extend_branch

end

%------------------------------------------------------------------------
% cosine (even-symmetric) reduction: since the AC/FCH equations have no
% explicit x-dependence and only involve even-order derivatives and even
% nonlinearities, they are equivariant under x -> -x, i.e. they map the
% subspace of even u (real, palindromic ck) to itself. That lets us work
% entirely with the reduced, real vector of independent cosine amplitudes
% c = ck(1:N/2+1), instead of the full redundant complex ck.

function ck_full = expand_cos(c, N, Nc)
% Build the full length-N ck from the reduced cosine coefficients c
% (length N/2+1). Grid.x is centered at x=0 but is a 1-indexed grid
% x_i = i*dx - L/2, so the true reflection point sits half a period off
% from the standard (0-indexed) FFT phase origin. That means an even u
% does NOT simply give a real, palindromic ck here -- ck picks up a known
% linear phase exp(-2*pi*i*k*(N/2-1)/N). We build the real/palindromic
% ("de-rotated") array first, then undo that phase to get the actual ck.

c = reshape(c, Nc(1), Nc(2));

ck_derot = zeros(N(1), N(2));
ck_derot(1:Nc(1), 1:Nc(2)) = c;

ck_derot(Nc(1)+1:N(1), 1:Nc(2)) = c(Nc(1)-1:-1:2, :);
ck_derot(:, Nc(2)+1:N(2)) = ck_derot(:, Nc(2)-1:-1:2);

kx = (0:N(1)-1);
ky = (0:N(2)-1)';

phase_x = exp(2i*pi*kx*(N(1)/2-1)/N(1));
phase_y = exp(2i*pi*ky*(N(2)/2-1)/N(2));

phase = phase_y * phase_x;

ck_full = ck_derot .* conj(phase);
end

function c = reduce_cos(ck_full, N, Nc)
% Restrict a full length-N ck (assumed to come from an even u) down to the
% reduced cosine coefficients c (length N/2+1). First undoes the grid's
% phase offset (see expand_cos) so the result is real & palindromic, then
% averages the two mirrorexd entries for robustness (e.g. against roundoff).
kx = (0:N(1)-1);
ky = (0:N(2)-1)';

phase_x = exp(2i*pi*kx*(N(1)/2-1)/N(1));
phase_y = exp(2i*pi*ky*(N(2)/2-1)/N(2));

phase = phase_y * phase_x;

D = ck_full .* phase;

A = D(1:Nc(1), 1:Nc(2));

% Reflection indices.
ix = 1:Nc(1);
iy = 1:Nc(2);

ix(2:Nc(1)-1) = N(1) - (2:Nc(1)-1) + 2;
iy(2:Nc(2)-1) = N(2) - (2:Nc(2)-1) + 2;

% Average reflected x coefficients.
if Nc(1) > 2
    idx = 2:Nc(1)-1;
    A(idx,:) = 0.5*(A(idx,:) + D(ix(idx),1:Nc(2)));
end

% Average reflected y coefficients.
if Nc(2) > 2
    idx = 2:Nc(2)-1;
    A(:,idx) = 0.5*(A(:,idx) + D(1:Nc(1),iy(idx)));
end

c = real(A(:));

end

function C = reduce_cos_batch(CK, N, Nc)
% Vectorized reduction for a batch of Fourier coefficient arrays.

% B = size(CK,3);

kx = 0:N(1)-1;
ky = (0:N(2)-1)';

phase_x = exp(2i*pi*kx*(N(1)/2-1)/N(1));
phase_y = exp(2i*pi*ky*(N(2)/2-1)/N(2));

phase = reshape(phase_y*phase_x,N(1),N(2),1);

D = CK .* phase;

A = D(1:Nc(1),1:Nc(2),:);

ix = 1:Nc(1);
iy = 1:Nc(2);

if Nc(1)>2
    ix(2:Nc(1)-1) = N(1)-(2:Nc(1)-1)+2;
    idx = 2:Nc(1)-1;

    A(idx,:,:) = 0.5*(A(idx,:,:) + D(ix(idx),1:Nc(2),:));
end

if Nc(2)>2
    iy(2:Nc(2)-1) = N(2)-(2:Nc(2)-1)+2;
    idx = 2:Nc(2)-1;

    A(:,idx,:) = 0.5*(A(:,idx,:) + D(1:Nc(1),iy(idx),:));
end

C = reshape(real(A),Nc(1)*Nc(2),[]);

end



function mu_c = compute_mu_cos(c, eps, eta, Grid)
ck_full = expand_cos(c, Grid.N, Grid.Nc);
mu_full = compute_mu(ck_full, eps, eta, Grid);
mu_c = reduce_cos(mu_full, Grid.N, Grid.Nc);
end

function Hv = fch_jac_cos(u, vc, eps, eta, Grid)
v_full = expand_cos(vc, Grid.N, Grid.Nc);
Hv_full = fch_jac(u, v_full, eps, eta, Grid);
Hv = reduce_cos(Hv_full, Grid.N, Grid.Nc);
end

function J = fch_jac_matrix_cos(u, c, eps, eta, Grid)
% Dense (N/2+1)x(N/2+1) real Jacobian in the reduced cosine basis.

Nc = prod(Grid.Nc);
J = zeros(Nc, Nc);

% e = zeros(Nc, 1);

batch_size = 32;

for first = 1:batch_size:Nc

    last = min(first+batch_size-1,Nc);
    idx = first:last;

    V = zeros(Nc,numel(idx));
    V(sub2ind(size(V),idx,1:numel(idx))) = 1;

    J(:,idx) = fch_jac_cos_batch(u,V,eps,eta,Grid);

end

end

function J = fch_jac_matrix_cos_old(u, c, eps, eta, Grid)
% Dense (N/2+1)x(N/2+1) real Jacobian in the reduced cosine basis.

Nc = prod(Grid.Nc);
J = zeros(Nc, Nc);

e = zeros(Nc, 1);


for i = 1:Nc
    e(i) = 1;
    col = fch_jac_cos(u, e, eps, eta, Grid);
    J(:,i) = col(:);
    e(i) = 0;
end

end

function CK = expand_cos_batch(C, N, Nc)
% C: Nc(1)*Nc(2) x batch_size
% CK: N(1) x N(2) x batch_size

B = size(C,2);

C = reshape(C, Nc(1), Nc(2), B);

CK = zeros(N(1), N(2), B);

CK(1:Nc(1),1:Nc(2),:) = C;

CK(Nc(1)+1:N(1),1:Nc(2),:) = ...
    C(Nc(1)-1:-1:2,:,:);

CK(:,Nc(2)+1:N(2),:) = ...
    CK(:,Nc(2)-1:-1:2,:);

kx = 0:N(1)-1;
ky = (0:N(2)-1)';

phase_x = exp(2i*pi*kx*(N(1)/2-1)/N(1));
phase_y = exp(2i*pi*ky*(N(2)/2-1)/N(2));

phase = reshape(phase_y * phase_x, N(1), N(2), 1);

CK = CK .* conj(phase);

end

function HV = fch_jac_cos_batch(u, VC, eps, eta, Grid)
% Batched Jacobian-vector products in the reduced cosine basis.
%
% VC: Nc x B
% HV: Nc x B

N = prod(Grid.N);
Nx = Grid.N(1);
Ny = Grid.N(2);
Nc = Grid.Nc;
B = size(VC,2);

% Expand all directions simultaneously.
V = expand_cos_batch(VC, Grid.N, Nc);

k2 = reshape(Grid.k2, Nx, Ny, 1);

[~,dF,d2F] = SSAV_FCH_helpers.compute_F(u,0);
d3F = 6*u;

% Base-state chemical potential terms.
% what = -eps^2*Grid.k2.*fftn(u)/N;
what = -eps^2*Grid.k2.*fftn(u)/N - fftn(dF)/N;

% Base-state w.
w = real(ifft2(what*N));

% Batched perturbations in physical space.
dv = real(ifft2(N*V));

% First Jacobian operation.
dwhat = -eps^2*k2.*V ...
        - fft2(d2F.*dv)/N;

dw = real(ifft2(N*dwhat));

% Second Jacobian operation.
HV_full = -eps^2*k2.*dwhat ...
          - fft2(d2F.*dw)/N ...
          - fft2(d3F.*dv.*w)/N ...
          + eta*dwhat;

HV_full = -k2.*HV_full;

% Reduce back to independent cosine coefficients.
HV = reduce_cos_batch(HV_full, Grid.N, Nc);

end


function [c, resid_hist] = gradient_descent_cos(c, eps, eta, Grid, maxit, tol)
% Semi-implicit ("stabilized") gradient descent toward a genuine FCH
% steady state at fixed mass, starting from c.
%
% compute_mu_cos(c) is (proportional to) the mass-conserving L2 gradient
% of the energy. Its linear part is 6th-order in derivatives and very
% stiff (eps^4*k^6), so a naive explicit step would need an impractically
% tiny dtau for stability. Instead the linear part is treated implicitly
% -- cheap here since it's diag1onal in this cosine/Fourier basis -- while
% the nonlinear remainder is explicit, with adaptive step-size control:
% grow dtau on an accepted (energy-decreasing) step, halve and retry on
% any step that would increase the energy. The mean/mass mode needs no
% special handling: both the linear operator and compute_mu are
% structurally zero there (Grid.k2(1)=0), so c(1) never moves.

N = prod(Grid.N);
k2r = Grid.k2(1:Grid.Nc(1), 1:Grid.Nc(2));
k2r = k2r(:);
Lop = -eps^4*k2r.^3 + eta*eps^2*k2r.^2;

ck = expand_cos(c, Grid.N, Grid.Nc);
u = real(ifftn(ck * N));
E = SSAV_FCH_helpers.compute_E_direct(u, ck*N, eps, eta, Grid);

dtau = 1e-3;
dtau_min = 1e-8;
resid_hist = zeros(maxit, 1);

for it = 1:maxit
    mu = compute_mu_cos(c, eps, eta, Grid);
    resid_hist(it) = norm(mu(2:end));

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
    ck_trial = expand_cos(c_trial, Grid.N, Grid.Nc);
    u_trial = real(ifftn(ck_trial * N));
    E_trial = SSAV_FCH_helpers.compute_E_direct(u_trial, ck_trial*N, eps, eta, Grid);

    if E_trial > E
        dtau = dtau * 0.5;
        continue
    end

    c = c_trial;
    E = E_trial;
    dtau = dtau * 1.1;
end

end

function Jv = reduced_Jv(u, v, eps, eta, Grid)

vc = [0; v];
Jvc = fch_jac_cos(u, vc, eps, eta, Grid);
Jv = Jvc(2:end);

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
N = prod(Grid.N);
Nc = prod(Grid.Nc);
c(1) = m;
resid_hist = zeros(maxit, 1);

for it = 1:maxit
    mu = compute_mu_cos(c, eps, eta, Grid);
    resid_hist(it) = norm(mu(2:end));

    fprintf('Newton polish %3d: residual = %.3e\n', it, norm(mu(2:end)));

    if resid_hist(it) < tol
        resid_hist = resid_hist(1:it);
        return
    end

    ck = expand_cos(c, Grid.N, Grid.Nc);
    u = real(ifftn(ck * N));
    % % J = fch_jac_matrix_cos(u, c, eps, eta, Grid);
    % % A = J(2:end, 2:end);
    % % 
    % % dc = zeros(Nc, 1);
    % % dc(2:end) = J(2:end, 2:end) \ (-mu(2:end));
    % % 
    % % relres = norm(A*dc(2:end) + mu(2:end)) / max(norm(mu(2:end)), 1e-14);
    % % 
    % % fprintf('direct relres = %.3e\n', relres);
    % % 
    % % fprintf('  rcond(A) = %.3e, direct relres = %.3e\n', ...
    % %     rcond(A), relres);

    Jfun = @(v) reduced_Jv(u, v, eps, eta, Grid);

    k2 = Grid.k2(1:Grid.Nc(1), 1:Grid.Nc(2));
    k2 = k2(:);
    k2 = k2(2:end); 

    pdiag1 = eps^4*k2.^3 + eta*eps^2*k2.^2;

    % Avoid division by zero or excessively small entries.
    pmin = 1e-10 * max(1, max(abs(pdiag1)));
    pdiag1 = max(abs(pdiag1), pmin);

    Mfunc = @(v) v ./ pdiag1;
    
    [dc, flag, relres, iter] = gmres(Jfun, -mu(2:end), 50, 1e-8, 200, Mfunc);

    if flag ~= 0
        warning('Matrix-free GMRES did not converge');
    end

    % if flag ~= 0 || ~all(isfinite(dc))
    %     error('Matrix-free GMRES failed. Inspect the preconditioner and conditioning.');
    % end


    dc = [0; dc];

    c = c + dc;     % c(1) unchanged: dc(1) = 0, mass stays exactly m

end

end

function [data_out, u_vals_out, stability_out, diag1] = extend_branch(c, m, u, t, eps, eta, Grid, ...
    ds_start, ds_max, max_steps, m_bounds, loop_tol)
% Extends the FCH continuation branch from (c, m) along tangent t.
%
% Arclength is measured in the weighted norm
%   ||(c, m)||^2 = sum(Grid.wc .* c.^2) + m^2  =  mean((u - m)^2) + m^2
% (see Grid.wc in the main function), so ds is an RMS change of the profile.
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
%     4*loop_tol away from it first (so this can't fire immediately), or
%   - max_steps is reached, or
%   - ds underflows while retrying a rejected step (a robustness stop).
%
% stability_out is a 1xK vector aligned with data_out's columns: +1 where
% the point is linearly stable under the mass-conserving flow u_t = mu
% (every eigenvalue of the Jacobian, restricted to the mass-conserving
% subspace J(2:end,2:end), has negative real part), -1 otherwise. Just the
% sign is kept, not the full spectrum, to keep this cheap to store. See
% is_stable_cos below for how the eigenvalues are computed.
%
% diag1 is a struct with fields stop_reason (string), steps (number of
% accepted steps taken), and final_ds (ds at the point extension stopped).

N = prod(Grid.N);
Nc = prod(Grid.Nc);
ck = expand_cos(c, Grid.N, Grid.Nc);
E = SSAV_FCH_helpers.compute_E_direct(u, ck*N, eps, eta, Grid);
data_out = [m; E];
u_vals_out = cell(1, 1);
u_vals_out{1} = u;

J0 = fch_jac_matrix_cos(u, c, eps, eta, Grid);
stability_out = 2*is_stable_cos(J0, Grid) - 1;

% corrector settings (see the corrector loop below)
res_tol = 1e-8;         % converged when norm(b) < res_tol
newton_max = 12;        % iterations before the step is rejected
div_factor = 10;        % residual growth per iteration that counts as divergence
J_chord = J0;           % Jacobian at the last accepted point, reused by the corrector
e1 = zeros(1, Nc);      e1(1) = 1;

c0 = c;
m0 = m;
ds = ds_start;
ds_min = ds_start * 1e-2;
s_accum = 0;
left_start = false;     % set once the branch has moved 4*loop_tol away from the start
prev_disp = NaN;    % no reference displacement yet for the very first step

diag1.stop_reason = 'max_steps_reached';    % default; overridden below if
diag1.steps = max_steps;                    % something else stops us first
diag1.final_ds = ds;

for step = 1:max_steps

    accepted = false;
    while ~accepted

        % predictor
        c_curr = c + ds * t(1:end-1);
        m_curr = m + ds * t(end);
        ck_curr = expand_cos(c_curr, Grid.N, Grid.Nc);
        u_curr = real(ifftn(ck_curr * N));

        % corrector (Newton on the bordered system), with three safeguards:
        %   - A step counts as converged only when the RESIDUAL norm(b) is
        %     below res_tol. A small update alone is not enough: near a
        %     singular Jacobian the update can be small at a point that is
        %     not a solution.
        %   - It fails (-> rejected step and smaller ds, below) as soon as
        %     the residual or the update is non-finite, the residual grows
        %     by more than div_factor in one iteration, or newton_max
        %     iterations have been used.
        %   - The linear solve is direct, and the Jacobian is reused.
        %     Iterations start with J_chord, the dense Jacobian already
        %     built at the last accepted point (for its tangent and
        %     stability check), and it is rebuilt at the current iterate
        %     only when an iteration reduces the residual by less than 4x.
        %     The factorization of the bordered matrix is reused the same
        %     way, so most iterations cost one back-substitution.
        its = 0;
        success = false;
        bnorm_prev = inf;
        J = J_chord;
        Gdec = [];
        arc_row = (Grid.wc .* t(1:Nc))';    % arclength row: weighted tangent

        while true

            F_curr = compute_mu_cos(c_curr, eps, eta, Grid);

            b = [F_curr(2:end); c_curr(1) - m_curr;...
                arc_row * (c_curr - c) + (m_curr - m)*t(end) - ds];
            bnorm = norm(b);

            if ~isfinite(bnorm) || bnorm > div_factor * bnorm_prev
                break       % diverging: reject the step
            end

            if bnorm < res_tol
                success = true;
                break
            end

            if its >= newton_max
                break       % too slow: reject the step
            end

            % % % if its > 0 && bnorm > 0.25 * bnorm_prev
            % % %     % slow contraction with the reused Jacobian: refresh it
            % % %     J = fch_jac_matrix_cos(u_curr, c_curr, eps, eta, Grid);
            % % %     Gdec = [];
            % % % end
            % % % 
            % % % if isempty(Gdec)
            % % %     Gdec = decomposition([J(2:end,:), zeros(Nc-1, 1);
            % % %                           e1,          -1;
            % % %                           arc_row,     t(end)]);
            % % % % end
            % % % 
            % % % dx = Gdec \ (-b);

            Gfun = @(v) G_matrix_free( ...
                v, u_curr, eps, eta, Grid, Nc, t);

            % Check dimensions once before GMRES.
            vtest = zeros(Nc+1, 1);
            ytest = Gfun(vtest);
            assert(numel(ytest) == Nc+1, 'Incorrect Gfun output size.');

            % Reduced cosine modes, in the same column-major ordering as c.
            k2red = Grid.k2(1:Grid.Nc(1), 1:Grid.Nc(2));
            k2red = k2red(:);
            
            % Principal sixth-order FCH symbol.
            lambda = -eps^4*k2red.^3 + eta*eps^2*k2red.^2;
            
            % Regularize small eigenvalues while preserving their signs.
            scale = max(abs(lambda));
            floorval = max(1e-8*scale, 1e-12);
            
            lambda_reg = lambda;
            for j = 2:Nc
                if abs(lambda_reg(j)) < floorval
                    if lambda_reg(j) < 0
                        lambda_reg(j) = -floorval;
                    else
                        lambda_reg(j) = floorval;
                    end
                end
            end
            
            % Right preconditioner: approximate inverse of the bordered system.
            Pfun = @(r) apply_fch_prec(r, lambda_reg, t, Grid, Nc);
            
            % Solve G*P*y = -b, then recover dx = P*y.

            Afun = @(y) Gfun(Pfun(y));

            [y, flag, relres, iter] = gmres( ...
                Afun, -b, 50, 1e-7, 200);
            
            dx = Pfun(y);

            if flag ~= 0 || any(~isfinite(dx))
                fprintf('  GMRES failed: flag=%d, relres=%e\n', ...
                    flag, relres);
                break
            end

            if any(~isfinite(dx))
                break       % reject the step
            end

            alpha = 1;
            line_search_accepted = false;

            for ls = 1:12

                c_trial = c_curr + alpha*dx(1:end-1);
                m_trial = m_curr + alpha*dx(end);

                F_trial = compute_mu_cos(c_trial, eps, eta, Grid);

                % Force all vectors to columns
                F_trial = reshape(F_trial, [], 1);
                c_trial = reshape(c_trial, [], 1);
                c = reshape(c, [], 1);
                t = reshape(t, [], 1);
                wc = reshape(Grid.wc, [], 1);
                
                % Individual residuals
                b_physics = F_trial(2:end);
                b_mass = c_trial(1) - m_trial;
                
                b_arc = sum(wc .* t(1:Nc) .* (c_trial-c)) ...
                        + (m_trial-m)*t(end)-ds;
                
                % Explicit dimension checks
                assert(size(b_physics,1) == Nc-1);
                assert(size(b_physics,2) == 1);
                assert(isscalar(b_mass));
                assert(isscalar(b_arc));
                
                % Assemble
                b_trial = zeros(Nc+1,1);
                b_trial(1:Nc-1) = b_physics;
                b_trial(Nc) = b_mass;
                b_trial(Nc+1) = b_arc;

                if norm(b_trial) < norm(b)
                    line_search_accepted = true;
                    break
                end

                alpha = alpha / 2;
            end

            if ~line_search_accepted
                fprintf('  Line search failed at Newton iteration %d\n', its);
                break
            end

            % c_curr = c_curr + dx(1:end-1);
            % m_curr = m_curr + dx(end);
            c_curr = c_trial;
            m_curr = m_trial;
            u_curr = real(ifftn(expand_cos(c_curr, Grid.N, Grid.Nc) * N));

            bnorm_prev = bnorm;
            its = its + 1;

        end

        % Jacobian at the converged point: reused below for the tangent
        % update and the stability check, and as J_chord for the next step.
        J_final = [];
        if success
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
        this_disp = sqrt(sum(Grid.wc .* (c_curr - c).^2) + (m_curr - m)^2);
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
                diag1.stop_reason = 'ds_underflow';
                diag1.steps = step - 1;
                diag1.final_ds = ds;
                return
            end
        end
    end

    % commit the accepted step
    c = c_curr;
    ck = expand_cos(c, Grid.N, Grid.Nc);
    u = u_curr;
    m = m_curr;
    J_chord = J_final;
    s_accum = s_accum + ds;

    E = SSAV_FCH_helpers.compute_E_direct(u, ck * N, eps, eta, Grid);
    data_out = [data_out, [m; E]];
    u_vals_out{end+1} = u; %#ok<AGROW>

    % Stability at this point: reuse J_final (already built at this exact
    % (u, c) by the corrector above) instead of rebuilding it.
    is_stable = is_stable_cos(J_final, Grid);
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
        diag1.stop_reason = 'm_bound_reached';
        diag1.steps = step;
        diag1.final_ds = ds;
        break
    end

    % stop if the branch has closed into a loop back to the start:
    % same weighted norm as the steps, and armed once the branch has moved
    % 4*loop_tol away (a fixed arclength threshold before testing can be
    % longer than one lap, so a short loop ran around more than once)
    dist_to_start = sqrt(sum(Grid.wc .* (c - c0).^2) + (m - m0)^2);
    left_start = left_start || dist_to_start > 4*loop_tol;
    if left_start
        if dist_to_start < loop_tol
            fprintf('  extend_branch: closed loop detected at step %d (dist to start=%e, s=%.4f)\n', ...
                step, dist_to_start, s_accum);
            diag1.stop_reason = 'closed_loop';
            diag1.steps = step;
            diag1.final_ds = ds;
            break
        end
    end

    % update tangent -- reuse J_final (built at this exact (u, c) by the
    % corrector above) instead of rebuilding the Jacobian a third time.
    Jm = zeros(Nc-1, 1);
    Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
    A = [J_final(2:end,:), Jm; Fm; [Grid.wc.*t(1:end-1); t(end)]'];
    t = A \ [zeros(Nc, 1); 1];
    t = t / sqrt(sum(Grid.wc .* t(1:end-1).^2) + t(end)^2);

end

if strcmp(diag1.stop_reason, 'max_steps_reached')
    diag1.final_ds = ds;
end

end

function tf = is_stable_cos(J, Grid)
% Linear stability from the reduced Jacobian J: true when every eigenvalue
% of the mass-conserving block J(2:end,2:end) is negative.
%
% In the reduced cosine basis J = -K * W^-1 * H, with K = diag(k^2),
% W = diag(Grid.wc) (how many entries of the full ck each reduced
% coefficient stands for) and H the symmetric energy Hessian. J itself is
% not symmetric, but with D = sqrt(W/K) the similar matrix
%   S = D * J * D^-1 = -D^-1 * H * D^-1
% is. It has the same eigenvalues as J, they are exactly real, and eig
% uses its symmetric solver (about 2x faster at 64x64). The explicit
% 0.5*(S + S') only removes roundoff: S is symmetric to ~1e-17 already.

k2r = Grid.k2(1:Grid.Nc(1), 1:Grid.Nc(2));
k2r = k2r(:);

D = sqrt(Grid.wc(2:end) ./ k2r(2:end));

S = (D .* J(2:end, 2:end)) ./ D';
S = 0.5*(S + S');

tf = all(eig(S) < 0);

end

function y = G_matrix_free(v, u, eps, eta, Grid, Nc, t)

    dc = v(1:Nc);
    dm = v(end);

    % Physics equations
    q = fch_jac_cos(u, dc, eps, eta, Grid);
    y_physics = q(2:end);

    % Mass constraint
    y_mass = dc(1) - dm;

    % Weighted arclength constraint
    y_arc = sum(Grid.wc .* t(1:Nc) .* dc) ...
            + t(end)*dm;

    y = [y_physics; y_mass; y_arc];

end

function z = apply_block_prec(r, Pphys, Nc)
    z = zeros(size(r));
    z(1:Nc-1) = Pphys(r(1:Nc-1));
    z(Nc:end) = r(Nc:end);
end

function z = apply_fch_prec(r, lambda, t, Grid, Nc)

    r = r(:);
    z = zeros(Nc+1,1);

    wc = Grid.wc(:);

    % Approximate physics solve
    dc = zeros(Nc,1);

    dc(2:Nc) = r(1:Nc-1) ./ lambda(2:Nc);

    % Border residuals
    rmass = r(Nc);
    rarc  = r(Nc+1);

    % Weighted contribution from modes 2:Nc
    tangent_contribution = ...
        sum(wc(2:Nc) .* t(2:Nc) .* dc(2:Nc));

    % Solve the mass and arclength equations
    denom = wc(1)*t(1) + t(end);

    if abs(denom) < 1e-12
        error('Degenerate border preconditioner.');
    end

    dc(1) = (rarc - tangent_contribution ...
             + t(end)*rmass) / denom;

    dm = dc(1) - rmass;

    z(1:Nc) = dc;
    z(end) = dm;

end

