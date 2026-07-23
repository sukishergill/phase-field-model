% arc-length continuation for FCH
%
% Uses a cosine-only (even-symmetric) representation: u is assumed even
% about the domain center, so its Fourier coefficients ck are real and
% palindromic (ck(N+2-j) = ck(j)). The state is stored as the reduced,
% real vector c = ck(1:N/2+1) of independent cosine amplitudes rather than
% the full complex ck. This eliminates the translation mode by
% construction -- du/dx is odd, so translation simply isn't a direction
% available within the even/cosine subspace -- instead of relying on
% bordering to paper over it after the fact.

eps = 0.1;      eta = eps^2;


% generate grid
dim = 1;
L = 2*pi;
N = 2^7;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);
Nc = N/2 + 1;       % number of independent cosine (even) Fourier coefficients

% initialize values
m = 0.01;         % initial m

% u = compute_u(Grid, m, eps);
u = Results.uu{end};
u_vals = cell(1001,1);
u_vals{1} = u;
ck = fft(u) / N;
c = reduce_cos(ck, N);      % restrict to the cosine (even) subspace

E = compute_E(u, ck*N, eps, eta, Grid);       % initial E(u)

data = [m; E];

ds = 1E-4;

% initialize tangent vector

[J, J_FD] = fch_jac_matrix_cos(u, c, eps, eta, Grid);

% The raw (un-bordered) reduced Jacobian J is still singular in row 1: the
% mean/k=0 mode is identically zero because fch_jac/compute_mu end by
% multiplying through by -Grid.k2, and Grid.k2(1)=0. Solving J\Fm directly
% would still trigger a singular-matrix error. Instead build the same
% bordered (under-determined) system used later in the loop -- drop row 1
% and replace it with the linearized mass constraint tu(1) - tm = 0 -- and
% take its null vector directly as the initial tangent. There is no longer
% a separate translation mode to worry about, since that direction doesn't
% exist in this reduced (even-only) parameterization.
Jm = zeros(Nc-1, 1);
Fm_row = zeros(1, Nc+1);    Fm_row(1) = 1; Fm_row(end) = -1;
B = [J(2:end,:), Jm; Fm_row];

% null(B) relies on a rank tolerance that scales with the largest singular
% value; B's singular values here span many orders of magnitude (the FCH
% operator is 6th-order), so several small-but-distinct singular values can
% end up under that same tolerance and null() returns more than one column.
% Sidestep the ambiguity: take the right singular vector for the single
% smallest singular value directly.
[~, ~, V] = svd(B);
t = V(:, end);
if t(end) < 0
    t = -t;
end


% t = zeros(Nc+1, 1);      t(1) = 1;       t(end) = 1;

% t = initialize_t(m, ds, ds, eps, eta, Grid);
% t = rand(Nc+1, 1);
% t = [dck/ds; ds];

t = t / norm(t, 2);         % initial tangent vector

% i = 2;

for i = 1:1000

% compute predictor
c_curr = c + ds * t(1:end-1);
m_curr = m + ds * t(end);
ck_curr = expand_cos(c_curr, N);
u_curr = real(ifft(ck_curr * N));

%compute corrector
for j = 1:10

    F_curr = compute_mu_cos(c_curr, eps, eta, Grid);

    [J, ~] = fch_jac_matrix_cos(u_curr, c_curr, eps, eta, Grid);

    Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
    Jm = zeros(Nc-1, 1);

    G = [J(2:end,:), Jm; Fm; t'];

    b = [F_curr(2:end); c_curr(1) - m_curr;...
        (c_curr - c)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

    dx = G \ (-b);

    c_new = c_curr + dx(1:end-1);
    ck_new = expand_cos(c_new, N);
    u_new = real(ifft(ck_new * N));

    m_new = m_curr + dx(end);
    m_curr = m_new;
    c_curr = c_new;
    u_curr = u_new;

    % stopping criterion for Newton solver
    if norm(dx) < 1e-8 && norm(b) < 1e-8
        break
    end


end

c = c_curr;
ck = expand_cos(c, N);
u = u_curr;
m = m_curr;

E = compute_E(u, ck * N, eps, eta, Grid);

u_vals{i+1} = u;
data = [data, [m; E]];

% update tangent
[J, ~] = fch_jac_matrix_cos(u, c, eps, eta, Grid);
Jm = zeros(Nc-1, 1);
Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
A = [J(2:end,:), Jm; Fm; t'];
t = A \ [zeros(Nc, 1); 1];   t = t / norm(t, 2);

end


figure;
plot(data(1,:), data(2,:), 'Linewidth', 3)
% xlim([0, 1])
xlabel('$m$', 'Interpreter','latex')
ylabel('$E$', 'Interpreter', 'latex')
set(gca, 'FontSize', 30)
set(gca, 'TickLabelInterpreter', 'latex')



function u = compute_u(Grid, m, eps)

u = tanh((Grid.x + pi*(m + 1)/2) / eps) - ...
    tanh((Grid.x- pi*(m + 1)/2) / eps) - 1;

% u = m*ones(size(Grid.x));         % test case

end


function E = compute_E(u, uhat, eps, eta, Grid)

[F, dF, ~] = SSAV_FCH_helpers.compute_F(u, 0);

ux = real(ifftn(-1i*Grid.k.*uhat));
uxx = real(ifftn(-Grid.k.^2 .* uhat));

E = sum(0.5*(-eps^2*uxx + dF).^2 - eta*(eps^2*0.5*ux.^2 + F))*prod(Grid.d);

E = E / prod(Grid.L);

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

function t = initialize_t(m, dm, ds, eps, eta, Grid)

up = compute_u(Grid, m + dm, eps);
ckp = fft(up) / Grid.N;
Fp = compute_RHS(ckp, up, eps, eta, Grid);

um = compute_u(Grid, m - dm, eps);
ckm = fft(up) / Grid.N;
Fm = compute_RHS(ckm, um, eps, eta, Grid);

dFdm = 2*(Fp - Fm) / dm;

% dmds = dm / ds;

t = [dFdm; dm];
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
