% arclength continuation for Allen-Cahn
%
% Uses a cosine-only (even-symmetric) representation: u is assumed even
% about the domain center, so its Fourier coefficients ck are real and
% palindromic (ck(N+2-j) = ck(j)). The state is stored as the reduced,
% real vector c = ck(1:N/2+1) of independent cosine amplitudes rather than
% the full complex ck. This eliminates the translation mode by
% construction -- du/dx is odd, so translation simply isn't a direction
% available within the even/cosine subspace -- instead of relying on
% bordering to paper over it after the fact.

eps = 0.1;

% generate grid
dim = 1;
L = 2*pi;
N = 128;
Grid = SSAV_helpers.generate_Grid(L, N, dim);
Nc = N/2 + 1;       % number of independent cosine (even) Fourier coefficients


m = 0;
% u = compute_u(Grid, m, eps);
% AC_SSAV_ex1_1D;
u = Results.uu{end};
u_vals = cell(1001,1);
u_vals{1} = u;
ck = fft(u) / N;
c = reduce_cos(ck, N);      % restrict to the cosine (even) subspace

E = compute_E(u, eps, Grid);

data = [m; E];

ds = 1/1000;

[J, J_FD] = ac_jac_matrix_cos(u, c, eps, Grid);

% The raw (un-bordered) reduced Jacobian J is still singular in the k=0
% row for the same reason as before (AC deliberately treats the mean-mode
% equation as a free Lagrange-multiplier slot so m can be imposed
% directly), so we still border with the linearized mass constraint
% c(1) - m = 0 and take the resulting system's null vector as the initial
% tangent -- but there is no longer a translation mode to also worry
% about, since that direction doesn't exist in this reduced parameterization.
Jm = zeros(Nc-1, 1);
Fm_row = zeros(1, Nc+1);    Fm_row(1) = 1; Fm_row(end) = -1;
B = [J(2:end,:), Jm; Fm_row];

% null(B) relies on a rank tolerance that scales with the largest singular
% value, which can lump more than one small-but-distinct singular value in
% as "null" for an ill-conditioned matrix. Sidestep the ambiguity: take the
% right singular vector for the single smallest singular value directly.
[~, ~, V] = svd(B);
t = V(:, end);
if t(end) < 0
    t = -t;
end

% initialize tangent vector
% t = zeros(Nc+1, 1);      t(1) = 1;       t(end) = 1;
% t = initialize_t(m, 1/1000, ds, eps, Grid);
t = t / norm(t, 2);     % initial tangent vector;

for i = 1:1000

    % compute predictor
    c_curr = c + ds * t(1:end-1);
    m_curr = m + ds * t(end);
    ck_curr = expand_cos(c_curr, N);
    u_curr = real(ifft(ck_curr * N));

    for j = 1:10

        F_curr = compute_F_cos(c_curr, eps, Grid);

        [J, ~] = ac_jac_matrix_cos(u_curr, c_curr, eps, Grid);

        Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
        G = [J(2:end,:), zeros(Nc-1, 1); Fm; t'];

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
        if norm(dx) < 1e-8 && norm(b) < 1e-8
            break
        end

    end

    c = c_curr;
    ck = expand_cos(c, N);
    u = u_curr;
    m = m_curr;
    u_vals{i+1} = u;
    data = [data, [m_curr; ...
        compute_E(u, eps, Grid)]];

    % update tangent vector
    [J, ~] = ac_jac_matrix_cos(u, c, eps, Grid);
    Jm = zeros(Nc-1, 1);
    Fm = zeros(1, Nc+1);     Fm(1) = 1; Fm(end) = -1;
    A = [J(2:end,:), Jm; Fm; t'];
    t = A \ [zeros(Nc, 1); 1];   t = t / norm(t, 2);
end

figure;
plot(data(1,:), data(2,:), 'Linewidth', 3)
xlabel('$m$', 'Interpreter','latex')
ylabel('$E$', 'Interpreter', 'latex')
set(gca, 'FontSize', 30)
set(gca, 'TickLabelInterpreter', 'latex')


function u = compute_u(Grid, m, eps)


u = tanh((Grid.xx + pi*(m + 1)/2) / eps) - ...
    tanh((Grid.xx - pi*(m + 1)/2) / eps) - 1;
% u = m*ones(size(Grid.xx));         % test case

end

function Eu = compute_E(u, eps, Grid)
uhat = fft(u);
ux = real(ifftn(-1i*Grid.k.*uhat));
F = 0.25 * (u.^2 - 1).^2;


Eu = sum(0.5*eps^2*ux.^2 + F)*prod(Grid.d);
Eu = Eu / prod(Grid.L);

end

function F = compute_F(ck, eps, Grid)

u = real(ifft(Grid.N*ck));
F = -eps^2*ck .* Grid.k2 + ck - fft(u.^3)/Grid.N;

end

function t = initialize_t(m, dm, ds, eps, Grid)

up = compute_u(Grid, m + dm, eps);
ckp = fft(up) / Grid.N;
Fp = compute_F(ckp, eps, Grid);

um = compute_u(Grid, m - dm, eps);
ckm = fft(um) / Grid.N;
Fm = compute_F(ckm, eps, Grid);

dFdm = (Fp - Fm) / (2*dm);

% dmds = dm / ds;

t = [dFdm; 1];
end


function Hv = ac_jac(u, v, eps, Grid)
% Action of Jacobian on one vector

dv = real(ifft(Grid.N * v));

Hv = -eps^2*Grid.k2.*v + v - fft(3*u.^2.*dv)/Grid.N;

end

function [J, J_FD] = ac_jac_matrix(u, ck, eps, Grid)
% Constuct the dense Jacobian matrix in Fourier space

J = zeros(Grid.N);
J_FD = J;

e = zeros(size(u));
h = 1e-6;

for i = 1:Grid.N
    e(i) = 1;
    col = ac_jac(u, e, eps, Grid);
    J_FD(:,i) = (compute_F(ck+e*h, eps, Grid) - ...
        compute_F(ck-e*h, eps, Grid))/(2*h);
    J(:,i) = col(:);
    e(i) = 0;

end

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
% (length N/2+1). Grid.xx is centered at x=0 but is a 1-indexed grid
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

function Fc = compute_F_cos(c, eps, Grid)
ck_full = expand_cos(c, Grid.N);
F_full = compute_F(ck_full, eps, Grid);
Fc = reduce_cos(F_full, Grid.N);
end

function Hv = ac_jac_cos(u, vc, eps, Grid)
v_full = expand_cos(vc, Grid.N);
Hv_full = ac_jac(u, v_full, eps, Grid);
Hv = reduce_cos(Hv_full, Grid.N);
end

function [J, J_FD] = ac_jac_matrix_cos(u, c, eps, Grid)
% Dense (N/2+1)x(N/2+1) real Jacobian in the reduced cosine basis.

Nc = Grid.N/2 + 1;
J = zeros(Nc);
J_FD = J;

e = zeros(Nc, 1);
h = 1e-6;

for i = 1:Nc
    e(i) = 1;
    col = ac_jac_cos(u, e, eps, Grid);
    J_FD(:,i) = (compute_F_cos(c+e*h, eps, Grid) - ...
        compute_F_cos(c-e*h, eps, Grid))/(2*h);
    J(:,i) = col(:);
    e(i) = 0;
end

end
