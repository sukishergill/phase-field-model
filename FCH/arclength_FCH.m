% arc-length continuation for FCH
%

eps = 0.1;      eta = eps^2;


% generate grid
dim = 1;
L = 2*pi;
N = 2^7;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

% initialize values
m = 0.001;         % initial m

% compute u, F(u), F'(u), F''(u) to compute initial E
% u = compute_u(Grid, m, eps);
u = Results.uu{end};
uhat = fft(u);
ck = uhat / N;

E = compute_E(u, uhat, eps, eta, Grid);       % initial E(u)

data = [m; E];

ds = 1/1000;

% initialize tangent vector
t = zeros(N+1, 1);      t(1) = 1;       t(end) = 1;

% t = initialize_t(m, 1/1000, ds, eps, eta, Grid);
% t = rand(N+1, 1);

t = t / norm(t, 2);         % initial tangent vector

while m < 1

% compute predictor
ck_curr = ck + ds * t(1:end-1);
m_curr = m + ds * t(end);
u_curr = ifft(ck_curr * N, 'symmetric');

% compute corrector
for j = 1:10

    grad = compute_RHS(gradientinit(ck_curr), u_curr, eps, eta, Grid);
    F_curr = grad.x;
    J = grad.dx;

    Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
    G = [J(2:end,:), zeros(N-1, 1); Fm; t'];
       
    b = [F_curr(2:end); ck_curr(1) - m;...
        (ck_curr - ck)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

    dx = G \ (-b);
    
    ck_new = ck_curr + dx(1:end-1);
    u_new = ifft(ck_new * N, 'symmetric');

    m_new = m_curr + dx(end);
    m_curr = m_new;

    % stopping criterion for Newton solver
    if abs(m_new - m_curr) < 1e-5
        break
    end

    ck_curr = ck_new;
    u_curr = u_new;

end

ck = ck_curr;
u = u_curr;
m = m_curr;

E = compute_E(u, ck * N, eps, eta, Grid);

data = [data, [m_curr; E]];

% update tangent
grad = compute_RHS(gradientinit(ck), u, eps, eta, Grid);
J = grad.dx;
Jm = zeros(N-1, 1);
Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
A = [J(2:end,:), Jm; Fm; t'];
t = A \ [zeros(N, 1); 1];   t = t / norm(t, 2);

end


figure;
plot(data(1,:), data(2,:), 'Linewidth', 3)
xlim([0, 1])



function u = compute_u(Grid, m, eps)

u = tanh((Grid.x + pi*(m + 1)/2) / eps) - ...
    tanh((Grid.x - pi*(m + 1)/2) / eps) - 1;

% u = m*ones(size(Grid.x));         % test case

end


function E = compute_E(u, uhat, eps, eta, Grid)

[F, dF, ~] = SSAV_FCH_helpers.compute_F(u, 0);

ux = real(ifftn(-1i*Grid.k.*uhat));
uxx = real(ifftn(-Grid.k.^2 .* uhat));

E = sum(0.5*(-eps^2*uxx + dF).^2 - eta*(eps^2*0.5*ux.^2 + F))*prod(Grid.d);

E = E / prod(Grid.L);

end


function lu = fch_coslap(u, L)
%FCH_COSLAP  Cosine-basis Laplacian via even extension and FFT.
%   lu = FCH_COSLAP(u,L) applies the Neumann/cosine Laplacian to a 1D
%   vector or 2D array sampled on the cell-centered grid. L is either a
%   scalar domain length in 1D, or [Lx Ly] in 2D. Defaults are 2*pi.

    if nargin < 2 || isempty(L)
        L = 2*pi;
    end
        s = size(u);
        u = u(:).';
        N = numel(u);
        Lx = L(1);

        ue  = [u, fliplr(u)];
        kap = (pi/Lx) * [0:N-1, -N:-1];
        le  = real(ifft(-(kap.^2) .* fft(ue)));
        lu  = reshape(le(1:N), s);
end
%------------------

% Action of the Jacobian on one vector
%------------------
function Hv = fch_jac(u, v, eps, eta, L, W, dW, d2W, d3W)
%FCH_JAC  Matrix-free FCH Hessian/Jacobian action.
%   Hv = FCH_JAC(u,v,eps,eta,L,W,dW,d2W,d3W) computes
%       Jac(u)v = (A+eta)A v - W'''(u) w v.

    if nargin < 5 || isempty(L),   L   = 2*pi;                    end
    if nargin < 6 || isempty(W),   W   = @(u) 0.25*(u.^2 - 1).^2; end
    if nargin < 7 || isempty(dW),  dW  = @(u) u.^3 - u;           end
    if nargin < 8 || isempty(d2W), d2W = @(u) 3*u.^2 - 1;         end
    if nargin < 9 || isempty(d3W), d3W = @(u) 6*u;                end

    lu = fch_coslap(u, L);
    w  = eps^2*lu - dW(u);
    A  = @(f) eps^2*fch_coslap(f, L) - d2W(u).*f;
    Av = A(v);
    Hv = A(Av) + eta*Av - d3W(u).*w.*v;
end

%------------------


% Build the full Jacobian
function J = fch_jac_matrix(u, eps, eta, L, varargin)
%FCH_JAC_MATRIX  Dense FCH Hessian/Jacobian matrix assembled by matvecs.
%   Intended for small 1D or small 2D problems.

    if nargin < 4 || isempty(L), L = 2*pi; end

    sz = size(u);
    n  = numel(u);
    J  = zeros(n, n);
    e  = zeros(sz);

    for j = 1:n
        e(j) = 1;
        col = fch_jac(u, e, eps, eta, L, varargin{:});
        J(:, j) = col(:);
        e(j) = 0;
    end
end
%------------------

function f = compute_RHS(ck, u, eps, eta, Grid)

f = eps^4 * Grid.k4 .* ck ...
    + eps^2 * Grid.k2 .* ck .* (fft(u.^3) - 1) ...
    + eps^2 * Grid.k2 .* ck.^1 .* (fft(3*u.^2 - 1)) ...
    - ck .* (fft(3*u.^5 - 4*u.^3) + 1) ...
    - eta * eps^2 * Grid.k2 .* ck ...
    - eta * ck .* (fft(u.^3) - 1);

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