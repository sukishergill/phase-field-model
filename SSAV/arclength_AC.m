% arclength continuation for Allen-Cahn
 
eps = 0.1;

% generate grid
dim = 1;
L = 2*pi;
N = 128;


Grid = SSAV_helpers.generate_Grid(L, N, dim);



m = 0;
u = compute_u(Grid, m, eps);
x = [u; m];
uhat = fft(u);
ck = uhat / N;

% f = @(ck) (-eps^2*ck .* Grid.k2 + ck);% ...
    %+ fft(ifft(ck) .* (ifft(ck * 128)).^3));

% E_vals = zeros(1,6);
E = compute_E(u, eps, Grid);

data = [m; E];

ds = 1/1024;

J = ac_jac_matrix(u, eps, L);

% initialize tangent vector
t = zeros(N+1, 1);      t(1) = 1;       t(end) = 1;
t = t / norm(t, 2);     % initial tangent vector;

for i = 1:1000

    % compute predictor
    ck_curr = ck + ds * t(1:end-1);
    m_curr = m + ds * t(end);
    % u_curr = u + ds * t(1:end-1);
    % ucurr_hat = fft(u_curr);
    % ucurr_hat = uhat + ds * t(1:end-1);
    % ucurr = ifft(ucurr_hat, 'symmetric');
    u_curr = ifft(ck_curr * N, 'symmetric');
    % % x0 = x + ds * t;
    % % x_curr = x0;
    % ck_curr = ck_0;
    % m_curr = m_0;

    for j = 1:10

        % u_curr = u0;
        % m_curr = m0;

        F_curr = compute_f(ck_curr, u_curr, eps, Grid);

        % J = ac_jac_matrix(u_curr, eps);
        grad = compute_f(gradientinit(ck), u_curr, eps, Grid);
        J = grad.dx;

        Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
        G = [J(2:end,:), zeros(N-1, 1); Fm; t'];
       
        b = [F_curr(2:end); ck_curr(1) - m;...
            (ck_curr - ck)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

        % b = [F_curr(2:end); ucurr_hat(1) - N*m_curr;...
        %     (ucurr_hat - uhat)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

        % cond(G)
        dx = G \ (-b);

        ck_new = ck_curr + dx(1:end-1);
        u_new = ifft(ck_new * N, 'symmetric');
        % unew_hat = ucurr_hat + dx(1:end-1);
        % u_new = ifft(unew_hat, 'symmetric');

        m_new = m_curr + dx(end);

        m_curr = m_new;
        if abs(m_new - m_curr) < 1e-5
            break
        end
        ck_curr = ck_new;
        u_curr = u_new;
        % ucurr_hat = unew_hat;
    end

    ck = ck_curr;
    u = u_curr;
    % uhat = ucurr_hat;
    m = m_curr;

    data = [data, [m_curr; ...
        compute_E(u, eps, Grid)]];

    % update tangent vector
    grad = compute_f(gradientinit(ck), u, eps, Grid);
    J = grad.dx;
    % J = ac_jac_matrix(u, eps);
    Jm = zeros(N-1, 1);
    Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
    A = [J(2:end,:), Jm; Fm; t'];
    t = A \ [zeros(N, 1); 1];   t = t / norm(t, 2);
end

figure;
plot(data(1,:), data(2,:), 'Linewidth', 3)


%------------------
function lu = ac_coslap(u, L)
%AC_COSLAP  Cosine-basis Laplacian via even extension and FFT.
%   lu = AC_COSLAP(u,L) applies the Neumann/cosine Laplacian to a 1D
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
function Hv = ac_jac(u, v, eps, L, W, dW, d2W)
%AC_JAC  Matrix-free AC Hessian/Jacobian action.
%   Hv = AC_JAC(u,v,eps,eta,L,W,dW,d2W,d3W) computes
%       Jac(u)v = (A+eta)A v - W'''(u) w v.

    if nargin < 4 || isempty(L),   L   = 2*pi;                    end
    if nargin < 5 || isempty(W),   W   = @(u) 0.25*(1 - u.^2).^2; end
    if nargin < 6 || isempty(dW),  dW  = @(u) u - u.^3;           end
    if nargin < 7 || isempty(d2W), d2W = @(u) 1 - 3*u.^2;         end

    % lu = fch_coslap(u, L);
    % w  = eps^2*lu - dW(u);
    % A  = @(f) eps^2*fch_coslap(f, L) - d2W(u).*f;
    A = @(f) -eps^2*ac_coslap(f, L) - dW(u);
    Av = A(v);
    % Hv = A(Av) + eta*Av - d3W(u).*w.*v;
    Hv = Av - d2W(u).*v;
end

%------------------


% Build the full Jacobian
function J = ac_jac_matrix(u, eps, L, varargin)
%AC_JAC_MATRIX  Dense ac Hessian/Jacobian matrix assembled by matvecs.
%   Intended for small 1D or small 2D problems.

    if nargin < 3 || isempty(L), L = 2*pi; end

    sz = size(u);
    n  = numel(u);
    J  = zeros(n, n);
    e  = zeros(sz);

    for j = 1:n
        e(j) = 1;
        col = ac_jac(u, e, eps, L, varargin{:});
        J(:, j) = col(:);
        e(j) = 0;
    end
end
%------------------

function u = compute_u(Grid, m, eps)


% u = tanh((Grid.xx + pi*(m + 1)/2) / eps) - ...
%     tanh((Grid.xx - pi*(m + 1)/2) / eps) - 1;
u = m*ones(size(Grid.xx));         % test case

end

function Eu = compute_E(u, eps, Grid)

uhat = fft(u);
ux = real(ifftn(-1i*Grid.k.*uhat));
F = 0.25 * (u.^2 - 1).^2;

Eu = sum(0.5*eps^2*ux.^2 + F)*prod(Grid.d);
Eu = Eu / prod(Grid.L);

end

function F = compute_F(ck, u, eps, Grid)

uhat = fft(u);
F = -eps^2*ck .* Grid.k2 + ck - ck .* fft(u.^3);

% F = eps^2 * ifft(-Grid.k2.*uhat, 'symmetric') + u.^3 - u;

end

function f = compute_f(ck, u, eps, Grid)

f = -eps^2*ck .* Grid.k2 - ck + ck.* fft(u.^3);

end

