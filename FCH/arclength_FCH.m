% arc-length continuation for FCH
%

eps = 0.1;      eta = eps^2;


% generate grid
dim = 1;
L = 2*pi;
N = 2^7;
Grid = SSAV_FCH_helpers.generate_Grid(L, N, dim);

% initialize values
m = 0;         % initial m

u = compute_u(Grid, m, eps);
% u = Results.uu{end};
u_vals = cell(1001,1);
u_vals{1} = u;
uhat = fft(u);
ck = uhat / N;

E = compute_E(u, uhat, eps, eta, Grid);       % initial E(u)

data = [m; E];

ds = 1E-4;

% initialize tangent vector

[J, ~] = fch_jac_matrix(u, ck, eps, Grid);
Fm = zeros(Grid.N, 1);
Fm(1) = -1;

tu = -J\Fm;

t = [tu; 1];

% t = zeros(N+1, 1);      t(1) = 1;       t(end) = 1;

% t = initialize_t(m, ds, ds, eps, eta, Grid);
% t = rand(N+1, 1);
% t = [dck/ds; ds];

t = t / norm(t, 2);         % initial tangent vector

% i = 2;

for i = 1:1000

% compute predictor
ck_curr = ck + ds * t(1:end-1);
m_curr = m + ds * t(end);
u_curr = real(ifft(ck_curr * N));

%compute corrector
for j = 1:10

    ck_curr(1) = real(ck_curr(1));
    m_curr = real(m_curr);
    F_curr = compute_mu(ck_curr, eps, eta, Grid);

    [J, ~] = fch_jac_matrix(u_curr, eps, eta, Grid);

    Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
    G = [J(2:end,:), zeros(N-1, 1); Fm; t'];

    b = [F_curr(2:end); ck_curr(1) - m_curr;...
        (ck_curr - ck)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

    dx = G \ (-b);

    ck_new = ck_curr + dx(1:end-1);
    u_new = real(ifft(ck_new * N));
    % u_new = idct(ck_new * sqrt(N));

    m_new = m_curr + dx(end);
    m_curr = m_new;
    ck_curr = ck_new;
    u_curr = u_new;

    % stopping criterion for Newton solver
    if norm(dx) < 1e-8 && norm(b) < 1e-8
        break
    end


end

ck = ck_curr;
u = u_curr;
m = m_curr;

E = compute_E(u, ck * N, eps, eta, Grid);

u_vals{i} = u;
i = i + 1;
data = [data, [m_curr; E]];

% update tangent
[J, ~] = fch_jac_matrix(u_curr, eps, eta, Grid);
Jm = zeros(N-1, 1);
Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
A = [J(2:end,:), Jm; Fm; t'];
t = A \ [zeros(N, 1); 1];   t = t / norm(t, 2);

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


% function lu = fch_coslap(u, L)
% %FCH_COSLAP  Cosine-basis Laplacian via even extension and FFT.
% %   lu = FCH_COSLAP(u,L) applies the Neumann/cosine Laplacian to a 1D
% %   vector or 2D array sampled on the cell-centered grid. L is either a
% %   scalar domain length in 1D, or [Lx Ly] in 2D. Defaults are 2*pi.
% 
%     if nargin < 2 || isempty(L)
%         L = 2*pi;
%     end
%         s = size(u);
%         u = u(:).';
%         N = numel(u);
%         Lx = L(1);
% 
%         ue  = [u, fliplr(u)];
%         % kap = (pi/Lx) * [0:N-1, -N:-1];
%         kap = (pi/Lx) * [0:N -N+1:-1];
%         le  = real(ifft(-(kap.^2) .* fft(ue)));
%         lu  = reshape(le(1:N), s);
% end
% %------------------
% 
% % Action of the Jacobian on one vector
% %------------------
% function Hv = fch_jac(u, v, eps, eta, L, W, dW, d2W, d3W)
% %FCH_JAC  Matrix-free FCH Hessian/Jacobian action.
% %   Hv = FCH_JAC(u,v,eps,eta,L,W,dW,d2W,d3W) computes
% %       Jac(u)v = (A+eta)A v - W'''(u) w v.
% 
%     if nargin < 5 || isempty(L),   L   = 2*pi;                    end
%     if nargin < 6 || isempty(W),   W   = @(u) 0.25*(u.^2 - 1).^2; end
%     if nargin < 7 || isempty(dW),  dW  = @(u) u.^3 - u;           end
%     if nargin < 8 || isempty(d2W), d2W = @(u) 3*u.^2 - 1;         end
%     if nargin < 9 || isempty(d3W), d3W = @(u) 6*u;                end
% 
%     lu = fch_coslap(u, L);
%     w  = eps^2*lu - dW(u);
%     A  = @(f) eps^2*fch_coslap(f, L) - d2W(u).*f;
%     Av = A(v);
%     Hv = A(Av) + eta*Av - d3W(u).*w.*v;
% end
% 
% %------------------

function Hv = fch_jac(u, v, eps, eta, Grid)
% Action of Jacobian on one vector

[~, dF, d2F] = SSAV_FCH_helpers.compute_F(u, 0);

dv = real(ifft(Grid.N * v));

w = -eps^2*Grid.k2.*v - fft(dF.*dv)/Grid.N;

Hv = -eps^2*Grid.k2.*w - w.*fft(d2F.*dv)/Grid.N + eta*w;

end


% Build the full Jacobian
function [J, J_FD] = fch_jac_matrix(u, eps, eta, Grid)
%FCH_JAC_MATRIX  Dense FCH Hessian/Jacobian matrix assembled by matvecs.
%   Intended for small 1D or small 2D problems.

J = zeros(Grid.N);
J_FD = J;

e = zeros(size(u));
h = 1e-6;

for i = 1:Grid.N
    e(i) = 1;
    col = fch_jac(u, e, eps, Grid);
    J_FD(:,i) = (compute_mu(ck+e*h, eps, eta, Grid) - ...
        compute_mu(ck-e*h, eps, eta, Grid))/(2*h);
    J(:,i) = col(:);
    e(i) = 0;

end
end
%------------------

function mu = compute_mu(ck, eps, eta, Grid)

u = real(ifft(Grid.N*ck));
lap_ck = -Grid.k2 .* ck;
[~, dF, d2F] = SSAV_FCH_helpers.compute_F(u, 0);

mu = eps^4 * Grid.k4 .* ck ...
    - eps^2 * (-Grid.k2) .* (fft(dF)/Grid.N) ...
    - eps^2 * lap_ck .* (fft(d2F)/Grid.N) ...
    + fft(dF.*d2F)/Grid.N ...
    + eta * eps^2 * lap_ck ...
    - eta * fft(dF)/Grid.N;

end

% function f = compute_RHS(ck, u, eps, eta, Grid)
% N = 128;
% f = eps^4 * Grid.k4 .* ck ...
%     + eps^2 * Grid.k2 .* ck .* (fft(u.^3)/N - 1) ...
%     + eps^2 * Grid.k2 .* ck.^2 .* (fft(3*u.^2 - 1)/N) ...
%     - ck .* (fft(3*u.^5 - 4*u.^3)/N + 1) ...
%     - eta * eps^2 * Grid.k2 .* ck ...
%     - eta * ck .* (fft(u.^3)/N - 1);
% 
% end

% function f = compute_RHS_v2(uhat, u, eps, eta, Grid)
% N = 128;
% f = eps^4 * ifft(Grid.k4 .* uhat, 'symmetric') ...
%     - eps^2 * ifft(-Grid.k2 .* fft(u.^3 -u),'symmetric')...
%     - eps^2 * ifft(-Grid.k2 .* uhat,'symmetric').* (3*u.^2 - 1) ...
%     + 3*u.^5 - 4*u.^3 + u ...
%     + eta * eps^2 * ifft(-Grid.k2 .* uhat, 'symmetric')...
%     - eta * (u.^3 - u);
% 
% end

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

% [~, F1, F2] = SSAV_FCH_helpers.compute_F(uf, 0);
% mu = eps^4*ifft(Grid.k4.*uf_hat, 'symmetric') ...
%     -eps^2 * ifft(-Grid.k2 .* fft(F1), 'symmetric') ...
%     -eps^2 * F2 .* ifft(-Grid.k2 .* uf_hat, 'symmetric') ...
%     + F1 .* F2 ...
%     + eta * eps^2 * ifft(-Grid.k2 .* uf_hat, 'symmetric') ...
%     - eta * F1;

