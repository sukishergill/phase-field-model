% arclength continuation for Allen-Cahn


eps = 0.1;

% generate grid
dim = 1;
L = 2*pi;
N = 128;
Grid = SSAV_helpers.generate_Grid(L, N, dim);


m = 0;
% u = compute_u(Grid, m, eps);
% AC_SSAV_ex1_1D;
u = Results.uu{end};
u_vals = cell(1001,1);
u_vals{1} = u;
uhat = fft(u);
ck = uhat / N;

% E_vals = zeros(1,6);
E = compute_E(u, eps, Grid);

data = [m; E];

ds = 1/1000;

[~, J] = ac_jac_matrix(u, ck, eps, Grid);
Fm = zeros(Grid.N, 1);
Fm(1) = -1;

tu = -J\Fm;

t = [tu; 1];

% initialize tangent vector
% t = zeros(N+1, 1);      t(1) = 1;       t(end) = 1;
% t = initialize_t(m, 1/1000, ds, eps, Grid);
t = t / norm(t, 2);     % initial tangent vector;

for i = 1:1000

    % compute predictor
    ck_curr = ck + ds * t(1:end-1);
    m_curr = m + ds * t(end);
    m_curr = real(m_curr);
    u_curr = real(ifft(ck_curr * N));

    for j = 1:10

        F_curr = compute_F(ck_curr, eps, Grid);

        [J, ~] = ac_jac_matrix(u_curr, ck_curr, eps, Grid);

        Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
        G = [J(2:end,:), zeros(N-1, 1); Fm; t'];
       
        b = [F_curr(2:end); ck_curr(1) - m_curr;...
            (ck_curr - ck)'*t(1:end-1) + (m_curr - m)*t(end) - ds];

        dx = G \ (-b);

        ck_new = ck_curr + dx(1:end-1);
        u_new = real(ifft(ck_new * N));

        m_new = m_curr + real(dx(end));

        m_curr = m_new;
        ck_curr = ck_new;
        u_curr = u_new;
        if norm(dx) < 1e-8 && norm(b) < 1e-8
            break
        end

    end

    ck = ck_curr;
    u = u_curr;
    m = m_curr;
    u_vals{i+1} = u;
    data = [data, [m_curr; ...
        compute_E(u, eps, Grid)]];

    % update tangent vector
    [J, ~] = ac_jac_matrix(u, ck, eps, Grid);
    Jm = zeros(N-1, 1);
    Fm = zeros(1, N+1);     Fm(1) = 1; Fm(end) = -1;
    A = [J(2:end,:), Jm; Fm; t'];
    t = A \ [zeros(N, 1); 1];   t = t / norm(t, 2);
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
