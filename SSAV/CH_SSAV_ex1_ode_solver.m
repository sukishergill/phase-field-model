% Solve CH using an ode solver
model = 1;
epsilon = 0.1;

dim = 2;
L = [2*pi, 2*pi];
Nx = 64;       Ny = Nx;
Grid = SSAV_helpers.generate_Grid(L, [Nx, Ny], dim, model);
k2 = Grid.k2;       k4 = Grid.k4;

Nt = 101;
tf = 1;
t_plt = linspace(0, tf, Nt);

u_0 = 0.05*sin(Grid.xx).*sin(Grid.yy);

uu = cell(101,1);
uu{1} = u_0;

u_fft = fft2(u_0);
u_fft_vec = u_fft(:);



ode_func = @(t, u_fft_vec) compute_rhs(t, u_fft_vec, epsilon, k2, k4, Nx, Ny);

opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);

[t_vals, u_fft_sol] = ode15s(ode_func, [0 tf], u_fft_vec, opts);

for i = 2:size(t_vals,1)
    uu{i} = ifft2(reshape(u_fft_sol(i,:), Ny, Nx), 'symmetric');
end

function CH_rhs = compute_rhs(t, u_fft_vec, epsilon, k2, k4, Nx, Ny)

u_fft = reshape(u_fft_vec, Ny, Nx);
u = ifft2(u_fft, 'symmetric');

% [~, f] = SSAV_helpers.compute_F(u, 1);
f = u.^3 - u;

f_fft = fft2(f);

CH_rhs = -epsilon^2*k4.*u_fft - k2.*f_fft;
CH_rhs(1, 1) = 0;

CH_rhs = CH_rhs(:);

end