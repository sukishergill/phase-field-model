% Checks for newton_polish_FCH. Prints PASS/FAIL per check.
%   1. mu is the gradient of the energy the time-stepper uses (compute_E)
%   2. Jv matches finite differences of mu, and J is symmetric
%   3. dense J (16^2): finite-difference columns match Jv columns
%   4. polishing a final state from the m-sweep (Ex5_test.mat)
%   5. a shifted, non-even copy of that state polishes to the shifted solution
%   6. stripes: one translation direction; constant state: none

clear
pass = @(ok, name) fprintf('%s  %s\n', ternary(ok, 'PASS', 'FAIL'), name);

Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = 0.7*Para.epsilon^2;     % eta1 ~= eta2 and tau ~= 0 exercise every term
Para.tau = 0.3;
Para.B = 1e3;

L = [2*pi, 2*pi];
N = [32, 32];
Grid = SSAV_FCH_helpers.generate_Grid(L, N, 2);

% smooth random field with no Nyquist content (compute_E zeroes the Nyquist
% mode of its Laplacian, the solver does not -- irrelevant for smooth u)
rng(1);
smooth = @(G) real(ifft2(fft2(randn(size(G.xx))) .* exp(-G.k2/8)));
u = 0.2 + 0.8*smooth(Grid)/max(abs(smooth(Grid)), [], 'all');
v = smooth(Grid);   v = v - mean(v(:));
w = smooth(Grid);   w = w - mean(w(:));

[~, info] = newton_polish_FCH(u, Grid, Para, struct('maxit', 0, 'verbose', false));

%% 1. mu = dE/du for the time-stepper's energy
% tau = 0 here: for tau ~= 0, compute_F's F1 is not dF/du, so no mu can be
% the exact gradient of compute_E + compute_G. (eta1 ~= eta2 still tested.)
ParaE = Para;   ParaE.tau = 0;
[~, infoE] = newton_polish_FCH(u, Grid, ParaE, struct('maxit', 0, 'verbose', false));
E_fun = @(u) energy_SAV(u, ParaE, Grid);
h = 1e-5;
dE_fd = (E_fun(u + h*v) - E_fun(u - h*v)) / (2*h);
dE_mu = mean(infoE.mu(u) .* v, 'all');          % compute_E is normalized by |Omega|
err1 = abs(dE_fd - dE_mu) / abs(dE_mu);
pass(err1 < 1e-6, sprintf('mu vs dE/du (compute_E + compute_G): rel err %.1e', err1));

%% 2. Jv vs finite differences of mu; symmetry
Jv_fd = (info.mu(u + h*v) - info.mu(u - h*v)) / (2*h);
err2 = norm(info.Jv(u, v) - Jv_fd, 'fro') / norm(Jv_fd, 'fro');
pass(err2 < 1e-6, sprintf('Jv vs finite-difference Jacobian: rel err %.1e', err2));

sym_err = abs(sum(w .* info.Jv(u, v), 'all') - sum(v .* info.Jv(u, w), 'all')) / ...
    abs(sum(w .* info.Jv(u, v), 'all'));
pass(sym_err < 1e-10, sprintf('J symmetric: rel err %.1e', sym_err));

%% 3. dense Jacobian at 16^2, column by column
G16 = SSAV_FCH_helpers.generate_Grid(L, [16 16], 2);
u16 = 0.2 + 0.5*cos(G16.xx).*cos(2*G16.yy) + 0.1*sin(G16.yy);
n16 = numel(u16);
J_an = zeros(n16);  J_fd = zeros(n16);
[~, info16] = newton_polish_FCH(u16, G16, Para, struct('maxit', 0, 'verbose', false));
for j = 1:n16
    e = zeros(size(u16));  e(j) = 1;
    J_an(:, j) = reshape(info16.Jv(u16, e), [], 1);
    J_fd(:, j) = reshape((info16.mu(u16 + h*e) - info16.mu(u16 - h*e)) / (2*h), [], 1);
end
err3 = norm(J_an - J_fd, 'fro') / norm(J_an, 'fro');
pass(err3 < 1e-6, sprintf('dense J (16^2): rel err %.1e, asymmetry %.1e', err3, ...
    norm(J_an - J_an', 'fro')/norm(J_an, 'fro')));

%% 4. polish a real final state from the m-sweep
S = load('Ex5_test.mat', 'final_u', 'Grids', 'm_vals');
ParaS.epsilon = 0.18;
ParaS.eta1 = ParaS.epsilon^2;
ParaS.eta2 = ParaS.eta1;
ParaS.tau = 0;

i_test = find(abs(S.m_vals - 0.3) < 1e-12);
u0 = S.final_u{i_test};
G64 = S.Grids{i_test};
fprintf('\n--- sweep state m = %.1f ---\n', S.m_vals(i_test));
[u_pol, inf4] = newton_polish_FCH(u0, G64, ParaS);
[ny, nx] = size(u_pol);
ix = [nx-1:-1:1, nx];  iy = [ny-1:-1:1, ny];
even_err = max(abs(u_pol - u_pol(:, ix)), [], 'all') + max(abs(u_pol - u_pol(iy, :)), [], 'all');
pass(inf4.converged && ~inf4.collapsed, 'sweep state converged, not collapsed');
pass(even_err < 1e-8, sprintf('polished state still even: err %.1e', even_err));
pass(norm(inf4.sigma(:, end))/sqrt(numel(u0)) < 1e-8, 'sigma -> 0');
pass(abs(mean(u_pol(:)) - mean(u0(:))) < 1e-12, 'mass preserved');

%% 5. shifted, non-even copy
shift = [7, 19];
fprintf('\n--- same state shifted by %s grid points (no longer even) ---\n', mat2str(shift));
[u_sh, inf5] = newton_polish_FCH(circshift(u0, shift), G64, ParaS);
sh_err = max(abs(u_sh - circshift(u_pol, shift)), [], 'all');
pass(inf5.converged, 'shifted state converged');
pass(sh_err < 1e-7, sprintf('shifted result = shifted polished state: max err %.1e', sh_err));
pass(abs(inf5.E - inf4.E) < 1e-10, sprintf('same energy: diff %.1e', inf5.E - inf4.E));

%% 6. stripes and constant state
% 1D steady state (also exercises the 1D path), then extend it in y
G1 = SSAV_FCH_helpers.generate_Grid(2*pi, 64, 1);
u1 = tanh((G1.x + 1)/0.18) - tanh((G1.x - 1)/0.18) - 1;
u1 = circshift(u1, 5);      % not even, so the 1D phase condition is active
fprintf('\n--- 1D bilayer ---\n');
[u1_pol, inf1] = newton_polish_FCH(u1, G1, ParaS);
pass(inf1.converged && size(inf1.T, 2) == 1, '1D: converged with 1 translation direction');

G2 = SSAV_FCH_helpers.generate_Grid([2*pi 2*pi], [64 64], 2);
u_str = repmat(u1_pol', 64, 1) + 1e-6*randn(64);     % x-stripes, slightly perturbed
fprintf('\n--- 2D stripes ---\n');
[~, inf6] = newton_polish_FCH(u_str, G2, ParaS);
pass(inf6.converged && size(inf6.T, 3) == 1, sprintf('stripes: converged with %d translation direction', size(inf6.T, 3)));

[~, inf7] = newton_polish_FCH(0.3*ones(64), G2, ParaS, struct('verbose', false));
pass(isempty(inf7.T) || size(inf7.T, 3) == 0, 'constant state: no translation directions');

%% helpers
function E = energy_SAV(u, Para, Grid)
% physical energy as the time-stepper computes it: compute_E with w^2 - B = G
[F, F1] = SSAV_FCH_helpers.compute_F(u, Para.tau);
G = SSAV_FCH_helpers.compute_G(u, Para, F, F1, Grid);
w = SSAV_FCH_helpers.compute_w(G, Para.B);
E = SSAV_FCH_helpers.compute_E(fft2(u), w, Para, Grid);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
