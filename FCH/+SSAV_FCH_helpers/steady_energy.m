function E = steady_energy(u, Para, Grid)
% Energy per unit volume consistent with steady_mu (|grad u|^2 integrated
% by parts to -u Lap u, exact for the spectral Laplacian). Differs from
% compute_E_direct at the ~1e-9 level because that zeroes the Nyquist mode
% of its derivatives.
[W, W1] = SSAV_FCH_helpers.compute_F(u, Para.tau);
eps2 = Para.epsilon^2;
lap = @(v) SSAV_FCH_helpers.spectral_lap(v, Grid);
g = eps2*lap(u) - W1;
E = mean(0.5*g.^2 + 0.5*Para.eta1*eps2*u.*lap(u) - Para.eta2*W, 'all');
end
