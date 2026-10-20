function w = steady_Jv(u, v, Para, Grid)
% Exact Jacobian of steady_mu at u applied to v (W''' = 6u for any tau).
% J is the Hessian of the energy: symmetric. v may be a stack of fields
% (N x B in 1D, Ny x Nx x B in 2D) to apply J to B vectors at once.
[~, W1, W2] = SSAV_FCH_helpers.compute_F(u, Para.tau);
eps2 = Para.epsilon^2;
lap = @(v) SSAV_FCH_helpers.spectral_lap(v, Grid);
g = eps2*lap(u) - W1;
dg = eps2*lap(v) - W2.*v;
w = eps2*lap(dg) - W2.*dg + Para.eta1*dg - 6*u.*v.*g + (Para.eta1 - Para.eta2)*W2.*v;
end
