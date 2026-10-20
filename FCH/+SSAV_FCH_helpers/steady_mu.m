function mu = steady_mu(u, Para, Grid)
% Chemical potential of the FCH energy used by FCH_BDF2_SAV_adap_dt,
%   E(u)  = int 1/2 (eps^2 Lap u - W'(u))^2 - eta1 eps^2/2 |grad u|^2 - eta2 W(u)
%   mu(u) = (eps^2 Lap - W''(u) + eta1)(eps^2 Lap u - W'(u)) + (eta1 - eta2) W'(u)
% with W' = F1, W'' = F2 from compute_F(u, Para.tau). Steady states of
% u_t = Lap mu at fixed mass solve mu = const.
[~, W1, W2] = SSAV_FCH_helpers.compute_F(u, Para.tau);
eps2 = Para.epsilon^2;
lap = @(v) SSAV_FCH_helpers.spectral_lap(v, Grid);
g = eps2*lap(u) - W1;
mu = eps2*lap(g) - W2.*g + Para.eta1*g + (Para.eta1 - Para.eta2)*W1;
end
