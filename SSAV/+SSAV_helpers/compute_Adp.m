function [A_dp, rel_err, F_tilde_new, num_fft] = ...
    compute_Adp(u_new, u, u_new_fft, F_tilde_old, F_tilde, dt_new, ...
    gamma, Para, Grid, num_fft)

[~, f] = SSAV_helpers.compute_F(u_new, Para.beta);

[F_tilde_new, num_fft] = SSAV_helpers.compute_F_tilde(u_new, ...
    u_new_fft, f, Grid, Para, num_fft);

u_p = SSAV_helpers.AM3_up(u, gamma, dt_new, F_tilde_old, ...
    F_tilde, F_tilde_new);

rel_err = sum((u_new - u_p).^2, 'all') / sum((u_p).^2, 'all');

A_dp = Para.rho_s * dt_new * (Para.err_tol_AM3 / rel_err) ^ (1/Para.p);


end