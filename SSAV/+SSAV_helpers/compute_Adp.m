function A_dp = compute_Adp(rel_err, dt, Para)

A_dp = Para.rho_s * dt * (Para.err_tol / rel_err) ^ (1/Para.p);

end