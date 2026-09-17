tau_vals = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
rate_tol_u = zeros(1,6);
for i = 1:6
    tau_1 = tau_vals(i);
    SSAV_convergence_study;
    u_err(i) = sqrt(sum((Results.uu{end} - u_exact).^2, 'all') ...
        / prod(Grid.N));
end
rate_tol_u = zeros(1,5);
for i = 1:5
    rate_tol_u(i) = log(u_err(i+1)/u_err(i))/log(tau_vals(i+1)/tau_vals(i));
end