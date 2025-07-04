function F_tilde = compute_F_tilde(u, u_fft, f, Grid, Para)

f_fft = fft2(f);

F_tilde = -Para.epsilon^2 * ifft2(Grid.k4 .* u_fft, 'symmetric') - ...
    ifft2(Grid.k2 .* f_fft, 'symmetric') + ...
    Para.alpha * Para.epsilon * (u - Para.m);

end