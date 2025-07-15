function [F_tilde, num_fft_out] = compute_F_tilde(u, u_fft, f, Grid,...
    Para, num_fft_in)

f_fft = fft2(f);

F_tilde = -Para.epsilon^2 * ifft2(Grid.k4 .* u_fft, 'symmetric') - ...
    ifft2(Grid.k2 .* f_fft, 'symmetric') + ...
    Para.alpha * Para.epsilon * (u - Para.m);

num_fft_out = num_fft_in + 1;
end