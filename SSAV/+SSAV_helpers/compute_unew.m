function [u_new, w_new, num_fft] = compute_unew(u_curr, u_curr_fft, ...
    u_prev, u_prev_fft, w_curr, w_prev,...
    dt, dt_new, Para, Grid, G, P1, num_fft, t)

u_snew = SSAV_helpers.compute_u_snew(u_curr, u_prev, dt_new/dt);      % extrapolation of u
    
[F, f] = SSAV_helpers.compute_F(u_snew, Para.beta);

int_F = sum(F(:)) * prod(Grid.d);

% w(u^{*,n+1})
% w_snew = SSAV_helpers.compute_w(int_F, Para.B);
w_snew = (1 + dt_new/dt)*w_curr - (dt_new/dt)*w_prev;

% H^{*,n+1}
H_snew = SSAV_helpers.compute_H(f, w_snew);

[a, b, c] =  SSAV_helpers.compute_dt_coeffs(dt, dt_new);

H2 = fftn(H_snew);          num_fft = num_fft + 1;

[r, r_hat] = SSAV_helpers.compute_r(u_snew, 2*u_curr_fft - u_prev_fft, u_curr,...
    u_prev, w_curr, w_prev, H_snew, H2, dt_new, a, b, c, Grid, Para, G, t + dt_new, dt_new);

% define P for BDF2 
P = a + P1;

r_hat = fftn(r_hat);        num_fft = num_fft + 1;
psi_r = P .\ r_hat;
psi_r = ifftn(psi_r, 'symmetric');

% H2 = fftn(H_snew);          num_fft = num_fft + 1;
psi_H = P .\ (G.*H2);
psi_H = ifftn(psi_H, 'symmetric');

innprod_Hu = SSAV_helpers.compute_ip(H_snew, psi_r, psi_H, Grid);

w_new = 0.5*innprod_Hu + r;      

u_new = 0.5*innprod_Hu * psi_H + psi_r;

end