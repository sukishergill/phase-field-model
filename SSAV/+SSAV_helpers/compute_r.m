function [r, r_hat] = compute_r(u_snew_fft, u_curr, u_prev, w_curr, w_prev, ... 
    H_snew, ~, a, b, c, Grid, Para, G)

% r = -(b*w + c*w_old)/a + 0.5 * sum(sum((H_snew .* (b*u + c*u_old)/a)))...
%     * Grid.dx*Grid.dy;
r = -(b*w_curr + c*w_prev)/a + 0.5 * sum((H_snew .* ...
    (b*u_curr + c*u_prev)/a), 'all') * prod(Grid.d);

rH = fftn(r .* H_snew);         

% num_fft_out = num_fft_in + 1;

r_tilde = -(b*u_curr + c*u_prev) - Para.S*(ifftn(G .* u_snew_fft, 'symmetric')) + ...
    ifftn(G .* rH, 'symmetric');

m_est = (Para.alpha*Para.epsilon^2)/(a*prod(Grid.L)) * ...
    sum(r_tilde, 'all') * prod(Grid.d);

r_hat = r_tilde + m_est;

end