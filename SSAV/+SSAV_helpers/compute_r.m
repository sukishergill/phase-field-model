function [r, r_hat] = compute_r(u_snew, u_snew_fft, u_curr, u_prev, w_curr, w_prev, ... 
    H_snew, H2, ~, a, b, c, Grid, Para, G, t, dt)

% r = -(b*w + c*w_old)/a + 0.5 * sum(sum((H_snew .* (b*u + c*u_old)/a)))...
%     * Grid.dx*Grid.dy;
r = -(b*w_curr + c*w_prev)/a + 0.5 * sum((H_snew .* ...
    (b*u_curr + c*u_prev)/a), 'all') * prod(Grid.d);

% rH = fftn(r .* H_snew);         num_fft = num_fft + 1;    

forcing = SSAV_helpers.compute_forcing(Para.epsilon, Grid.xx, Grid.yy, t);

r_tilde = -(b*u_curr + c*u_prev) - Para.S*(ifftn(G .* u_snew_fft, 'symmetric')) + ...
    r*ifftn(G .* H2, 'symmetric');% + forcing;

m_est = (Para.alpha*Para.epsilon^2)/(a*prod(Grid.L)) * ...
    sum(r_tilde, 'all') * prod(Grid.d);

% m_est = Para.OK /(a*prod(Grid.L)) * ...
%     sum(r_tilde, 'all') * prod(Grid.d);

r_hat = r_tilde + m_est;

end