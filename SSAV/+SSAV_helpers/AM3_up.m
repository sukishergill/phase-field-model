function u_AM3 = AM3_up(u, gamma, dt_new, ...
        F_tilde_old, F_tilde, F_tilde_new)

u_AM3 = u + (dt_new / (6*(1 + gamma))) * ((3 + 2*gamma)*F_tilde_new + ...
    (3 + gamma)*(1 + gamma)*F_tilde - (gamma^2)*F_tilde_old);

end