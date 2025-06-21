function [F, f] = compute_F(u, beta)
    
f = u.^3 - beta*u;

F = 0.25 * (u.^2 - beta).^2;

end