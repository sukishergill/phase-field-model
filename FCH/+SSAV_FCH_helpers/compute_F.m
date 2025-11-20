function [F, F1, F2] = compute_F(u, tau)

F = 0.5*(u + tau/3 + 1).^2 .* (0.5*(u + tau/3 - 1).^2 + ...
    (2*tau/3) * (u + tau/3 - 2));

F1 = u.^3 - (1 + tau^2/3)*u - 2*tau/3 * (1 + tau/3)*(1 - tau/3);

F2 = 3*u.^2 - (1 + tau^2/3);

end