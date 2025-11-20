function [a, b, c] = compute_dt_coeffs(dt, dt_new)

a = 1 / dt_new + 1 / (dt_new + dt);

b = - 1 / dt_new - 1 / dt;

c = 1 / dt - 1 / (dt_new + dt);

end