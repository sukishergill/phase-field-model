function x_tt = compute_deriv2(x0, x1, x2, dt0, dt1)
% computes the second derivative using backwards difference

a = 2 / (dt1 * (dt1 + dt0));

b = -2 / (dt1 * dt0);

c = -a - b;


x_tt = a * x2 + b * x1 + c * x0;

end