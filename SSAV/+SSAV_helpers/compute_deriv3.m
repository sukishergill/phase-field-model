function x_ttt = compute_deriv3(x, dt0, dt1, dt2)
% computes the third derivative using backwards difference

s1 = dt0;       s2 = dt0 + dt1;         s3 = s2 + dt2;

a = -6 / (s3 * (s1 - s3) * (s2 - s3));

b = 6 / (s2 * (s1 - s2) * (s2 - s3));

c = -6 / (s1 * (s1 - s2) * (s1 - s3));

d = -a - b - c;

x_ttt = a * x{4} + b * x{3} + c * x{2} + d * x{1};

end