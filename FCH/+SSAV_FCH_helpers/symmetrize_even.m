function [u, asym] = symmetrize_even(u)
% Project u onto functions even in x (and in y, for 2D) about the origin of
% the generate_Grid grid, x_j = L*j/N - L/2. Reflection x -> -x maps index j
% to N-j, with j = N (x = L/2 = -L/2 on the periodic grid) fixed.
%
% Works on 1D columns, 2D arrays, and stacks of 2D arrays (size [Ny Nx K],
% each page symmetrized). asym = max|u - even part| / max|u - mean(u)|,
% the relative asymmetry that was removed (per page for stacks).

[ny, nx, nk] = size(u);
ix = [nx-1:-1:1, nx];
iy = [ny-1:-1:1, ny];

if nx == 1
    ue = (u + u(iy, :, :)) / 2;
else
    ue = (u + u(:, ix, :) + u(iy, :, :) + u(iy, ix, :)) / 4;
end

asym = zeros(nk, 1);
for k = 1:nk
    d = u(:, :, k) - ue(:, :, k);
    scale = max(abs(u(:, :, k) - mean(u(:, :, k), 'all')), [], 'all');
    asym(k) = max(abs(d), [], 'all') / max(scale, realmin);
end
u = ue;
end
