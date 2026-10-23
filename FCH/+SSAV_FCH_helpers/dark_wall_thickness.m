function t = dark_wall_thickness(u, Grid)
% Mean thickness of the dark phase (u < 0) of a periodic 2D cell:
% t = 2 * (dark area) / (length of the u = 0 boundary in the cell).
% For a honeycomb or a stripe the dark phase is walls bounded by two
% boundary lines, so t is the wall (bilayer) thickness; for light spots in
% a thick dark matrix t is much larger.

d = Grid.L(1);
Ly = Grid.L(2);
x = Grid.x{1}';
y = Grid.x{2};
% pad periodically so contours crossing the cell edge are closed
up = repmat(u, 3, 3);
xp = [x - d, x, x + d];
yp = [y - Ly; y; y + Ly];
C = contourc(xp, yp, up, [0 0]);
Lb = 0;
k = 1;
while k < size(C, 2)
    n = C(2, k);
    pts = C(:, k+1:k+n);
    k = k + n + 1;
    mid = (pts(:, 1:end-1) + pts(:, 2:end)) / 2;   % count segments whose midpoint is in the cell
    inside = mid(1, :) >= x(1) - Grid.d(1)/2 & mid(1, :) < x(end) + Grid.d(1)/2 & ...
             mid(2, :) >= y(1) - Grid.d(2)/2 & mid(2, :) < y(end) + Grid.d(2)/2;
    seg = sqrt(sum(diff(pts, 1, 2).^2, 1));
    Lb = Lb + sum(seg(inside));
end
t = 2 * mean(u(:) < 0) * d * Ly / Lb;
end
