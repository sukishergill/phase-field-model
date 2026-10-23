function ratio = light_domain_facet_ratio(u, Grid)
% Shape of the light domain centred at the origin of a hex cell (width d,
% height ~ d*sqrt(3), light domains at the origin and the cell corners):
% ratio = r(30 deg) / r(0 deg), the distance from the centre to the u = 0
% boundary towards the triple junction (30 deg) over that towards the
% nearest neighbour (0 deg). A round spot gives ~1; a hexagonal cell with
% straight walls (honeycomb) gives 1/cos(30 deg) ~ 1.155. NaN if the
% centre is not light (u(0,0) <= 0) or no boundary is found.

d = Grid.L(1);
x = Grid.x{1};  y = Grid.x{2};
% periodic extension for interpolation across the cell edge
xe = [x - d; x; x + d];
ye = [y - Grid.L(2); y; y + Grid.L(2)];
ue = repmat(u, 3, 3);
F = griddedInterpolant({ye, xe}, ue, 'spline');
if F(0, 0) <= 0
    ratio = NaN;
    return
end
r = linspace(0, d/sqrt(3), 400);              % out to the triple junction distance
edge = @(th) boundary_distance(F, r, th);
ratio = edge(pi/6) / edge(0);
end

function rb = boundary_distance(F, r, th)
vals = F(r*sin(th), r*cos(th));
k = find(vals <= 0, 1);
if isempty(k) || k == 1
    rb = NaN;
    return
end
% linear interpolation of the zero crossing
rb = r(k-1) + (r(k) - r(k-1)) * vals(k-1) / (vals(k-1) - vals(k));
end
