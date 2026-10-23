function [u, G, info, kind, ratio] = hex_lattice_state(m, d, polarity, Para, h)
% Polished hexagonal lattice at mass m and spacing d on one exact hex cell
% (width d, height d*sqrt(3), domains at the origin and the corners, kept
% identical by the half-cell shift symmetry), grid spacing ~h.
%
%   polarity  'light' (u ~ +1 domains in u ~ -1) or 'dark'
%
%   kind   'lattice' if it stayed a hexagonal lattice; 'stripes' /
%          'uniform' / 'failed' if it relaxed into those or did not
%          converge. Light lattices turn continuously into honeycombs of
%          membranes as d grows; the spots / honeycomb split is a
%          convention made by the caller (by d: honeycomb for d >= 4).
%   ratio  light_domain_facet_ratio (diagnostic only; NaN for dark centres)

if nargin < 5, h = 0.08; end
ep = Para.epsilon;
layer = @(s) (1 - tanh(s / (sqrt(2)*ep))) / 2;
Ly = d*sqrt(3);
G = SSAV_FCH_helpers.generate_Grid([d, Ly], [2*ceil(d/(2*h)), 2*ceil(Ly/(2*h))], 2);
r1 = hypot(G.xx, G.yy);
ddx = G.xx - d/2;   ddx = ddx - d*round(ddx/d);
ddy = G.yy - Ly/2;  ddy = ddy - Ly*round(ddy/Ly);
r2 = hypot(ddx, ddy);
if strcmp(polarity, 'light')
    R = sqrt((1 + m)/2 * d * Ly / (2*pi));
    u = -1 + 2*(layer(r1 - R) + layer(r2 - R));
else
    R = sqrt((1 - m)/2 * d * Ly / (2*pi));
    u = 1 - 2*(layer(r1 - R) + layer(r2 - R));
end
u = u - mean(u(:)) + m;
half_shift = @(v) (v + circshift(v, size(v)/2)) / 2;
[u, info] = newton_polish_FCH(u, G, Para, struct('even', true, 'symmetry', half_shift, ...
    'verbose', false, 'maxit', 3000));

ratio = NaN;
vx = mean(std(u, 0, 2));  vy = mean(std(u, 0, 1));
if ~info.converged
    kind = 'failed';
elseif info.collapsed
    kind = 'uniform';
elseif min(vx, vy) < 0.05*max(vx, vy)
    kind = 'stripes';
else
    ratio = SSAV_FCH_helpers.light_domain_facet_ratio(u, G);
    kind = 'lattice';
end
end
