% Locate the edges of the m-windows where an "other" (m-sweep) state has
% lower energy than every simple pattern: n identical stripes (reference
% branches) and hexagonal light-spot lattices at their best spacing
% (coarse d scan, then fminbnd, then a stability check). Dark-spot
% lattices are left out: in these ranges they are > 1e-3 above.
% All comparisons use stable states only.
%
% Settings (override before running):
%   edge_m      masses to evaluate
%   d_range_fn  @(m) [d_min, d_max] for the hex spacing search
%   h           hex-cell grid spacing (0.08)

if ~exist('edge_m', 'var')
    edge_m = [0.255:0.005:0.285, 0.355:0.005:0.375, 0.460:0.005:0.520];
end
if ~exist('d_range_fn', 'var'), d_range_fn = @(m) [2.0, 4.5] + (m > 0.4)*[1.5, 2.5]; end
if ~exist('h', 'var'), h = 0.08; end

Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;

E_other = lowest_dE_at_m(load_sweep_branches(), edge_m, true);
Rf = load_reference_patterns();
E_stripe = lowest_dE_at_m(Rf.stripes, edge_m, true);

nm = numel(edge_m);
E_hex = nan(1, nm);
d_hex = nan(1, nm);
dr = arrayfun(@(m) {d_range_fn(m)}, edge_m);
parfor q = 1:nm
    [E_hex(q), d_hex(q)] = best_light_hex(edge_m(q), dr{q}, Para, h);
end

simple = min(E_stripe, E_hex);
gap = E_other - simple;                       % < 0: the other state is lowest
fprintf('   m      other         stripes       best hex (d)          other - best simple\n');
for q = 1:nm
    if isinf(E_other(q)), note = 'no stable other state here';
    elseif gap(q) < 0, note = 'OTHER LOWEST';
    elseif E_hex(q) < E_stripe(q), note = 'hex lowest';
    else, note = 'stripes lowest';
    end
    fprintf('  %.3f  %+.7f  %+.7f  %+.7f (%.2f)  %+.2e  %s\n', edge_m(q), E_other(q), E_stripe(q), ...
        E_hex(q), d_hex(q), gap(q), note);
end

% edges: sign changes of the gap between neighbouring masses (within each
% contiguous run of evaluated m), or where the other state stops existing
fprintf('\nwindow edges:\n');
for q = 1:nm - 1
    if edge_m(q+1) - edge_m(q) > 0.0051, continue, end      % not neighbours
    g1 = gap(q);  g2 = gap(q+1);
    if isfinite(g1) && isfinite(g2) && sign(g1) ~= sign(g2)
        me = edge_m(q) - g1*(edge_m(q+1) - edge_m(q))/(g2 - g1);
        fprintf('  gap changes sign at m ~ %.4f\n', me);
    elseif xor(isfinite(g1), isfinite(g2))
        fprintf('  other state starts/ends between m = %.3f and %.3f (branch fold)\n', edge_m(q), edge_m(q+1));
    end
end
save('window_edges.mat', 'edge_m', 'E_other', 'E_stripe', 'E_hex', 'd_hex', 'gap', 'h');

function [E, d_best] = best_light_hex(m, d_range, Para, h)
% lowest stable light-spot hex lattice energy (dE) at mass m over spacing d
dgrid = d_range(1):0.1:d_range(2);
e = inf(size(dgrid));
st = false(size(dgrid));
for j = 1:numel(dgrid)
    [e(j), st(j)] = hex_dE(m, dgrid(j), Para, h, true);
end
e(~st) = inf;
[E, j] = min(e);
d_best = dgrid(j);
if ~isfinite(E), d_best = NaN; return, end
% refine the spacing around the best stable grid point
lo = max(d_range(1), d_best - 0.1);
hi = min(d_range(2), d_best + 0.1);
d_opt = fminbnd(@(d) hex_dE(m, d, Para, h, false), lo, hi, optimset('TolX', 1e-3));
[e_chk, stab] = hex_dE(m, d_opt, Para, h, true);   % keep the refinement only if stable
if stab && e_chk < E
    E = e_chk;
    d_best = d_opt;
end
end

function [dE, stable] = hex_dE(m, d, Para, h, check_stability)
% light-spot hex lattice at spacing d and mass m on one exact hex cell
ep = Para.epsilon;
layer = @(s) (1 - tanh(s / (sqrt(2)*ep))) / 2;
Ly = d*sqrt(3);
G = SSAV_FCH_helpers.generate_Grid([d, Ly], [2*ceil(d/(2*h)), 2*ceil(Ly/(2*h))], 2);
r1 = hypot(G.xx, G.yy);
ddx = G.xx - d/2;   ddx = ddx - d*round(ddx/d);
ddy = G.yy - Ly/2;  ddy = ddy - Ly*round(ddy/Ly);
r2 = hypot(ddx, ddy);
R = sqrt((1 + m)/2 * d * Ly / (2*pi));
u = -1 + 2*(layer(r1 - R) + layer(r2 - R));
u = u - mean(u(:)) + m;
half_shift = @(v) (v + circshift(v, size(v)/2)) / 2;
[u, info] = newton_polish_FCH(u, G, Para, struct('even', true, 'symmetry', half_shift, ...
    'verbose', false, 'maxit', 3000));
dE = info.dE;
vx = mean(std(u, 0, 2));  vy = mean(std(u, 0, 1));
if ~info.converged || info.collapsed || min(vx, vy) < 0.05*max(vx, vy)
    dE = inf;                                  % not a converged spot lattice
end
stable = true;
if check_stability && isfinite(dE)
    [~, ~, si] = stability_FCH(u, G, Para, struct('even_only', true, 'verbose', false));
    stable = si.stable;
end
end
