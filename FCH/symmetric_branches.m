function [hex, stripe] = symmetric_branches(d, Para, opts)
% Stripe and hexagonal branches with spacing d, each with its symmetry
% enforced, for 0 <= m < 1 (E(-u) = E(u), so m < 0 is the mirror image).
%
%   hex     hexagonal states with nearest-neighbour spacing d, on one exact
%           hex cell (d x d*sqrt(3)), continued by continuation_FCH
%           restricted to the full hexagonal (D6) symmetry
%           (SSAV_FCH_helpers.hex_symmetry), both directions down to m = 0,
%           through folds and through the uniform state. Two runs: one
%           started from dark spots (u ~ -1 in u ~ +1) and one from light
%           spots (u ~ +1 in u ~ -1) at small m, since the two families
%           need not lie on one closed branch for every d; hex.run says
%           which run a point is from. Each point is labelled by the sign
%           of the pattern at the lattice sites: 'spots' where u there is
%           below the mean (dark spots), 'network' where it is above (light
%           domains on the lattice sites, dark walls/triangles on the
%           honeycomb sites).
%   stripe  stripes of period d as a 1D problem (u independent of y), on
%           the same x grid, both directions down to m = 0.
%
% stable: no negative eigenvalue of the Hessian on even, mean-zero
% perturbations of the 2D hex cell (n_neg_full for hex; recomputed on the
% cell for the stripes, so it includes bending them or breaking them up).
%
% opts (all optional)
%   h          target grid spacing                          (0.08)
%   m0_hex     starting masses to try for the dark spots, in order, until
%              hex_lattice_state returns a lattice     ([0.57 0.45 0.7 0.35])
%   m0_light   starting masses to try for the light spots   ([0.15 0.3 0.05])
%   m0_stripe  starting mass for the stripes                (0.25)
%   ds_max     continuation step limit                      (0.05)
%   store_u    keep the profiles (single precision)        (true)
%
% Each output has per-point arrays m, E, amp, stable, fold, cross_uniform,
% (hex) n_neg, n_neg_full, part ('spots' / 'network'), (stripe) n_neg,
% n_neg_cell; u (cell, if store_u); plus d, i_start, m0, stop_reason, and
% ok (false if no branch could be started).

if nargin < 3, opts = struct(); end
def = struct('h', 0.08, 'm0_hex', [0.57 0.45 0.7 0.35], 'm0_light', [0.15 0.3 0.05], 'm0_stripe', 0.25, ...
    'ds_max', 0.05, 'store_u', true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end
end
copts = struct('even', true, 'stop_on_uniform', false, 'ds_max', opts.ds_max, ...
    'max_steps', 3000, 'm_range', [0 0.999], 'verbose', false);

% ---- hex: one run from dark spots, one from light spots ----
hex = struct('ok', false, 'd', d);
Grid = [];
runs = {};
starts = {'dark', opts.m0_hex; 'light', opts.m0_light};
for r = 1:2
    kind = '';
    for m0 = starts{r, 2}
        [u0, Grid, ~, kind] = SSAV_FCH_helpers.hex_lattice_state(m0, d, starts{r, 1}, Para, opts.h);
        if strcmp(kind, 'lattice'), break, end
    end
    if ~strcmp(kind, 'lattice'), continue, end
    S = SSAV_FCH_helpers.hex_symmetry(Grid);
    o = copts;  o.symmetry = S;  o.restrict = true;
    o.direction = 1;   b_up = continuation_FCH(u0, Grid, Para, o);
    o.direction = -1;  b_down = continuation_FCH(u0, Grid, Para, o);
    b = join_branches(b_down, b_up, true);
    b.run = r * ones(1, numel(b.m));
    b.m0 = m0;
    runs{end+1} = b; %#ok<AGROW>
end
if ~isempty(runs)
    f = {'m', 'E', 'dE', 'amp', 'n_neg', 'n_neg_full', 'fold', 'cross_uniform', 'u', 'run'};
    for i = 1:numel(f)
        hex.(f{i}) = cellfun(@(b) b.(f{i}), runs, 'UniformOutput', false);
        hex.(f{i}) = [hex.(f{i}){:}];
    end
    hex.i_start = cellfun(@(b) b.i_start, runs);
    hex.m0 = cellfun(@(b) b.m0, runs);
    hex.stop_reason = cellfun(@(b) b.stop_reason, runs, 'UniformOutput', false);
    hex.stable = hex.n_neg_full == 0;
    % sign of the pattern at the lattice site x = y = 0 (grid point N/2)
    ic = Grid.N(1)/2;  jc = Grid.N(2)/2;
    centre = cellfun(@(u) double(u(jc, ic)), hex.u);
    hex.part = repmat({'spots'}, 1, numel(hex.m));
    hex.part(centre > hex.m) = {'network'};
    hex.ok = true;
    if ~opts.store_u, hex = rmfield(hex, 'u'); end
end

% ---- stripes, on the same x grid (a fresh cell if the hex start failed) ----
if isempty(Grid)
    Grid = SSAV_FCH_helpers.generate_Grid([d, d*sqrt(3)], ...
        [2*ceil(d/(2*opts.h)), 2*ceil(d*sqrt(3)/(2*opts.h))], 2);
end
G1 = SSAV_FCH_helpers.generate_Grid(Grid.L(1), Grid.N(1), 1);
layer = @(s) (1 - tanh(s / (sqrt(2)*Para.epsilon))) / 2;
m0 = opts.m0_stripe;
w = (1 - m0) / 2 * d;
u1 = 1 - 2*layer(abs(G1.x) - w/2);
u1 = u1 - mean(u1) + m0;
o = copts;
o.direction = 1;   s_up = continuation_FCH(u1, G1, Para, o);
o.direction = -1;  s_down = continuation_FCH(u1, G1, Para, o);
stripe = join_branches(s_down, s_up, true);
Z = SSAV_FCH_helpers.even_basis([Grid.N(2), Grid.N(1)]);
stripe.n_neg_cell = zeros(size(stripe.m));
for k = 1:numel(stripe.m)
    stripe.n_neg_cell(k) = n_neg_even(repmat(double(stripe.u{k}(:))', Grid.N(2), 1), Z, Grid, Para);
end
stripe.stable = stripe.n_neg_cell == 0;
stripe.m0 = m0;
stripe.d = d;
stripe.ok = stripe.amp(stripe.i_start) > 0.05;
if ~opts.store_u, stripe = rmfield(stripe, 'u'); end

end

function b = join_branches(b_down, b_up, store_u)
% reversed downward branch followed by the upward one (start point once)
f = {'m', 'E', 'dE', 'amp', 'n_neg', 'n_neg_full', 'fold', 'cross_uniform'};
if store_u, f{end+1} = 'u'; end
for i = 1:numel(f)
    a = b_down.(f{i});
    c = b_up.(f{i});
    b.(f{i}) = [fliplr(a), c(2:end)];
end
if store_u
    b.u = cellfun(@single, b.u, 'UniformOutput', false);
end
b.i_start = numel(b_down.m);
b.stop_reason = {b_down.stop_reason, b_up.stop_reason};
end

function n_neg = n_neg_even(u, Z, Grid, Para)
% number of negative eigenvalues of the Hessian on the even, mean-zero
% perturbations of the cell (columns of Z)
sz = size(u);
n = numel(u);
ne = size(Z, 2);
JZ = zeros(n, ne);
for c0 = 1:64:ne
    cols = c0:min(c0 + 63, ne);
    JZ(:, cols) = reshape(SSAV_FCH_helpers.steady_Jv(u, reshape(Z(:, cols), [sz, numel(cols)]), ...
        Para, Grid), n, []);
end
Je = Z' * JZ;
n_neg = sum(eig((Je + Je') / 2) < -1e-6);
end
