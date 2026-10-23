% Reference ("simple pattern") branches in m, for comparison with the
% m-sweep branches (continue_m_sweep.m). Energies are per unit area, so
% they compare directly.
%
%   stripes   n identical stripes on L = 4*pi, n = 1..8, computed on one
%             period (length 4*pi/n, one stripe, grid spacing ~ the
%             sweep's): a y-independent 2D state has the 1D energy per
%             area. On the full 4*pi line the branches wander off into
%             unequal-stripe states (mass exchange between stripes), so
%             the one-period cell is what keeps the stripes identical.
%             Stability: 1D, one-period, even perturbations only.
%   hex spots ideal triangular lattice of u ~ -1 spots, computed on one
%             rectangular cell (width d, height ~ d*sqrt(3), two spots)
%             with the sweep's grid spacing; the two spots are kept
%             identical by averaging with the half-cell shift. A hexagonal
%             lattice cannot tile the square 4*pi box exactly; this is the
%             ideal-lattice energy per area. Stability: cell-periodic even
%             perturbations, including ones that make the spots unequal.
%
% Settings (set any of these before running the script to override):
%   ref_parts      {'stripes', 'spots'}: which families to compute; add
%                  'spots_light' for hex lattices of light spots (u ~ +1
%                  in u ~ -1) on exact hex cells, spacings ref_light_d,
%                  started at ref_light_m(d) (saved as spots, with labels);
%                  add 'honeycomb' for light lattices with d >= 4 at
%                  spacings ref_honey_d, each started at every mass in
%                  ref_honey_m0 (coexisting states), saved as honeycomb
%   ref_m_range    continuation m range ([-0.95 0.95])
%   ref_ds_max     continuation ds_max ([] = continuation_FCH default)
%   ref_max_steps  continuation max_steps per branch (3000)
%   ref_out        output file ('reference_patterns.mat'); any
%                  reference_patterns*.mat is merged by load_reference_patterns
% Hex branches run in parallel (parfor) when a pool is available.

Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;
ep = Para.epsilon;
L = 4*pi;
dx = 2*pi/64;                     % grid spacing of the 128^2, L = 4*pi sweep
layer = @(s) (1 - tanh(s / (sqrt(2)*ep))) / 2;   % ~1 inside (s < 0), ~0 outside

if ~exist('ref_parts', 'var'), ref_parts = {'stripes', 'spots'}; end
if ~exist('ref_ds_max', 'var'), ref_ds_max = []; end
if ~exist('ref_max_steps', 'var'), ref_max_steps = 3000; end
if ~exist('ref_out', 'var'), ref_out = 'reference_patterns.mat'; end
if ~exist('ref_light_d', 'var'), ref_light_d = [1.3 1.45 1.6 1.8 2.0 2.3 2.6 3.0 3.35]; end
if ~exist('ref_light_m', 'var'), ref_light_m = @(d) 0.118 + (d > 1.85)*(0.35 - 0.118); end
if ~exist('ref_honey_d', 'var'), ref_honey_d = 4.0:0.5:7.0; end
if ~exist('ref_honey_m0', 'var'), ref_honey_m0 = [0.40 0.50]; end
if ~exist('ref_m_range', 'var'), ref_m_range = [-0.95 0.95]; end

copts = struct('even', true, 'max_steps', ref_max_steps, 'verbose', false, 'm_range', ref_m_range);
if ~isempty(ref_ds_max), copts.ds_max = ref_ds_max; end
dirs = [1 -1];
fprintf('reference patterns: %s, ds_max = %s, max_steps = %d -> %s\n', strjoin(ref_parts, ' + '), ...
    mat2str(ref_ds_max), ref_max_steps, ref_out);

%% stripes
n_stripes = 1:8;
if ~ismember('stripes', ref_parts), n_stripes = []; end
w = 0.5;                                          % initial stripe width
stripes = cell(numel(n_stripes), 2);
for a = 1:numel(n_stripes)
    n = n_stripes(a);
    p = L / n;                                    % period = cell length
    Nc = 2*round(p / (2*dx));                     % even point count, spacing ~ dx
    G1 = SSAV_FCH_helpers.generate_Grid(p, Nc, 1);
    u = 1 - 2*layer(abs(G1.x) - w/2);             % one stripe centred at x = 0
    for d = 1:2
        o = copts;  o.direction = dirs(d);
        stripes{a, d} = continuation_FCH(u, G1, Para, o);
    end
    fprintf('stripes n = %d: m in [%.3f, %.3f], %d + %d points, folds %d (%s | %s)\n', n, ...
        min([stripes{a,1}.m stripes{a,2}.m]), max([stripes{a,1}.m stripes{a,2}.m]), ...
        numel(stripes{a,1}.m), numel(stripes{a,2}.m), sum(stripes{a,1}.fold) + sum(stripes{a,2}.fold), ...
        stripes{a,1}.stop_reason, stripes{a,2}.stop_reason);
end

%% hexagonal spot lattices
Nx_cells = [16 20 24 32];                         % cell widths d = Nx*dx: 1.57, 1.96, 2.36, 3.14
if ~ismember('spots', ref_parts), Nx_cells = []; end
R = 0.6;                                          % initial spot radius
ncell = numel(Nx_cells);
spot_d = Nx_cells * dx;
cells = cell(ncell, 1);                           % grid, initial state, symmetry per cell
for a = 1:ncell
    Nx = Nx_cells(a);
    d = spot_d(a);
    Ny = 2*round(d*sqrt(3) / (2*dx));             % even, same spacing as x
    Ly = Ny * dx;
    G = SSAV_FCH_helpers.generate_Grid([d, Ly], [Nx, Ny], 2);
    % spots at the origin and at the cell corners (+-d/2, +-Ly/2): even in x and y
    r1 = hypot(G.xx, G.yy);
    ddx = G.xx - d/2;   ddx = ddx - d*round(ddx/d);
    ddy = G.yy - Ly/2;  ddy = ddy - Ly*round(ddy/Ly);
    r2 = hypot(ddx, ddy);
    cells{a} = struct('G', G, 'u', 1 - 2*layer(r1 - R) - 2*layer(r2 - R), ...
        'half_shift', @(v) (v + circshift(v, size(v)/2)) / 2, ...   % spot at origin <-> corners (any resolution)
        'Nx', Nx, 'Ny', Ny, 'aspect', Ly/d);
end
spots = cell(ncell, 2);
parfor job = 1:2*ncell
    [a, dd] = ind2sub([ncell, 2], job);
    o = copts;  o.direction = dirs(dd);  o.symmetry = cells{a}.half_shift;
    spots{job} = continuation_FCH(cells{a}.u, cells{a}.G, Para, o);
end
spots = reshape(spots, ncell, 2);
for a = 1:ncell
    fprintf('hex spots d = %.2f (%dx%d, aspect %.4f vs sqrt(3) = %.4f): m in [%.3f, %.3f], %d + %d points, folds %d (%s | %s)\n', ...
        spot_d(a), cells{a}.Nx, cells{a}.Ny, cells{a}.aspect, sqrt(3), min([spots{a,1}.m spots{a,2}.m]), ...
        max([spots{a,1}.m spots{a,2}.m]), numel(spots{a,1}.m), numel(spots{a,2}.m), ...
        sum(spots{a,1}.fold) + sum(spots{a,2}.fold), spots{a,1}.stop_reason, spots{a,2}.stop_reason);
end

spot_label = arrayfun(@(d) sprintf('d=%.2f', d), spot_d, 'UniformOutput', false);

%% hexagonal lattices of light spots (exact hex cells, any d)
if ismember('spots_light', ref_parts)
    h = 0.08;                                     % target grid spacing
    nl = numel(ref_light_d);
    lcells = cell(nl, 1);
    for a = 1:nl
        d = ref_light_d(a);
        Ly = d*sqrt(3);
        G = SSAV_FCH_helpers.generate_Grid([d, Ly], [2*ceil(d/(2*h)), 2*ceil(Ly/(2*h))], 2);
        r1 = hypot(G.xx, G.yy);
        ddx = G.xx - d/2;   ddx = ddx - d*round(ddx/d);
        ddy = G.yy - Ly/2;  ddy = ddy - Ly*round(ddy/Ly);
        r2 = hypot(ddx, ddy);
        m0 = ref_light_m(d);
        R = sqrt((1 + m0)/2 * d * Ly / (2*pi));   % light-phase area fraction (1+m)/2
        u = -1 + 2*(layer(r1 - R) + layer(r2 - R));
        lcells{a} = struct('G', G, 'u', u - mean(u(:)) + m0, ...
            'half_shift', @(v) (v + circshift(v, size(v)/2)) / 2);
    end
    light = cell(nl, 2);
    parfor job = 1:2*nl
        [a, dd] = ind2sub([nl, 2], job);
        o = copts;  o.direction = dirs(dd);  o.symmetry = lcells{a}.half_shift;
        light{job} = continuation_FCH(lcells{a}.u, lcells{a}.G, Para, o);
    end
    light = reshape(light, nl, 2);
    for a = 1:nl
        fprintf('hex light spots d = %.2f (start m %.3f): m in [%.3f, %.3f], %d + %d points, folds %d (%s | %s)\n', ...
            ref_light_d(a), ref_light_m(ref_light_d(a)), min([light{a,1}.m light{a,2}.m]), ...
            max([light{a,1}.m light{a,2}.m]), numel(light{a,1}.m), numel(light{a,2}.m), ...
            sum(light{a,1}.fold) + sum(light{a,2}.fold), light{a,1}.stop_reason, light{a,2}.stop_reason);
    end
    spots = [spots; light];
    spot_d = [spot_d, ref_light_d];
    spot_label = [spot_label, arrayfun(@(d) sprintf('d=%.2f light', d), ref_light_d, 'UniformOutput', false)];
end

%% honeycombs: light lattices with d >= 4, several starting masses each
honeycomb = cell(0, 2);
honey_d = [];
honey_label = {};
if ismember('honeycomb', ref_parts)
    [Dg, Mg] = meshgrid(ref_honey_d, ref_honey_m0);
    starts = [Dg(:), Mg(:)];                      % one row per (d, m0)
    ns = size(starts, 1);
    hstate = cell(ns, 1);
    parfor a = 1:ns                               % starting states, polished at (d, m0)
        [u0, G0] = SSAV_FCH_helpers.hex_lattice_state(starts(a, 2), starts(a, 1), 'light', Para, 0.098);
        hstate{a} = struct('u', u0, 'G', G0);
    end
    honeycomb = cell(ns, 2);
    half_shift = @(v) (v + circshift(v, size(v)/2)) / 2;
    parfor job = 1:2*ns
        [a, dd] = ind2sub([ns, 2], job);
        o = copts;  o.direction = dirs(dd);  o.symmetry = half_shift;  o.polish_first = false;
        honeycomb{job} = continuation_FCH(hstate{a}.u, hstate{a}.G, Para, o);
    end
    honeycomb = reshape(honeycomb, ns, 2);
    honey_d = starts(:, 1)';
    honey_label = arrayfun(@(a) sprintf('d=%.2f (from m %.2f)', starts(a, 1), starts(a, 2)), 1:ns, 'UniformOutput', false);
    for a = 1:ns
        fprintf('honeycomb d = %.2f from m = %.2f: m in [%.3f, %.3f], %d + %d points, folds %d (%s | %s)\n', ...
            starts(a, 1), starts(a, 2), min([honeycomb{a,1}.m honeycomb{a,2}.m]), max([honeycomb{a,1}.m honeycomb{a,2}.m]), ...
            numel(honeycomb{a,1}.m), numel(honeycomb{a,2}.m), sum(honeycomb{a,1}.fold) + sum(honeycomb{a,2}.fold), ...
            honeycomb{a,1}.stop_reason, honeycomb{a,2}.stop_reason);
    end
end

save(ref_out, 'stripes', 'n_stripes', 'spots', 'spot_d', 'spot_label', 'honeycomb', 'honey_d', 'honey_label', ...
    'Nx_cells', 'Para', 'copts', '-v7.3');
fprintf('Saved %s\n', ref_out);
