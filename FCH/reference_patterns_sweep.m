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
%   ref_parts      {'stripes', 'spots'}: which families to compute
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

copts = struct('even', true, 'max_steps', ref_max_steps, 'verbose', false, 'm_range', [-0.95 0.95]);
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
        'half_shift', @(v) (v + circshift(v, [Ny/2, Nx/2])) / 2, ...   % spot at origin <-> corners
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

save(ref_out, 'stripes', 'n_stripes', 'spots', 'spot_d', 'Nx_cells', 'Para', 'copts', '-v7.3');
fprintf('Saved %s\n', ref_out);
