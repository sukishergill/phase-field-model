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
% Results are saved to reference_patterns.mat.

Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;
ep = Para.epsilon;
L = 4*pi;
dx = 2*pi/64;                     % grid spacing of the 128^2, L = 4*pi sweep
layer = @(s) (1 - tanh(s / (sqrt(2)*ep))) / 2;   % ~1 inside (s < 0), ~0 outside

copts = struct('even', true, 'max_steps', 3000, 'verbose', false, 'm_range', [-0.95 0.95]);
dirs = [1 -1];

%% stripes
n_stripes = 1:8;
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
R = 0.6;                                          % initial spot radius
spots = cell(numel(Nx_cells), 2);
spot_d = zeros(size(Nx_cells));
for a = 1:numel(Nx_cells)
    Nx = Nx_cells(a);
    d = Nx * dx;
    Ny = 2*round(d*sqrt(3) / (2*dx));             % even, same spacing as x
    Ly = Ny * dx;
    spot_d(a) = d;
    G = SSAV_FCH_helpers.generate_Grid([d, Ly], [Nx, Ny], 2);
    % spots at the origin and at the cell corners (+-d/2, +-Ly/2): even in x and y
    r1 = hypot(G.xx, G.yy);
    ddx = G.xx - d/2;   ddx = ddx - d*round(ddx/d);
    ddy = G.yy - Ly/2;  ddy = ddy - Ly*round(ddy/Ly);
    r2 = hypot(ddx, ddy);
    u = 1 - 2*layer(r1 - R) - 2*layer(r2 - R);
    half_shift = @(v) (v + circshift(v, [Ny/2, Nx/2])) / 2;   % spot at origin <-> corners
    for dd = 1:2
        o = copts;  o.direction = dirs(dd);  o.symmetry = half_shift;
        spots{a, dd} = continuation_FCH(u, G, Para, o);
    end
    fprintf('hex spots d = %.2f (%dx%d, aspect %.4f vs sqrt(3) = %.4f): m in [%.3f, %.3f], folds %d (%s | %s)\n', ...
        d, Nx, Ny, Ly/d, sqrt(3), min([spots{a,1}.m spots{a,2}.m]), max([spots{a,1}.m spots{a,2}.m]), ...
        sum(spots{a,1}.fold) + sum(spots{a,2}.fold), spots{a,1}.stop_reason, spots{a,2}.stop_reason);
end

save('reference_patterns.mat', 'stripes', 'n_stripes', 'spots', 'spot_d', 'Nx_cells', 'Para', 'copts', '-v7.3');
fprintf('Saved reference_patterns.mat\n');
